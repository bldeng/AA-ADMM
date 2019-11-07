// Original work Copyright (c) 2017, University of Minnesota
// Modified work Copyright 2019, Yue Peng
//
// ADMM-Elastic Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list
//    of conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF MINNESOTA, DULUTH OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Solver.hpp"
#include "MCL/MicroTimer.hpp"
#include <fstream>
#include <unordered_set>
#include <unordered_map>

using namespace admm;
using namespace Eigen;

Solver::Solver() : initialized(false) {
    m_constraints = std::make_shared<ConstraintSet>( ConstraintSet() );
}

void Solver::step(){
    if( m_settings.verbose > 0 ){
        std::cout << "\nSimulating with dt: " <<
                     m_settings.timestep_s << "s..." << std::flush;
    }

    mcl::MicroTimer t;

    // Other const-variable short names and runtime data
    const int dof = static_cast<int>(m_x.rows());
    const int n_nodes = dof/3;
    const double dt = m_settings.timestep_s;
    const int n_energyterms = static_cast<int>(energyterms.size());
    const int n_threads = std::min(n_energyterms, omp_get_max_threads());
    m_runtime = RuntimeData(); // reset

    // Take an explicit step to get predicted node positions
    // with simple forces (e.g. wind).
    for( size_t i=0; i<ext_forces.size(); ++i ){
        ext_forces[i]->project( dt, m_x, m_v, m_masses );
    }

    // Add gravity
    if( std::abs(m_settings.gravity)>0 ){
        for( int i=0; i<n_nodes; ++i ){
            if (positive_pin(i) > 0)// Free points
            {
                m_v[i*3+1] += dt*m_settings.gravity;
            }
        }
    }

    int count = 0;
    for( int i=0; i<n_nodes; ++i ){
        if (positive_pin(i) > 0)// Free points
        {
            m_x_free.segment<3>(3*count) = m_x.segment<3>(i*3) + dt * m_v.segment<3>(i*3);
            count++;
        }
    }

    // Storing some of these variables as members can reduce allocation overhead
    // and speed up run time.

    // Position without elasticity/constraints
    VecX x_bar = m_x_free;
    VecX M_xbar = m_M * x_bar;
    VecX curr_x = x_bar; // Temperorary x used in optimization

    // Initialize ADMM vars
    VecX C_fix = m_C * m_x_pin;
    VecX A_xbar = m_D*x_bar;
    VecX curr_z = solver_W_inv * (A_xbar - C_fix);//m_D*x_bar;
    VecX last_z = curr_z;

    VecX curr_u = VecX::Zero( curr_z.rows() );
    VecX solver_termB = VecX::Zero( m_x_free.rows() );
    VecX gradient = VecX::Zero( curr_z.rows() );

    //Anderson acceleration
    AndersonAcceleration accelerator;
    int z_size = static_cast<int>(curr_z.rows());

    VecX dual_residual_Vec = VecX::Zero( z_size );
    VecX prim_residual_Vec = VecX::Zero( z_size );
    VecX energy_sum = VecX::Zero( n_energyterms );
    double prim_residual = 0.0, comb_residual = 0.0, prev_prim_residual = 1e+20, eps = 1e-20;

    t.reset();
    // Update x for initialization
    solver_termB.noalias() = M_xbar + solver_Dt_Wt_W * ( solver_W * curr_z + C_fix - curr_u );
    m_linsolver->solve( curr_x, solver_termB );

    //Update z
#pragma omp parallel for num_threads(n_threads)
    for( size_t i=0; i<energyterms.size(); ++i ){
        energyterms[i]->update_z( m_D, solver_W_inv, solver_W, curr_x, curr_z, curr_u, C_fix);
    }

    VecX default_z = curr_z, default_x = curr_x, default_u = curr_u;

    accelerator.init(m_settings.Anderson_m, z_size, curr_z);

    m_runtime.initialization_ms += t.elapsed_ms();

    int reset_num = 0;
    // Run a timestep
    int s_i = 0;
    for( ; s_i < m_settings.admm_iters; ++s_i ){

        t.reset();
        if (m_settings.acceleration_type == Settings::ANDERSON)
        {
            //Compute Gradient g(x)
#pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->get_all_gradient( curr_z, gradient );
            }

            curr_u = solver_W_inv * gradient;
        }
        else
        {
            // Update u
#pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->update_u( m_D, solver_W, curr_x, curr_z, curr_u, C_fix);
            }
        }

        m_runtime.local_ms += t.elapsed_ms();

        // Update x
        t.reset();
        solver_termB.noalias() = M_xbar + solver_Dt_Wt_W * ( solver_W * curr_z + C_fix - curr_u );
        m_runtime.inner_iters += m_linsolver->solve( curr_x, solver_termB );
        m_runtime.global_ms += t.elapsed_ms();

        //Compute residual
        t.reset();
        prim_residual_Vec = m_D * curr_x - solver_W *curr_z - C_fix;
        prim_residual = prim_residual_Vec.norm();


        // If reject last accelerating solver
        if ((m_settings.acceleration_type == Settings::ANDERSON) && (prev_prim_residual < prim_residual))
        {
            reset_num = reset_num + 1;
            // Reset curr_u, curr_x, curr_z
            curr_u = default_u;
            curr_x = default_x;
            curr_z = default_z;
            accelerator.replace(curr_z);

            // Update u
            #pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->update_u( m_D, solver_W, curr_x, curr_z, curr_u, C_fix);
            }

            // Update x
            solver_termB.noalias() = M_xbar + solver_Dt_Wt_W * ( solver_W * curr_z + C_fix - curr_u );
            m_runtime.inner_iters += m_linsolver->solve( curr_x, solver_termB );

            //Compute residual
            prim_residual_Vec = m_D * curr_x - solver_W *curr_z - C_fix;
            prim_residual = prim_residual_Vec.norm();
        }

        prev_prim_residual = prim_residual;

        m_runtime.acceleration_ms += t.elapsed_ms();

        // Get combined residual
        last_z = curr_z;

        t.reset();
        //Update z using Anderson acceleration
        if (m_settings.acceleration_type == Settings::ANDERSON)
        {
            default_x = curr_x;
            default_u = curr_u;

            //Update default z
#pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->update_z( m_D, solver_W_inv, solver_W, curr_x, default_z, curr_u, C_fix);
            }

            accelerator.compute(curr_z, default_z);
        }
        else
        {
            //Update z
#pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->update_z( m_D, solver_W_inv, solver_W, curr_x, curr_z, curr_u, C_fix);
            }

        }
        m_runtime.local_ms += t.elapsed_ms();

        // Calc combined residual for drawing figures, actually not need it in application.
        if (m_settings.acceleration_type == Settings::ANDERSON)
        {
            VecX comb_x = curr_x, comb_z = curr_z;
            // Update x for initialization
            solver_termB.noalias() = M_xbar + solver_Dt_Wt_W * ( solver_W * default_z + C_fix - curr_u );
            m_linsolver->solve( comb_x, solver_termB );

            //Update z
        #pragma omp parallel for num_threads(n_threads)
            for( size_t i=0; i<energyterms.size(); ++i ){
                energyterms[i]->update_z( m_D, solver_W_inv, solver_W, comb_x, comb_z, curr_u, C_fix);
            }

            dual_residual_Vec = solver_W * (comb_z - default_z);
            prim_residual_Vec = m_D * comb_x - solver_W *comb_z - C_fix;
            comb_residual = dual_residual_Vec.squaredNorm() + prim_residual_Vec.squaredNorm();

        } else {
            dual_residual_Vec = solver_W * (curr_z - last_z);
            prim_residual_Vec = m_D * curr_x - solver_W *curr_z - C_fix;
            comb_residual = dual_residual_Vec.squaredNorm() + prim_residual_Vec.squaredNorm();
        }


        //Save step time
        step_prim_residual.push_back( prim_residual );
        step_comb_residual.push_back( comb_residual );
        m_runtime.step_time.push_back(m_runtime.local_ms+m_runtime.global_ms+m_runtime.acceleration_ms);

        // End iteration
        if (comb_residual < eps)
        {
            break;
        }
    } // end solver loop

    std::cout << "reset number = " << reset_num << std::endl;

    VecX actual_x = m_S_free * curr_x + m_S_fix * m_x_pin;
    m_v.noalias() = ( actual_x - m_x ) * ( 1.0 / dt );
    m_x = actual_x;

    // Output run time
    if( m_settings.verbose > 0 ){ m_runtime.print(m_settings); }

    save();
} // end timestep iteration

double Solver::Calc_function(const VecX &curr_x)
{
    VecX x_div_x_bar = curr_x - m_x_free;
    double energy = x_div_x_bar.transpose() * m_M * x_div_x_bar;
    return 0.5 * energy;
}

double Solver::Calc_function(MatrixXX &all_s, const VecX &curr_x, const VecX &prev_x, int row)
{
    VecX x_div_x_bar = curr_x - m_x_free;

    if (row > -1) {
        VecX F = curr_x - prev_x;

        all_s.row(row) = F.transpose();
    }
    double energy = x_div_x_bar.transpose() * m_M * x_div_x_bar;
    return 0.5 * energy;
}

void Solver::reset_fix_free_S_matrix(){

    const int dof = static_cast<int>(m_x.rows());
    const int pin_dof = 3 * static_cast<int>(m_constraints->pins.size());
    const int free_dof = dof - pin_dof;

    positive_pin.setOnes(static_cast<int>(dof/3));

    SparseMat S_fix(dof, pin_dof), S_free(dof, free_dof);
    Eigen::VectorXi nnz = Eigen::VectorXi::Ones( dof ); // non zeros per column
    S_fix.reserve(nnz);
    S_free.reserve(nnz);

    int count = 0;
    std::map<int,Vec3>::iterator pinIter = m_constraints->pins.begin();
    for( ; pinIter != m_constraints->pins.end(); ++pinIter ){
        //Construct Selection Matrix to select x
        int pin_id = pinIter->first;
        S_fix.coeffRef(3*pin_id, 3*count) = 1.0;
        S_fix.coeffRef(3*pin_id+1, 3*count+1) = 1.0;
        S_fix.coeffRef(3*pin_id+2, 3*count+2) = 1.0;
        positive_pin(pin_id) = 0;
        count++;
        //m_pin_energies[ pinIter->first ] = std::make_shared<SpringPin>( SpringPin(pinIter->first,pinIter->second) );
        //energyterms.emplace_back( m_pin_energies[ pinIter->first ] );
    }
    std::cout << " pin count = " << count << std::endl;

    count = 0;
    for( int i=0; i<positive_pin.size(); ++i )
    {
        if( positive_pin(i) == 1 )//Free point
        {
            S_free.coeffRef(3*i, 3*count) = 1.0;
            S_free.coeffRef(3*i+1, 3*count+1) = 1.0;
            S_free.coeffRef(3*i+2, 3*count+2) = 1.0;
            count++;
        }
    }

    m_S_fix = S_fix;
    m_S_free = S_free;

}

void Solver::set_pins( const std::vector<int> &inds, const std::vector<Vec3> &points ){

    int n_pins = static_cast<int>(inds.size());
    const int dof = static_cast<int>(m_x.rows());
    bool pin_in_place = (points.size() != inds.size());
    if( (dof == 0 && pin_in_place) || (pin_in_place && points.size() > 0) ){
        throw std::runtime_error("**Solver::set_pins Error: Bad input.");
    }

    if (m_x_pin.size() == 0)
    {
        m_x_pin.resize(n_pins*3);
    }

    m_constraints->pins.clear();
    for( size_t i=0; i<inds.size(); ++i ){
        int idx = inds[i];
        int pin_id = static_cast<int>(i)*3;
        if( pin_in_place ){
            m_constraints->pins[idx] = m_x_pin.segment<3>(pin_id);
        } else {
            m_constraints->pins[idx] = points[i];
            m_x_pin.segment<3>(pin_id) = points[i]; // We will not operate the pinned points.
        }
    }

    // If we're using energy based hard constraints, the pin locations may change
    // but which vertex is pinned may NOT change (aside from setting them to
    // active or inactive). So we need to do some extra work here.
    if( initialized ){
        //Update m_S_fix and m_S_free
        reset_fix_free_S_matrix();
    } // end set energy-based pins
}

void Solver::add_obstacle( std::shared_ptr<PassiveCollision> obj ){
    m_constraints->collider->add_passive_obj(obj);
}

void Solver::add_dynamic_collider( std::shared_ptr<DynamicCollision> obj ){
    m_constraints->collider->add_dynamic_obj(obj);
}

bool Solver::initialize( const Settings &settings_ ){
    using namespace Eigen;
    m_settings = settings_;

    mcl::MicroTimer t;
    const int dof = static_cast<int>(m_x.rows());
    if( m_settings.verbose > 0 ){ std::cout << "Solver::initialize: " << std::endl; }

    if( m_settings.timestep_s <= 0.0 ){
        std::cerr << "\n**Solver Error: timestep set to " << m_settings.timestep_s <<
                     "s, changing to 1/24s." << std::endl;
        m_settings.timestep_s = 1.0/24.0;
    }
    if( !( m_masses.rows()==dof && dof>=3 ) ){
        std::cerr << "\n**Solver Error: Problem with node data!" << std::endl;
        return false;
    }
    if( m_v.rows() != dof ){ m_v.resize(dof); }

    // Clear previous runtime stuff settings
    m_v.setZero();


    ///////////////////////////////////////////////
    // If we want energy-based constraints, set them up now.
    int pin_size = static_cast<int>(m_constraints->pins.size());
    int free_dof = dof - 3*pin_size;
    std::cout << "Energy term size = " << energyterms.size() << std::endl;
    std::cout << "pin size = " << pin_size << std::endl;

    Eigen::VectorXi pos_pin = Eigen::VectorXi::Ones( static_cast<int>(dof/3) );
    positive_pin = pos_pin;
    m_x_free.resize(free_dof);
    m_x_free.setZero();

    reset_fix_free_S_matrix();
    std::cout << "Initialized free and fix matrix." << std::endl;
    // end create energy based hard constraints
    ///////////////////////////////////////////////


    // Set up the selector matrix (D) and weight (W) matrix
    std::vector<Eigen::Triplet<double> > triplets;
    std::vector<double> weights;
    for(size_t i = 0; i < energyterms.size(); ++i){
        energyterms[i]->get_reduction( triplets, weights );
    }

    // Create the Selector+Reduction matrix
    m_W_diag = Eigen::Map<VecX>(&weights[0], static_cast<int>(weights.size()));
    int n_D_rows = static_cast<int>(weights.size());
    m_D.resize( n_D_rows, dof );
    m_D.setZero();
    m_D.setFromTriplets( triplets.begin(), triplets.end() );

    //    m_Dt = m_D.transpose();

    // Compute mass matrix
    SparseMat M( free_dof, free_dof ), Inv_M( free_dof, free_dof ); // Inv_M is just for test.
    Eigen::VectorXi nnz = Eigen::VectorXi::Ones( free_dof ); // non zeros per column
    M.reserve(nnz);
    Inv_M.reserve(nnz);
    double eps = 1e-6;
    int count = 0;
    for( int i=0; i<static_cast<int>(dof/3); ++i )
    {
        if (positive_pin(i) > 0)//free point
        {
            M.coeffRef(3*count,3*count) = m_masses[3*i];
            M.coeffRef(3*count+1,3*count+1) = m_masses[3*i+1];
            M.coeffRef(3*count+2,3*count+2) = m_masses[3*i+2];
            if (m_masses[3*i] > eps) {
                Inv_M.coeffRef(3*count,3*count) = 1.0/m_masses[3*i];
                Inv_M.coeffRef(3*count+1,3*count+1) = 1.0/m_masses[3*i+1];
                Inv_M.coeffRef(3*count+2,3*count+2) = 1.0/m_masses[3*i+2];
            } else {
                //                Inv_M.coeffRef(count,count) = 100000000.0;
                Inv_M.coeffRef(3*count,3*count) = 100000000.0;
                Inv_M.coeffRef(3*count+1,3*count+1) = 100000000.0;
                Inv_M.coeffRef(3*count+2,3*count+2) = 100000000.0;
            }
            count++;
        }

    }

    // Set global matrices
    SparseMat W( n_D_rows, n_D_rows ), W_inv( n_D_rows, n_D_rows );
    W.reserve(n_D_rows);
    W_inv.reserve(n_D_rows);
    for( int i=0; i<n_D_rows; ++i ){ W.coeffRef(i,i) = m_W_diag[i]; W_inv.coeffRef(i, i) = 1.0/m_W_diag[i]; }
    const double dt2 = (m_settings.timestep_s*m_settings.timestep_s);

    m_C = - W * m_D * m_S_fix;
    m_D = W * m_D * m_S_free;
    m_Dt = m_D.transpose();
    solver_Dt_Wt_W = dt2 * m_Dt;
    solver_termA = M + SparseMat(solver_Dt_Wt_W * m_D);

    solver_W = W;
    solver_W_inv = W_inv;
    dt2_ = dt2;
    m_M = M;

    // Set up the linear solver
    m_linsolver = std::make_shared<LDLTSolver>( LDLTSolver() );

    // If we haven't set a global solver, make one:
    if( !m_linsolver ){ throw std::runtime_error("What happened to the global solver?"); }
    if( m_settings.constraint_w > 0.0 ){ m_constraints->constraint_w = m_settings.constraint_w; }
    m_linsolver->update_system( solver_termA );

    // Make sure they don't have any collision obstacles
    if( m_constraints->collider->passive_objs.size() > 0 ||
            m_constraints->collider->dynamic_objs.size() > 0 ){
        throw std::runtime_error("**Solver::add_obstacle Error: No collisions with LDLT solver");
    }

    // All done
    if( m_settings.verbose >= 1 ){
        std::cout << m_x.size()/3 << " nodes, " << energyterms.size() << " energy terms" << std::endl;
    }
    initialized = true;
    return true;

} // end init


void Solver::save_matrix( const std::string &filename ){

    std::cout << "Saving matrix (" << solver_termA.rows() << "x" <<
                 solver_termA.cols() << ") to " << filename << std::endl;
    std::ofstream(filename.c_str()) << solver_termA;
}


template<typename T> void myclamp( T &val, T min, T max ){ if( val < min ){ val = min; } if( val > max ){ val = max; } }
bool Solver::Settings::parse_args( int argc, char **argv ){

    // Check args with params
    for( int i=1; i<argc-1; ++i ){
        std::string arg( argv[i] );
        std::stringstream val( argv[i+1] );
        if( arg == "-help" || arg == "--help" || arg == "-h" ){ help(); return true; }
        else if( arg == "-dt" ){ val >> timestep_s; }
        else if( arg == "-v" ){ val >> verbose; }
        else if( arg == "-it" ){ val >> admm_iters; }
        else if( arg == "-g" ){ val >> gravity; }
        //        else if( arg == "-ls" ){ val >> linsolver; }
        else if( arg == "-ck" ){ val >> constraint_w; }
        else if( arg == "-a" ){ int acc; val >> acc; (acc == 0) ? (acceleration_type = Solver::Settings::NOACC) : (acceleration_type = Solver::Settings::ANDERSON); }
        else if( arg == "-am" ){ val >> Anderson_m; acceleration_type = Solver::Settings::ANDERSON; }
        else if( arg == "-ab" ){ val >> beta; }
    }

    // Check if last arg is one of our no-param args
    std::string arg( argv[argc-1] );
    if( arg == "-help" || arg == "--help" || arg == "-h" ){ help(); return true; }

    return false;

} // end parse settings args

void Solver::Settings::help(){
    std::stringstream ss;
    ss << "\n==========================================\nArgs:\n" <<
          "\t-dt: time step (s)\n" <<
          "\t-v: verbosity (higher -> show more)\n" <<
          "\t-it: # admm iters\n" <<
          "\t-g: gravity (m/s^2)\n" <<
          "\t-ls: linear solver (0=LDLT, 1=NCMCGS, 2=UzawaCG) \n" <<
          "\t-ck: constraint weights (-1 = auto) \n" <<
          "\t-a: acceleration type (0=NoAcc, 1=Anderson) \n" <<
          "\t-am: anderson window size (>0, int) \n" <<
          "==========================================\n";
    printf( "%s", ss.str().c_str() );
}

void Solver::RuntimeData::print( const Settings &settings ){
    std::cout << "\nTotal global step: " << global_ms << "ms";
    std::cout << "\nTotal local step: " << local_ms << "ms";
    std::cout << "\nTotal acceleration step: " << acceleration_ms << "ms";
    std::cout << "\nTotal Initialization time: " << initialization_ms << "ms";
    std::cout << "\nAvg global step: " << global_ms/double(settings.admm_iters) << "ms";
    std::cout << "\nAvg local step: " << local_ms/double(settings.admm_iters) << "ms";
    std::cout << "\nAvg acceleration step: " << acceleration_ms/double(settings.admm_iters) << "ms";
    std::cout << "\nAvg Initialization step: " << initialization_ms/double(settings.admm_iters) << "ms";
    std::cout << "\nADMM Iters: " << settings.admm_iters;
    std::cout << "\nAvg Inner Iters: " << float(inner_iters) / float(settings.admm_iters);
    std::cout << "\nAnderson M: " << settings.Anderson_m;
    std::cout << std::endl;
}

