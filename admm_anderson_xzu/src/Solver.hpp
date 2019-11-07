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

#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H 1

#include "ConstraintSet.hpp"
#include "EnergyTerm.hpp"
#include "SpringEnergyTerm.hpp"
#include "ExplicitForce.hpp"
#include "LinearSolver.hpp"
#include "AndersonAcceleration.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>

namespace admm {

// The main solver
class Solver {
public:
	typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
	typedef Eigen::Matrix<double,3,1> Vec3;
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMat;

	// Solver settings
	struct Settings {
        enum AccelationType
        {
          NOACC = 0,
          ANDERSON = 1,
        };
		bool parse_args( int argc, char **argv ); // parse from terminal args. Returns true if help()
		void help();		// -help	print details, parse_args returns true if used
		double timestep_s;	// -dt <flt>	timestep in seconds
		int verbose;		// -v <int>	terminal output level (higher=more)
		int admm_iters;		// -it <int>	number of admm solver iterations
		double gravity;		// -g <flt>	force of (-y) gravity
        //int linsolver;		// -ls <int>	0=LDLT, 1=NCMCGS, 2=UzawaCG
		double constraint_w;	// -ck <flt>	constraint weights (-1 = auto)
        int Anderson_m;
        double beta;
        AccelationType acceleration_type;
        Settings() : timestep_s(1.0/30.0), verbose(1), admm_iters(500),
            gravity(-9.8), constraint_w(-1), Anderson_m(2),
            beta(1.0),
            acceleration_type(AccelationType::NOACC) {}
	};

	// RuntimeData struct used for logging.
	// Add timings are per time step.
	struct RuntimeData {
		double global_ms; // total ms for global solver
		double local_ms; // total ms for local solver
        double acceleration_ms; // total ms for Anderson acceleration
        double initialization_ms; // total ms for collision update/detection
		int inner_iters; // total global step iterations
        std::vector<double> step_time;
        RuntimeData() : global_ms(0), local_ms(0), acceleration_ms(0), initialization_ms(0), inner_iters(0) {}
		void print( const Settings &settings );
	};

	Solver();

	// Per-node (x3) data (for x, y, and z)
	VecX m_x; // node positions, scaled x3
	VecX m_v; // node velocities, scaled x3
	VecX m_masses; // node masses, scaled x3
	std::vector<int> surface_inds; // indices of surface vertices

	std::vector< std::shared_ptr<ExplicitForce> > ext_forces; // external/explicit forces
	std::vector< std::shared_ptr<EnergyTerm> > energyterms; // minimized (implicit)

	// Adds nodes to the Solver.
	// Returns the current total number of nodes after insert.
	// Assumes m is scaled x3 (i.e. 3 values per node).
	template <typename T>
	int add_nodes( T *x, T *m, int n_verts );

    void reset_fix_free_S_matrix();

	// Pins vertex indices to the location indicated. If the points
	// vector is empty (or not the same size as inds), vertices are pinned in place.
	virtual void set_pins( const std::vector<int> &inds,
		const std::vector<Vec3> &points = std::vector<Vec3>() );

	// An obstacle is a passive collision object.
	virtual void add_obstacle( std::shared_ptr<PassiveCollision> obj );

	// A dynamic obstacle has vertices in m_x and is updated every frame
	virtual void add_dynamic_collider( std::shared_ptr<DynamicCollision> obj );

	// Returns true on success.
	virtual bool initialize( const Settings &settings_=Settings() );

	// Performs a Solver step
	virtual void step();

	// Returns the runtime data from the last time step.
	virtual const RuntimeData &runtime_data(){ return m_runtime; }

	// Outputs solver_termA to a file (for debugging/analysis).
	virtual void save_matrix( const std::string &filename );

	// Returns the current settings
	const Settings &settings(){ return m_settings; }

    void save()
    {
        std::string file;
        if (m_settings.acceleration_type )
            file = "./result/residual-"+std::to_string(m_settings.Anderson_m)+".txt";
        else
            file = "./result/residual-no.txt";

        std::ofstream ofs;
        ofs.open(file, std::ios::out | std::ios::ate);
        if (!ofs.is_open())
        {
            std::cout << "Cannot open: " << file << std::endl;
        }

        ofs << std::setprecision(16);
        for (size_t i = 0; i < step_prim_residual.size(); i++)
        {
            ofs << m_runtime.step_time[i] << '\t' << step_prim_residual[i] << '\t' << step_comb_residual[i] << std::endl;
        }

        ofs.close();
        step_prim_residual.resize(0);
        step_comb_residual.resize(0);
        m_runtime.step_time.resize(0);
    }

    static bool load(const char *file_name, const char *file_name2, VecX &curr_z, VecX &curr_u, VecX &last_z, VecX &curr_x)
    {
        std::ifstream ifile(file_name), ifile2(file_name2);
        if (!ifile.is_open()){
            std::cerr << "Unable to open file " << file_name << std::endl;
            return false;
        }
        if (!ifile2.is_open()){
            std::cerr << "Unable to open file " << file_name2 << std::endl;
            return false;
        }

        int n_values = 0;
        if(!(ifile >> n_values)){
            std::cerr << "Error parsing the number of values" << std::endl;
            return false;
        }

        if((n_values <= 0) || (n_values != curr_z.size())){
            std::cerr << "Error: invalid number or values" << std::endl;
            return false;
        }

        ifile.precision(16);

        double val1 = 0, val2 = 0, val3 = 0;
        for(int i = 0; i < n_values; ++ i){
            if(ifile >> val1 >> val2 >> val3){
                curr_z(i) = val1;
                curr_u(i) = val2;
                last_z(i) = val3;
            }
            else{
                std::cerr << "Error parsing distance value at position " << i << std::endl;
                return false;
            }
        }

        n_values = 0;
        if(!(ifile2 >> n_values)){
            std::cerr << "Error parsing the number of values from file 2" << std::endl;
            return false;
        }

        if((n_values <= 0) || (n_values != curr_x.size())){
            std::cerr << "Error: invalid number or values from file 2" << std::endl;
            return false;
        }

        ifile2.precision(16);

        for(int i = 0; i < n_values; ++ i){
            if(ifile2 >> val3){
                curr_x(i) = val3;
            }
            else{
                std::cerr << "Error parsing distance value at position " << i << std::endl;
                return false;
            }
        }

        return true;
    }

    void reset_settings(Settings::AccelationType flag_is_acc)
    {
        m_settings.acceleration_type = flag_is_acc;
        std::cout << "m_settings.acceleration_type = " << ((m_settings.acceleration_type == Settings::ANDERSON) ? "anderson." : "no acceleration.") << std::endl;
    }

    bool isNan(double fN)
    {
        return !(fN==fN);
    }

    double Calc_function(const VecX &curr_x);

    double Calc_function(MatrixXX &all_s, const VecX &curr_x, const VecX &prev_x, int row);

protected:
	Settings m_settings; // copied from init
	RuntimeData m_runtime; // reset each iteration
	bool initialized; // has init been called?

	// Solver used in the global step
	std::shared_ptr<LinearSolver> m_linsolver;
	std::shared_ptr<ConstraintSet> m_constraints;
	std::unordered_map<int, std::shared_ptr<SpringPin> > m_pin_energies;

	// Global matrices
    SparseMat m_D, m_Dt, m_S_fix, m_S_free, m_C; // reduction matrix
	VecX m_W_diag; // diagonal of the weight matrix

	// Solver variables computed in initialize
	SparseMat solver_termA;
    SparseMat solver_Dt_Wt_W;

    SparseMat solver_W, solver_W_inv, m_M;
    std::vector<double> step_prim_residual, step_comb_residual;
    double dt2_ = 0.0;

    //For pin
    Eigen::VectorXi positive_pin;

    VecX m_x_pin; // node positions, scaled x3
    VecX m_x_free;
}; // end class Solver


template <typename T>
int Solver::add_nodes( T *x, T *m, int n_verts ){
	int prev_n = m_x.size();
	int n_verts3 = n_verts*3;
	m_x.conservativeResize(prev_n+n_verts3);
	m_v.conservativeResize(prev_n+n_verts3);
	m_masses.conservativeResize(prev_n+n_verts3);
	for( int i=0; i<n_verts3; ++i ){
		int idx = prev_n+i;
		m_x[idx] = x[i];
		m_v[idx] = 0.0;
		m_masses[idx] = m[i];
	}
	return (prev_n+n_verts3)/3;
}


} // end namespace admm

#endif




