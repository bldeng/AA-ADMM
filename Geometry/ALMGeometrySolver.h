//  BSD 3-Clause License
//
//  Copyright (c) 2019, Bailin Deng
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
//  * Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef ALMGEOMETRYSOLVER_H
#define ALMGEOMETRYSOLVER_H

#include "SolverCommon.h"
#include "Constraint.h"
#include "SPDSolver.h"
#include "LinearRegularization.h"
#include "OMPHelper.h"
#include "AndersonAcceleration.h"
#include "Parameters.h"
#include <vector>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <limits>
#include <deque>
#include <iomanip>
#include <set>

template<unsigned int N>
class ALMGeometrySolver {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;
  typedef MatrixT<Eigen::Dynamic, N> MatrixXN;

 public:
  ALMGeometrySolver()
      : n_hard_constraints_(0),
        n_soft_constraints_(0),
        penalty_parameter_(1.0),
        solver_initialized_(false) {
  }

  ~ALMGeometrySolver() {
    for (int i = 0; i < static_cast<int>(hard_constraints_.size()); ++i) {
      if (hard_constraints_[i]) {
        delete hard_constraints_[i];
      }
    }

    for (int i = 0; i < static_cast<int>(soft_constraints_.size()); ++i) {
      if (soft_constraints_[i]) {
        delete soft_constraints_[i];
      }
    }
  }

  bool setup_ADMM(int n_points, Scalar penalty_param,
                  SPDSolverType spd_solver_type = LDLT_SOLVER) {
    Timer timer;
    Timer::EventID t_begin = timer.get_time();

    penalty_parameter_ = penalty_param;
    n_hard_constraints_ = static_cast<int>(hard_constraints_.size());
    n_soft_constraints_ = static_cast<int>(soft_constraints_.size());

    assert(n_points != 0);
    assert((n_hard_constraints_ + n_soft_constraints_) != 0);

    current_x_.setZero(N, n_points);
    new_x_.setZero(N, n_points);

    std::vector<Triplet> triplets;
    int idO = 0;
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->add_constraint(false, triplets, idO);
    }

    // Set up full global update matrix
    ColMajorSparseMatrix D_hard(idO, n_points);
    D_hard.setFromTriplets(triplets.begin(), triplets.end());
    D_hard.makeCompressed();
    rho_D_hard_t_ = D_hard.transpose() * penalty_parameter_;
    ColMajorSparseMatrix global_mat = rho_D_hard_t_ * D_hard;
    int z_cols = D_hard.rows();
    std::cout << "z_cols = " << z_cols*N << std::endl;

    z_hard_.setZero(N, z_cols);
    current_u_.setZero(N, z_cols);
    new_u_.setZero(N, z_cols);
    Dx_hard_.setZero(N, z_cols);

    idO = 0;
    triplets.clear();
    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->add_constraint(true, triplets, idO);
    }
    std::cout << "n_soft_constraints_ = " << n_soft_constraints_ << std::endl;

    ColMajorSparseMatrix D_soft(idO, n_points);
    D_soft.setFromTriplets(triplets.begin(), triplets.end());
    D_soft.makeCompressed();
    D_soft_t_ = D_soft.transpose();
    global_mat += D_soft_t_ * D_soft;
    Dx_soft_.setZero(N, D_soft.rows());
    z_soft_.setZero(N, D_soft.rows());

    // Set up regularization terms
    rhs_fixed_.setZero(n_points, N);
    x_system_rhs_.setZero(n_points, N);
    MatrixXN regularization_rhs;
    ColMajorSparseMatrix L;
    if (regularization_.get_regularization_system(n_points, L,
                                                  regularization_rhs)) {
      ColMajorSparseMatrix Lt = L.transpose();
      global_mat += Lt * L;
      rhs_fixed_ = Lt * regularization_rhs;
    }

    if (spd_solver_type == LDLT_SOLVER) {
      SPD_solver_ = std::make_shared<SimplicialLDLTSolver>();
    } else {
      SPD_solver_ = std::make_shared<SimplicialLLTSolver>();
    }

    SPD_solver_->initialize(global_mat);
    if (SPD_solver_->info() != Eigen::Success) {
      std::cerr << "Error: SPD solver initialization failed" << std::endl;
      return false;
    }

    std::cout << "predecomposition time = "
              << timer.elapsed_time(t_begin, timer.get_time()) << std::endl;

    solver_initialized_ = true;

    return true;
  }

  void solve_ADMM(const MatrixNX &init_x, Scalar rel_residual_eps, int max_iter,
                  int Anderson_m) {
    if (!solver_initialized_) {
      std::cerr << "Error: solver not initialized yet" << std::endl;
    }

    bool accel = Anderson_m > 0;

    init_variables(init_x);
    int iter_count = 0;
    Scalar residual_eps = rel_residual_eps * rel_residual_eps *  z_hard_.cols() * z_hard_.cols() * 2;
    Scalar prev_residual = std::numeric_limits<Scalar>::max();
    Scalar current_residual = 0;
    bool accept = false, reset = false;
    bool end_iteration = false;
    AndersonAcceleration *aa = NULL;

    if(accel){
      aa = new AndersonAcceleration(Anderson_m, current_x_.size() + current_u_.size(),
    		  current_x_.size() + current_u_.size());
      aa->init(current_u_, current_x_);
    }

    Scalar constraint_err_;
    OMP_PARALLEL
    {
      ADMM_compute_Dx_soft(current_x_);
    }
    soft_constraints_[0]->project(Dx_soft_, z_soft_, true, &constraint_err_);
    std::cout << "Init energy = " << constraint_err_ << std::endl;

    Timer timer;
    Timer::EventID t_begin = timer.get_time();

    while (!end_iteration) {
      OMP_PARALLEL
      {
        ADMM_compute_Dx_hard(current_x_);
        ADMM_compute_Dx_soft(current_x_);
        OMP_SINGLE
		{
        	prev_Dx_hard_ = Dx_hard_;
		}

        ADMM_z_update();

        ADMM_x_update();

        ADMM_compute_Dx_hard(new_x_);

        ADMM_u_update();

        OMP_SINGLE
        {
          current_residual = get_combined_residual();

          accept = (!accel) || reset || current_residual < prev_residual;

           if(accept) {
            default_x_ = new_x_;
            default_u_ = new_u_;
            iter_count ++;
            //std::cout << "Iteration " << iter_count << ":   ";
            //std::cout << "combined residual: " << current_residual << ", threshold: " << residual_eps;

//            if(reset){
//              std::cout << ", AA was reset";
//            }
//            std::cout << std::endl;

            function_values_.push_back(current_residual);
            Timer::EventID t_iter = timer.get_time();
            elapsed_time_.push_back(timer.elapsed_time(t_begin, t_iter));

            prev_residual = current_residual;
            reset = false;

            if(accel){
              aa->compute(new_u_, new_x_, current_u_, current_x_);
            }
            else{
              current_u_ = new_u_;
              current_x_ = new_x_;
            }
          }
           else{
             current_u_ = default_u_;
             current_x_ = default_x_;
             reset = true;

             if(accel){
               aa->reset(current_u_, current_x_);
             }
           }

//           if(iter_count >= max_iter || current_residual < residual_eps){
//             end_iteration = true;
//           }
           if(iter_count >= max_iter){
             end_iteration = true;
           }


        }
      }
    }
    //Scalar constraint_err_;

    OMP_PARALLEL
    {
      ADMM_compute_Dx_soft(current_x_);
    }

    soft_constraints_[0]->project(Dx_soft_, z_soft_, true, &constraint_err_);

    std::cout << "final energy = " << constraint_err_ << std::endl;

    if (aa) {
      delete aa;
    }
  }

  const MatrixNX& get_solution() {
    return default_x_;
  }

  void add_hard_constraint(Constraint<N>* constraint) {
    hard_constraints_.push_back(constraint);
  }

  void add_soft_constraint(Constraint<N>* constraint) {
    soft_constraints_.push_back(constraint);
  }

  void add_closeness(int idx, Scalar weight, const VectorN &target_pt) {
    regularization_.add_closeness(idx, weight, target_pt);
  }

  void add_uniform_laplacian(const std::vector<int> &indices, Scalar weight) {
    regularization_.add_uniform_laplacian(indices, weight);
  }

  void add_laplacian(const std::vector<int> &indices,
                     const std::vector<Scalar> coefs, Scalar weight) {
    regularization_.add_laplacian(indices, coefs, weight);
  }

  void add_relative_uniform_laplacian(const std::vector<int> &indices,
                                      Scalar weight,
                                      const MatrixNX &ref_points) {
    regularization_.add_relative_uniform_laplacian(indices, weight, ref_points);
  }

  void add_relative_laplacian(const std::vector<int> &indices,
                              const std::vector<Scalar> coefs, Scalar weight,
                              const MatrixNX &ref_points) {
    regularization_.add_relative_laplacian(indices, coefs, weight, ref_points);
  }

  void output_iteration_history(SolverType solver_type) {
    if (solver_type == AA_SOLVER)
      assert(function_values_.size() == Anderson_reset_.size());

    assert(function_values_.size() == elapsed_time_.size());
    int n_iter = function_values_.size();
    for (int i = 0; i < n_iter; ++i) {
      std::cout << "Iteration " << i << ": ";
      std::cout << std::setprecision(6) << elapsed_time_[i] << " secs, ";
      std::cout << " target value " << std::setprecision(16)
                << function_values_[i];
      if ((solver_type == AA_SOLVER) && (Anderson_reset_[i])) {
        std::cout << " (reject accelerator)";
      }

      std::cout << std::endl;
    }

    std::cout << std::endl;
  }

  void save(int Anderson_m)
  {
      std::string file;
      if (Anderson_m > 0)
          file = "./result/residual-"+std::to_string(Anderson_m)+".txt";
      else
          file = "./result/residual-no.txt";

      std::ofstream ofs;
      ofs.open(file, std::ios::out | std::ios::ate);
      if (!ofs.is_open())
      {
          std::cout << "Cannot open: " << file << std::endl;
      }

      ofs << std::setprecision(16);
      for (size_t i = 0; i < elapsed_time_.size(); i++)
      {
          ofs << elapsed_time_[i] << '\t' << function_values_[i] << std::endl;
      }

      ofs.close();
  }

 protected:
  MatrixNX current_x_, current_u_, new_x_, new_u_;
  MatrixNX default_x_, default_u_;
  MatrixNX Dx_hard_, Dx_soft_, z_hard_, z_soft_;
  MatrixNX prev_Dx_hard_;

  // Constraints and regularization terms
  std::vector<Constraint<N>*> soft_constraints_;
  std::vector<Constraint<N>*> hard_constraints_;
  int n_hard_constraints_, n_soft_constraints_;

  Scalar penalty_parameter_;

  LinearRegularization<N> regularization_;

  // Data structures used for direct solve update
  std::shared_ptr<SPDSolver> SPD_solver_;

  MatrixXN x_system_rhs_;
  ColMajorSparseMatrix rho_D_hard_t_, D_soft_t_;
  MatrixXN rhs_fixed_;

 public:
  // History of iterations
  std::deque<bool> Anderson_reset_;
  std::vector<double> function_values_;
  std::vector<double> elapsed_time_;

 private:
  bool solver_initialized_;

  void clear_iteration_history() {
    Anderson_reset_.clear();
    function_values_.clear();
    elapsed_time_.clear();
  }

  void init_variables(const MatrixNX &init_x) {
    current_x_ = init_x;
    default_x_ = init_x;
    current_u_.setZero();
    default_u_.setZero();
  }

  void ADMM_compute_Dx_hard(const MatrixNX &x) {
    OMP_FOR
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->apply_transform(x, Dx_hard_);
    }
  }

  void ADMM_compute_Dx_soft(const MatrixNX &x){
    OMP_FOR
    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->apply_transform(x, Dx_soft_);
    }
  }

  void ADMM_z_update() {
    OMP_SINGLE
    {
      Dx_hard_ += current_u_;
    }

    OMP_FOR
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->project(Dx_hard_, z_hard_);
    }

    OMP_FOR
    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->project(Dx_soft_, z_soft_, true);
    }
  }

  void ADMM_x_update() {
    OMP_FOR
    for (int i = 0; i < int(N); ++i) {
      x_system_rhs_.col(i) = rhs_fixed_.col(i)
          + rho_D_hard_t_ * (z_hard_.row(i) - current_u_.row(i)).transpose()
          + D_soft_t_ * z_soft_.row(i).transpose();
      new_x_.row(i) = SPD_solver_->solve(x_system_rhs_.col(i)).transpose();
    }
  }

  void ADMM_u_update() {
    OMP_SINGLE
    {
      new_u_ = current_u_ + Dx_hard_ - z_hard_;
    }
  }

  Scalar get_combined_residual() {
    return (Dx_hard_ - z_hard_).squaredNorm() + (Dx_hard_ - prev_Dx_hard_).squaredNorm();
  }

};

#endif	// ALMGEOMETRYSOLVER_H
