//  BSD 3-Clause License
//
//  Copyright (c) 2018, Bailin Deng, Yue Peng
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

#ifndef GEOMETRYSOLVER_H
#define GEOMETRYSOLVER_H

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
class GeometrySolver {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;
  typedef MatrixT<Eigen::Dynamic, N> MatrixXN;

 public:
  GeometrySolver()
      : current_x_(NULL),
        default_x_(NULL),
        current_u_(NULL),
        default_u_(NULL),
        n_hard_constraints_(0),
        n_soft_constraints_(0),
        penalty_parameter_(1.0),
        ADMM_solver_initialized_(false) {
  }

  ~GeometrySolver() {
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

    x_block_1_.setZero(N, n_points);
    x_block_2_.setZero(N, n_points);

    std::vector<Triplet> triplets;
    int idO = 0;
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->add_constraint(false, triplets, idO);
    }

    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->add_constraint(false, triplets, idO);
    }

    int z_cols = idO;
    z_.setZero(N, z_cols);
    u_block_1_.setZero(N, z_cols);
    u_block_2_.setZero(N, z_cols);
    Dx_.setZero(N, z_cols);
    Dx_plus_u_.setZero(N, z_cols);

    // Set up full global update matrix
    ColMajorSparseMatrix D(idO, n_points);
    D.setFromTriplets(triplets.begin(), triplets.end());
    D.makeCompressed();
    rho_Dt_ = D.transpose() * penalty_parameter_;
    ColMajorSparseMatrix global_mat = rho_Dt_ * D;

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

    ADMM_solver_initialized_ = true;

    return true;
  }

  void solve_ADMM(const MatrixNX &init_x, Scalar rel_residual_eps, int max_iter,
                  int Anderson_m) {
    if (!ADMM_solver_initialized_) {
      std::cerr << "Error: solver not initialized yet" << std::endl;
    }

    clear_iteration_history();
    ADMM_init_variables(init_x);
    int iter_count = 0;
    Scalar residual_eps = rel_residual_eps * z_.cols();
    Scalar prev_residual = std::numeric_limits<Scalar>::max();
    Scalar current_residual = 0;
    bool need_reset = false;
    bool end_iteration = false;
    AndersonAcceleration* aa = NULL;
    if (Anderson_m > 0) {
      aa = new AndersonAcceleration(Anderson_m,
                                    current_x_->size() + current_u_->size(),
                                    current_u_->size());
      aa->init(*current_u_, *current_x_);
    }

    Timer timer;
    Timer::EventID t_begin = timer.get_time();

    while (!end_iteration) {
      OMP_PARALLEL
      {
        ADMM_z_update();

        OMP_SINGLE
        {
          current_residual = get_ADMM_residual();
          need_reset = Anderson_m > 0 && current_residual > prev_residual;
        }

        if (need_reset) {
          OMP_SINGLE
          {
            std::cout << "AA reset" << std::endl;
            std::swap(current_x_, default_x_);
            std::swap(current_u_, default_u_);
            aa->replace(*current_u_, *current_x_);
          }

          ADMM_compute_Dx(*current_x_);

          ADMM_z_update();

          OMP_SINGLE
          {
            current_residual = get_ADMM_residual();
          }
        }

        OMP_SINGLE
        {
          iter_count++;
          //end_iteration = current_residual <= residual_eps
           //   || iter_count >= max_iter;
          end_iteration = iter_count >= max_iter;

          function_values_.push_back(current_residual);
          Timer::EventID t_iter = timer.get_time();
          elapsed_time_.push_back(timer.elapsed_time(t_begin, t_iter));

          std::cout << "Iteration " << iter_count << ",  residual: "
                    << current_residual << ", threshold: " << residual_eps
                    << std::endl;
        }

        if (!end_iteration) {
          OMP_SINGLE
          {
            prev_residual = current_residual;
          }

          ADMM_x_update();

          ADMM_compute_Dx(*default_x_);

          ADMM_u_update();

          OMP_SINGLE
          {
            Scalar primal_residual = get_ADMM_residual();
            std::cout << "Primal residual: " << primal_residual << std::endl;
          }

          OMP_SINGLE
          {
            if (aa) {
              aa->compute(*default_u_, *default_x_, *current_u_, *current_x_);
            } else {
              std::swap(default_u_, current_u_);
              std::swap(default_x_, current_x_);
            }
          }

          ADMM_compute_Dx(*current_x_);
        }
      }
    }

    if (aa) {
      delete aa;
    }
  }

  const MatrixNX& get_solution() {
    return *current_x_;
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
  MatrixNX x_block_1_, x_block_2_, u_block_1_, u_block_2_;	// Vertex coordinate matrices
  MatrixNX *current_x_, *default_x_, *current_u_, *default_u_;
  MatrixNX Dx_, Dx_plus_u_, z_;

  // Constraints and regularization terms
  std::vector<Constraint<N>*> soft_constraints_;
  std::vector<Constraint<N>*> hard_constraints_;
  int n_hard_constraints_, n_soft_constraints_;

  Scalar penalty_parameter_;

  LinearRegularization<N> regularization_;

  // Data structures used for direct solve update
  std::shared_ptr<SPDSolver> SPD_solver_;

  MatrixXN x_system_rhs_;
  ColMajorSparseMatrix rho_Dt_;
  MatrixXN rhs_fixed_;

 public:
  // History of iterations
  std::deque<bool> Anderson_reset_;
  std::vector<double> function_values_;
  std::vector<double> elapsed_time_;

 private:
  bool ADMM_solver_initialized_;

  void clear_iteration_history() {
    Anderson_reset_.clear();
    function_values_.clear();
    elapsed_time_.clear();
  }

  void ADMM_init_variables(const MatrixNX &init_x) {
    x_block_1_ = init_x;
    x_block_2_ = init_x;
    u_block_1_.setZero();
    u_block_2_.setZero();

    current_x_ = &x_block_1_;
    default_x_ = &x_block_2_;
    current_u_ = &u_block_1_;
    default_u_ = &u_block_2_;

    OMP_PARALLEL
    {
      ADMM_compute_Dx(*current_x_);

      ADMM_z_update();

      ADMM_x_update();

      ADMM_compute_Dx(*default_x_);

      ADMM_u_update();
    }

    (*current_x_) = (*default_x_);
    (*current_u_) = (*default_u_);
  }

  void ADMM_compute_Dx(const MatrixNX &x) {
    OMP_FOR
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->apply_transform(x, Dx_);
    }

    OMP_FOR
    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->apply_transform(x, Dx_);
    }
  }

  void ADMM_z_update() {
    OMP_SINGLE
    {
      Dx_plus_u_ = Dx_ + (*current_u_);
    }

    OMP_FOR
    for (int i = 0; i < n_hard_constraints_; ++i) {
      hard_constraints_[i]->project(Dx_plus_u_, z_);
    }

    OMP_FOR
    for (int i = 0; i < n_soft_constraints_; ++i) {
      soft_constraints_[i]->project_and_combine(Dx_plus_u_, penalty_parameter_,
                                                z_);
    }
  }

  void ADMM_x_update() {
    OMP_FOR
    for (int i = 0; i < int(N); ++i) {
      x_system_rhs_.col(i) = rhs_fixed_.col(i)
          + rho_Dt_ * (z_.row(i) - current_u_->row(i)).transpose();
      default_x_->row(i) = SPD_solver_->solve(x_system_rhs_.col(i)).transpose();
    }
  }

  void ADMM_u_update() {
    OMP_SINGLE
    {
      (*default_u_) = (*current_u_) + Dx_ - z_;
    }
  }

  Scalar get_ADMM_residual() {
    return (Dx_ - z_).norm();
  }
};

#endif	// GEOMETRYSOLVER_H
