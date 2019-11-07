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

#ifndef ANDERSONACCELERATION_H_
#define ANDERSONACCELERATION_H_

#include "Types.h"
#include <cassert>
#include <algorithm>

class AndersonAcceleration {
 public:
  AndersonAcceleration(int m, int total_dim, int effective_dim)
      : m_(m),
        total_dim_(total_dim),
        effective_dim_(effective_dim),
        iter_(-1),
        col_idx_(-1) {
    assert(m_ > 0);
    current_u_.setZero(total_dim_);
    current_effective_F_.setZero(effective_dim_);
    prev_effective_dF_.setZero(effective_dim_, m_);
    effective_dF_scale_.setZero(m_);
    G_.setZero(total_dim_);
    prev_dG_.setZero(total_dim_, m_);
    M_.setZero(m_, m_);
    theta_.setZero(m_);

  }

  template<typename Derived>
  void replace(const Eigen::PlainObjectBase<Derived> &u) {
    current_u_ = Eigen::Map<const VectorX>(u.data(), u.size());
  }

  // The first argument must be the effective variable
  template<typename Derived1, typename Derived2>
  void replace(const Eigen::PlainObjectBase<Derived1> &u1,
               const Eigen::PlainObjectBase<Derived2> &u2) {
    current_u_.head(u1.size()) = Eigen::Map<const VectorX>(u1.data(),
                                                           u1.size());
    current_u_.tail(u2.size()) = Eigen::Map<const VectorX>(u2.data(),
                                                           u2.size());
  }

  template<typename Derived>
  void reset(const Eigen::PlainObjectBase<Derived> &u) {
    current_u_ = Eigen::Map<const VectorX>(u.data(), u.size());
    iter_ = 0;
    col_idx_ = 0;
  }

  // The first argument must be the effective variable
  template<typename Derived1, typename Derived2>
  void reset(const Eigen::PlainObjectBase<Derived1> &u1,
               const Eigen::PlainObjectBase<Derived2> &u2) {
    current_u_.head(u1.size()) = Eigen::Map<const VectorX>(u1.data(),
                                                           u1.size());
    current_u_.tail(u2.size()) = Eigen::Map<const VectorX>(u2.data(),
                                                           u2.size());

    iter_ = 0;
    col_idx_ = 0;
  }

  template<typename Derived>
  void compute(const Eigen::PlainObjectBase<Derived> &g,
               Eigen::PlainObjectBase<Derived> &accel_u) {
    G_ = Eigen::Map<const VectorX>(g.data(), g.size());
    compute_impl();
    Eigen::Map<VectorX>(accel_u.data(), accel_u.size()) = current_u_;
  }

  // The first argument must be the effective variable
  template<typename Derived1, typename Derived2>
  void compute(const Eigen::PlainObjectBase<Derived1> &g1,
               const Eigen::PlainObjectBase<Derived2> &g2,
               Eigen::PlainObjectBase<Derived1> &accel_u1,
               Eigen::PlainObjectBase<Derived2> &accel_u2) {
    G_.head(g1.size()) = Eigen::Map<const VectorX>(g1.data(), g1.size());
    G_.tail(g2.size()) = Eigen::Map<const VectorX>(g2.data(), g2.size());
    compute_impl();
    Eigen::Map<VectorX>(accel_u1.data(), accel_u1.size()) = current_u_.head(
        accel_u1.size());
    Eigen::Map<VectorX>(accel_u2.data(), accel_u2.size()) = current_u_.tail(
        accel_u2.size());
  }

  template<typename Derived>
  void init(const Eigen::PlainObjectBase<Derived> &init_u) {
    assert(int(init_u.size()) == total_dim_);
    current_u_ = Eigen::Map<const VectorX>(init_u.data(), init_u.size());
    iter_ = 0;
    col_idx_ = 0;
  }

  // The first argument must be the effective variable
  template<typename Derived1, typename Derived2>
  void init(const Eigen::PlainObjectBase<Derived1> &init_u1,
            const Eigen::PlainObjectBase<Derived2> &init_u2) {
    assert(int(init_u1.size() + init_u2.size()) == total_dim_);
    current_u_.head(init_u1.size()) = Eigen::Map<const VectorX>(init_u1.data(),
                                                                init_u1.size());
    current_u_.tail(init_u2.size()) = Eigen::Map<const VectorX>(init_u2.data(),
                                                                init_u2.size());
    iter_ = 0;
    col_idx_ = 0;
  }

 private:
  int m_;   // Number of previous iterates used for Anderson Acceleration
  int total_dim_;  // Total dimension of all varaibles
  int effective_dim_;  // Dimension of variables used for computing the linear combination coefficients
  int iter_;  // Iteration count since initialization
  int col_idx_;  // Index for history matrix column to store the next value

  VectorX current_u_;
  VectorX current_effective_F_;
  MatrixXX prev_effective_dF_;
  VectorX effective_dF_scale_;  // The scaling factor for each column of prev_effective_dF_
  VectorX G_;    // Storage for the latest value of g
  MatrixXX prev_dG_;
  MatrixXX M_;    // Normal equations matrix for the computing theta
  VectorX theta_;  // theta value computed from normal equations
  Eigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;

  void compute_impl() {
    assert(iter_ >= 0);
    current_effective_F_ = G_.head(effective_dim_)
        - current_u_.head(effective_dim_);

    if (iter_ == 0) {
      prev_effective_dF_.col(0) = -current_effective_F_;
      prev_dG_.col(0) = -G_;
      current_u_ = G_;
    } else {
      prev_effective_dF_.col(col_idx_) += current_effective_F_;
      prev_dG_.col(col_idx_) += G_;

      Scalar eps = 1e-14;
      Scalar scale = std::max(eps, prev_effective_dF_.col(col_idx_).norm());
      effective_dF_scale_(col_idx_) = scale;
      prev_effective_dF_.col(col_idx_) /= scale;

      int m_k = std::min(m_, iter_);

      if (m_k == 1) {
        theta_(0) = 0;
        Scalar dF_sqrnorm = prev_effective_dF_.col(col_idx_).squaredNorm();
        M_(0, 0) = dF_sqrnorm;
        Scalar dF_norm = std::sqrt(dF_sqrnorm);

        if (dF_norm > eps) {
          // compute theta = (dF * F) / (dF * dF)
          theta_(0) = (prev_effective_dF_.col(col_idx_) / dF_norm).dot(
              current_effective_F_ / dF_norm);
        }
      } else {
        // Update the normal equation matrix, for the column and row corresponding to the new dF column
        VectorX new_inner_prod = (prev_effective_dF_.col(col_idx_).transpose()
            * prev_effective_dF_.block(0, 0, effective_dim_, m_k)).transpose();
        M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
        M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

        // Solve normal equation
        cod_.compute(M_.block(0, 0, m_k, m_k));
        theta_.head(m_k) = cod_.solve(
            prev_effective_dF_.block(0, 0, effective_dim_, m_k).transpose()
                * current_effective_F_);
      }

      // Use rescaled theata to compute new u
      current_u_ = G_
          - prev_dG_.block(0, 0, total_dim_, m_k)
              * (theta_.head(m_k).array()
                  / effective_dF_scale_.head(m_k).array()).matrix();

      col_idx_ = (col_idx_ + 1) % m_;
      prev_effective_dF_.col(col_idx_) = -current_effective_F_;
      prev_dG_.col(col_idx_) = -G_;
    }

    iter_++;
  }
};

#endif /* ANDERSONACCELERATION_H_ */
