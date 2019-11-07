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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Types.h"
#include "TriMeshAABB.h"
#include <igl/svd3x3.h>
#include <igl/AABB.h>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

#ifndef M_PI
/** Defining pi.*/
#define M_PI 3.14159265358979323846
#endif

template<unsigned int N>
class Constraint {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;

  // The type of transform before projection
  enum InvarianceTransformType {
    MEAN_CENTERING,
    SUBTRACT_FIRST,
    IDENTITY,
  };

 public:
  Constraint(const std::vector<int> &idI, Scalar weight,
             InvarianceTransformType transform_type)
      : idI_(Eigen::Map<const Eigen::VectorXi>(idI.data(), idI.size())),
        weight_(std::sqrt(weight)),
        transform_type_(transform_type),
        idO_(-1) {
  }

  virtual ~Constraint() {
  }

  void apply_transform(const MatrixNX &positions,
                       MatrixNX &transformed_positions) {
    if (transform_type_ == SUBTRACT_FIRST) {
      VectorN first_pt = positions.col(idI_[0]);

      for (int i = 1; i < static_cast<int>(idI_.size()); ++i) {
        transformed_positions.col(idO_ + i - 1) = positions.col(idI_[i])
            - first_pt;
      }
    } else {
      for (int i = 0; i < static_cast<int>(idI_.size()); ++i) {
        transformed_positions.col(idO_ + i) = positions.col(idI_[i]);
      }

      if (transform_type_ == MEAN_CENTERING) {
        VectorN mean_vector = transformed_positions.block(
            0, idO_, N, num_transformed_points()).rowwise().mean();
        transformed_positions.block(0, idO_, N, num_transformed_points())
            .colwise() -= mean_vector;
      }
    }
  }

  virtual void project(const MatrixNX& transformed_positions,
                       MatrixNX& projections, bool weighted = false,
                       Scalar* squared_deviation = NULL) {
    int num_block_cols = num_transformed_points();
    Eigen::Map<const MatrixNX> input_block(&(transformed_positions(0, idO_)), N,
                                           num_block_cols);
    Eigen::Map<MatrixNX> output_block(&(projections(0, idO_)), N,
                                      num_block_cols);
    project_impl(input_block, output_block);

    if (squared_deviation) {
      (*squared_deviation) = (input_block - output_block).squaredNorm() * 0.5;
    }

    if (weighted) {
      output_block *= weight_;
      if (squared_deviation) {
        (*squared_deviation) *= weight_ * weight_;
      }
    }
  }

  virtual void project_and_combine(const MatrixNX& transformed_positions,
                                   Scalar penalty_param,
                                   MatrixNX& projections) {
    project(transformed_positions, projections, false, NULL);
    int num_block_cols = num_transformed_points();
    Eigen::Map<const MatrixNX> input_block(&(transformed_positions(0, idO_)), N,
                                           num_block_cols);
    Eigen::Map<MatrixNX> output_block(&(projections(0, idO_)), N,
                                      num_block_cols);
    Scalar w = weight_ * weight_;
    Scalar a = penalty_param / (w + penalty_param);
    output_block = input_block * a + output_block * (1 - a);
  }

  virtual void add_constraint(bool weighted, std::vector<Triplet> &triplets,
                              int &idO) {
    idO_ = idO;
    int n_idx = static_cast<int>(idI_.size());
    Scalar w = weighted ? weight_ : Scalar(1.0);

    if (transform_type_ == MEAN_CENTERING) {
      Scalar coef1 = (1.0 - 1.0 / n_idx) * w;
      Scalar coef2 = -w / n_idx;
      for (int i = 0; i < n_idx; ++i) {
        for (int j = 0; j < n_idx; ++j) {
          triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
        }

        idO++;
      }
    } else if (transform_type_ == SUBTRACT_FIRST) {
      for (int i = 1; i < n_idx; ++i) {
        triplets.push_back(Triplet(idO, idI_[0], -w));
        triplets.push_back(Triplet(idO, idI_[i], w));
        idO++;
      }
    } else {
      for (int i = 0; i < n_idx; ++i) {
        triplets.push_back(Triplet(idO++, idI_[i], w));
      }
    }
  }

  int num_indices() const {
    return idI_.size();
  }

  int num_transformed_points() const {
    if (transform_type_ == SUBTRACT_FIRST) {
      return num_indices() - 1;
    } else {
      return num_indices();
    }
  }

 protected:

  Eigen::VectorXi idI_;

  Scalar weight_;

  InvarianceTransformType transform_type_;

  int idO_;

  // Implementation of the projection, each subclass should override this method
  virtual void project_impl(const Eigen::Map<const MatrixNX>& transformed_pts,
                            Eigen::Map<MatrixNX>& projection) {
    projection = transformed_pts;
  }

  Scalar clamp(Scalar val, Scalar min_value, Scalar max_value) {
    return std::min(std::max(min_value, val), max_value);
  }
};

template<unsigned int N>
class EdgeLengthConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;
  using typename Constraint<N>::MatrixNX;

 public:
  EdgeLengthConstraint(int idx1, int idx2, Scalar weight, Scalar target_length)
      : Constraint<N>(std::vector<int>( { idx1, idx2 }), weight,
                      Constraint<N>::SUBTRACT_FIRST),
        target_length_(target_length) {
  }

  virtual ~EdgeLengthConstraint() {
  }

 protected:
  virtual void project_impl(const Eigen::Map<const MatrixNX>& transformed_pts,
                            Eigen::Map<MatrixNX>& projection) {
    projection.col(0) = transformed_pts.col(0).normalized() * target_length_;
  }

 private:
  Scalar target_length_;
};

template<unsigned int N>
class AngleConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;
  using typename Constraint<N>::MatrixNX;
  using Constraint<N>::clamp;

 public:
  AngleConstraint(int tip_idx, int side_idx1, int side_idx2, Scalar weight,
                  Scalar min_radian, Scalar max_radian)
      : Constraint<N>(std::vector<int>( { tip_idx, side_idx1, side_idx2 }),
                      weight, Constraint<N>::SUBTRACT_FIRST),
        min_angle_(std::max(Scalar(0), min_radian)),
        max_angle_(std::min(Scalar(M_PI), max_radian)) {
    min_angle_cos_ = clamp(std::cos(min_angle_), Scalar(-1), Scalar(1));
    max_angle_cos_ = clamp(std::cos(max_angle_), Scalar(-1), Scalar(1));
    assert(max_angle_cos_ <= min_angle_cos_);
  }

  virtual ~AngleConstraint() {
  }

 protected:
  virtual void project_impl(const Eigen::Map<const MatrixNX>& transformed_pts,
                            Eigen::Map<MatrixNX>& projection) {
    projection = transformed_pts;

    VectorN v1 = transformed_pts.col(0), v2 = transformed_pts.col(1);
    Scalar epsilon = 1e-14;
    Scalar v1_sqrnorm = v1.squaredNorm(), v2_sqrnorm = v2.squaredNorm();
    Scalar v1_norm = v1.norm(), v2_norm = v2.norm();
    VectorN unit_v1 = v1.normalized(), unit_v2 = v2.normalized();

    // cosine value of the angle gamma between v1 and v2
    Scalar cos_gamma = clamp(unit_v1.dot(unit_v2), Scalar(-1), Scalar(1));

    // Proceed only when the current angle lies outside the target range, and v1, v2 are not colinear
    if ((Scalar(1) - std::abs(cos_gamma) > epsilon)
        && (cos_gamma > min_angle_cos_ || cos_gamma < max_angle_cos_)) {
      Scalar gamma = std::acos(cos_gamma);

      // Angle eta: the sum of displacement angles from v1 and v2
      Scalar eta =
          cos_gamma > min_angle_cos_ ?
              (min_angle_ - gamma) : (gamma - max_angle_);
      eta = (std::max)(eta, 0.0);

      // Compute angle theta between v1 and its projection, and angle phi between v2 and its projection
      Scalar theta = Scalar(0.5)
          * std::atan2(v2_sqrnorm * std::sin(2 * eta),
                       v1_sqrnorm + v2_sqrnorm * std::cos(2 * eta));
      theta = (std::max)(0.0, (std::min)(eta, theta));
      Scalar phi = eta - theta;

      // Compute unit vectors that are coplanar with v1, v2, and orthogonal to one of them.
      // They form orthogonal frames with v1 and v2 respectively, within which we compute the projection using the above angles
      VectorN unit_v3 = (unit_v2 - unit_v1 * cos_gamma).normalized(), unit_v4 =
          (unit_v1 - unit_v2 * cos_gamma).normalized();

      // Determine if v1, v2 should move away from each other or towards each other
      if (cos_gamma > min_angle_cos_) {
        unit_v3 *= -1.0;
        unit_v4 *= -1.0;
      }

      projection.col(0) =
          (unit_v1 * std::cos(theta) + unit_v3 * std::sin(theta))
              * (v1_norm * std::cos(theta));
      projection.col(1) = (unit_v2 * std::cos(phi) + unit_v4 * std::sin(phi))
          * (v2_norm * std::cos(phi));
    }
  }

 private:
  Scalar min_angle_, max_angle_;
  Scalar min_angle_cos_, max_angle_cos_;
};

// Use this class for dynamic handle constraint: handle position is updated using set_target_pos
template<unsigned int N>
class ClosenessConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;

 public:
  ClosenessConstraint(int idx, Scalar weight, const VectorN &target_pos)
      : Constraint<N>(std::vector<int>( { idx }), weight,
                      Constraint<N>::IDENTITY),
        target_pos_(target_pos) {
  }

  virtual ~ClosenessConstraint() {
  }

  void set_target_pos(const VectorN &target_pos) {
    target_pos_ = target_pos;
  }

 protected:
  virtual void proj_impl(const Eigen::Map<const Matrix3X>&,
                         Eigen::Map<Matrix3X>& projection) {
    projection.col(0) = target_pos_;
  }

 private:
  VectorN target_pos_;
};

class PointToRefSurfaceConstraint : public Constraint<3> {
 public:
  PointToRefSurfaceConstraint(int pt_idx, Scalar weight,
                              const std::shared_ptr<TriMeshAABB> &aabb)
      : Constraint(std::vector<int>( { pt_idx }), weight, IDENTITY),
        AABB_tree_(aabb) {
  }

  virtual ~PointToRefSurfaceConstraint() {
  }

 protected:
  virtual void project_impl(const Eigen::Map<const MatrixNX>& transformed_pts,
                            Eigen::Map<MatrixNX>& projection) {
    Vector3 closest_point;
    AABB_tree_->get_closest_point(transformed_pts.col(0), closest_point);
    projection.col(0) = closest_point;
  }

 private:
  std::shared_ptr<TriMeshAABB> AABB_tree_;
};

class ReferenceSurfceConstraint : public Constraint<3> {
 public:
  ReferenceSurfceConstraint(int n_points, Scalar weight,
                            const Matrix3X &ref_surface_vtx,
                            const Eigen::Matrix3Xi &ref_surface_faces)
      : Constraint(std::vector<int>(), weight, IDENTITY) {
    n_points_ = n_points;
    idI_.resize(n_points);

    for (int i = 0; i < n_points; ++i) {
      idI_(i) = i;
    }

    sqrD_.setZero(n_points);
    I_.setZero(n_points);
    C_.setZero(n_points, 3);

    ref_surface_vtx_ = ref_surface_vtx.transpose();
    ref_surface_faces_ = ref_surface_faces.transpose();
    AABB_tree_.init(ref_surface_vtx_, ref_surface_faces_);
  }

  virtual ~ReferenceSurfceConstraint() {
  }

  // Perform projection, and return the weighted squared distance from the transformed points to the projection points
  virtual void project_impl(const Eigen::Map<const Matrix3X>& transformed_pts,
                            Eigen::Map<Matrix3X>& projection) {
    pos_ = transformed_pts.transpose();
    AABB_tree_.squared_distance(ref_surface_vtx_, ref_surface_faces_, pos_,
                                sqrD_, I_, C_);
    projection = C_.transpose();
  }

 private:
  int n_points_;
  MatrixX3 pos_;
  MatrixX3 ref_surface_vtx_;
  Eigen::MatrixX3i ref_surface_faces_;
  VectorX sqrD_;
  Eigen::VectorXi I_;
  MatrixX3 C_;
  igl::AABB<MatrixX3, 3> AABB_tree_;
};

class PlaneConstraint : public Constraint<3> {
 public:
  PlaneConstraint(const std::vector<int> &idI, Scalar weight)
      : Constraint(idI, weight, MEAN_CENTERING) {
  }

  virtual ~PlaneConstraint() {
  }

 protected:
  virtual void project_impl(const Eigen::Map<const MatrixNX>& transformed_pts,
                            Eigen::Map<MatrixNX>& projection) {
    Eigen::JacobiSVD<Matrix3X, Eigen::FullPivHouseholderQRPreconditioner> jSVD;
    jSVD.compute(transformed_pts, Eigen::ComputeFullU);
    Vector3 best_fit_normal = jSVD.matrixU().col(2).normalized();
    projection = transformed_pts
        - best_fit_normal * (best_fit_normal.transpose() * transformed_pts);
  }
};

#endif // CONSTRAINT_H
