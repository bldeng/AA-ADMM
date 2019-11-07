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
//#include "OMPHelper.h"

class AndersonAcceleration
{
public:
    AndersonAcceleration()
    :m_(-1), dim_(-1), iter_(-1), col_idx_(-1) {}

//    void replace(const VectorX &g, const VectorX &d)
//    {
//        current_u_.head(dim_head_)=g;
//        current_u_.tail(dim_tail_)=d;
//    }

    void replace(const VectorX &g)
    {
        current_u_=g;
    }

//    void compute(VectorX &g, VectorX &d)
//    {
//        assert(iter_ >= 0);

//        //Timer timer;
//        //Timer::EventID before_energy = timer.get_time();
//        VectorX G(dim_);
//        G.head(dim_head_)=g;
//        G.tail(dim_tail_)=d;
//        //Timer::EventID after_energy = timer.get_time();
//        //std::cout << "Anderson time 1: " << timer.elapsed_time(before_energy, after_energy) << " seconds" << std::endl;

//        current_F_ = G - current_u_;

//        if(iter_ == 0)
//        {
//            prev_dF_.col(0) = -current_F_;
//            prev_dG_.col(0) = -G;
//            current_u_ = G;
//        }
//        else
//        {
//            //before_energy = timer.get_time();
//            prev_dF_.col(col_idx_) += current_F_;
//            prev_dG_.col(col_idx_) += G;

//            Scalar eps = 1e-14;
//            Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
//            dF_scale_(col_idx_) = scale;
//            prev_dF_.col(col_idx_) /= scale;
//            //after_energy = timer.get_time();
//            //std::cout << "Anderson time scale: " << timer.elapsed_time(before_energy, after_energy) << " seconds" << std::endl;


//            int m_k = std::min(m_, iter_);

//            if(m_k == 1)
//            {
//                theta_(0) = 0;
//                Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
//                M_(0, 0) = dF_sqrnorm;
//                Scalar dF_norm = std::sqrt(dF_sqrnorm);

//                if(dF_norm > eps){
//                    // compute theta = (dF * F) / (dF * dF)
//                    theta_(0) = (prev_dF_.col(col_idx_)/dF_norm).dot(current_F_/dF_norm);
//                }
//            }
//            else
//            {
//                // Update the normal equation matrix, for the column and row corresponding to the new dF column
//                VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose() * prev_dF_.block(0, 0, dim_, m_k)).transpose();
//                M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
//                M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

//                // Solve normal equation
//                cod_.compute(M_.block(0, 0, m_k, m_k));
//                theta_.head(m_k) = cod_.solve(prev_dF_.block(0, 0, dim_, m_k).transpose() * current_F_);
//            }

//            //before_energy = timer.get_time();
//            // Use rescaled theata to compute new u
//            // prev_dF_ have already been scaled
//            current_u_ = (1 - beta_) * (current_u_  - prev_dF_.block(0, 0, dim_, m_k) * ((theta_.head(m_k).array()).matrix()))
//                    + beta_ * (G  - prev_dG_.block(0, 0, dim_, m_k) * ((theta_.head(m_k).array() / dF_scale_.head(m_k).array()).matrix()));
//            //after_energy = timer.get_time();
//            //std::cout << "Anderson time current_u: " << timer.elapsed_time(before_energy, after_energy) << " seconds" << std::endl;

//            col_idx_ = (col_idx_ + 1) % m_;
//            prev_dF_.col(col_idx_) = -current_F_;
//            prev_dG_.col(col_idx_) = -G;
//        }

//        iter_++;
//        //before_energy = timer.get_time();
//        g = current_u_.head(dim_head_);
//        d = current_u_.tail(dim_tail_);
//        //after_energy = timer.get_time();
//        //std::cout << "Anderson time 2: " << timer.elapsed_time(before_energy, after_energy) << " seconds" << std::endl;
//        //return current_u_;
//    }

    void compute(VectorX &curr_g, const VectorX &g)
    {
        assert(iter_ >= 0);

        VectorX G(dim_);
        G = g;

        current_F_ = G - current_u_;

        if(iter_ == 0)
        {
            prev_dF_.col(0) = -current_F_;
            prev_dG_.col(0) = -G;
            current_u_ = G;
        }
        else
        {
            prev_dF_.col(col_idx_) += current_F_;
            prev_dG_.col(col_idx_) += G;

            Scalar eps = 1e-14;
            Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
            dF_scale_(col_idx_) = scale;
            prev_dF_.col(col_idx_) /= scale;


            int m_k = std::min(m_, iter_);

            if(m_k == 1)
            {
                theta_(0) = 0;
                Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
                M_(0, 0) = dF_sqrnorm;
                Scalar dF_norm = std::sqrt(dF_sqrnorm);

                if(dF_norm > eps){
                    // compute theta = (dF * F) / (dF * dF)
                    theta_(0) = (prev_dF_.col(col_idx_)/dF_norm).dot(current_F_/dF_norm);
                }
            }
            else
            {
                // Update the normal equation matrix, for the column and row corresponding to the new dF column
                VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose() * prev_dF_.block(0, 0, dim_, m_k)).transpose();
                M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
                M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

                // Solve normal equation
                cod_.compute(M_.block(0, 0, m_k, m_k));
                theta_.head(m_k) = cod_.solve(prev_dF_.block(0, 0, dim_, m_k).transpose() * current_F_);
            }

            current_u_ = G  - prev_dG_.block(0, 0, dim_, m_k) * ((theta_.head(m_k).array() / dF_scale_.head(m_k).array()).matrix());

            col_idx_ = (col_idx_ + 1) % m_;
            prev_dF_.col(col_idx_) = -current_F_;
            prev_dG_.col(col_idx_) = -G;
        }

        iter_++;

        curr_g = current_u_;
    }


//    void compute_type_one(VectorX &g, VectorX &d)
//    {
//        assert(iter_ >= 0);

//        VectorX G(dim_);
//        G.head(dim_head_)=g;
//        G.tail(dim_tail_)=d;

//        current_F_ = G - current_u_;

//        if(iter_ == 0)
//        {
//            prev_dF_.col(0) = -current_F_;
//            prev_dG_.col(0) = -G;
//            current_u_ = G;
//        }
//        else
//        {
//            //before_energy = timer.get_time();
//            prev_dF_.col(col_idx_) += current_F_;
//            prev_dG_.col(col_idx_) += G;

//            Scalar eps = 1e-14;
//            Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
//            Scalar scale_G = std::max(eps, prev_dG_.col(col_idx_).norm());
//            dF_scale_(col_idx_) = scale;
//            prev_dF_.col(col_idx_) /= scale;
//            dG_scale_(col_idx_) = scale_G;
//            prev_dG_.col(col_idx_) /= scale_G;


//            int m_k = std::min(m_, iter_);

//            if(m_k == 1)
//            {
//                theta_(0) = 0;
//                Scalar dF_norm = prev_dF_.col(col_idx_).norm();
//                Scalar dG_norm = prev_dG_.col(col_idx_).norm();
//                M_(0, 0) = dF_norm * dG_norm;

//                if((dF_norm > eps) && (dG_norm > eps)){
//                    // compute theta = (dG * F) / (dG * dF)
//                    theta_(0) = (prev_dG_.col(col_idx_)/dG_norm).dot(current_F_/dF_norm);
//                }
//            }
//            else
//            {
//                // Update the normal equation matrix, for the column and row corresponding to the new dF column
//                M_.block(col_idx_, 0, 1, m_k) = prev_dG_.col(col_idx_).transpose() * prev_dF_.block(0, 0, dim_, m_k);
//                M_.block(0, col_idx_, m_k, 1) = (prev_dF_.col(col_idx_).transpose() * prev_dG_.block(0, 0, dim_, m_k)).transpose();

//                // Solve normal equation
//                cod_.compute(M_.block(0, 0, m_k, m_k));
//                theta_.head(m_k) = cod_.solve(prev_dG_.block(0, 0, dim_, m_k).transpose() * current_F_);
//            }

//            // Use rescaled theata to compute new u
//            // prev_dF_ have already been scaled
//            current_u_ = (1 - beta_) * (current_u_  - prev_dF_.block(0, 0, dim_, m_k) * ((theta_.head(m_k).array()).matrix()))
//                    + beta_ * (G  - prev_dG_.block(0, 0, dim_, m_k) * ((dG_scale_.head(m_k).array() * theta_.head(m_k).array() / dF_scale_.head(m_k).array()).matrix()));

//            col_idx_ = (col_idx_ + 1) % m_;
//            prev_dF_.col(col_idx_) = -current_F_;
//            prev_dG_.col(col_idx_) = -G;
//        }

//        iter_++;

//        g = current_u_.head(dim_head_);
//        d = current_u_.tail(dim_tail_);

//    }

    // m: number of previous iterations used
    // d: dimension of variables
    // u0: initial variable values
    void init(int m, int d, const VectorX &g0)
    {
        assert(m > 0);
        m_ = m;
        dim_ = d;
        current_u_.resize(d);
        current_F_.resize(d);
        prev_dG_.resize(d, m);
        prev_dF_.resize(d, m);
        M_.resize(m, m);
        theta_.resize(m);
        dF_scale_.resize(m);
        dG_scale_.resize(m);
        current_u_=g0;
        iter_ = 0;
        col_idx_ = 0;
    }

    //    void init(int m, int d, int d_head, double beta, const VectorX &g0, const VectorX &d0)
    //    {
    //        assert(m > 0);
    //        m_ = m;
    //        dim_ = d;
    //        dim_head_ = d_head;
    //        dim_tail_ = d - d_head;
    //        beta_ = beta;
    //        current_u_.resize(d);
    //        current_F_.resize(d);
    //        prev_dG_.resize(d, m);
    //        prev_dF_.resize(d, m);
    //        M_.resize(m, m);
    //        theta_.resize(m);
    //        dF_scale_.resize(m);
    //        dG_scale_.resize(m);
    //        current_u_.head(dim_head_)=g0;
    //        current_u_.tail(dim_tail_)=d0;
    //        iter_ = 0;
    //        col_idx_ = 0;
    //    }

private:
    VectorX current_u_;
    VectorX current_x_;
    VectorX current_F_;
    MatrixXX prev_dG_;
    MatrixXX prev_dF_;
    MatrixXX M_;		// Normal equations matrix for the computing theta
    VectorX	theta_;	// theta value computed from normal equations
    VectorX dF_scale_, dG_scale_;		// The scaling factor for each column of prev_dF
    Eigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;

    int m_;		// Number of previous iterates used for Andreson Acceleration
    int dim_;	// Dimension of variables
    int iter_;	// Iteration count since initialization
    int col_idx_;	// Index for history matrix column to store the next value
//    int dim_head_; // Dimension of variables G
//    int dim_tail_; // Dimension of variables D
};


#endif /* ANDERSONACCELERATION_H_ */
