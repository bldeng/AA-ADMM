//  BSD 3-Clause License
//
//  Copyright (c) 2019, Yue Peng
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

#ifndef ADMM_COLLISIONENERGYTERM_H
#define ADMM_COLLISIONENERGYTERM_H 1

#include "EnergyTerm.hpp"
#include "ConstraintSet.hpp"

namespace admm {
//
//	It's an infinitely hard hard spring, with energy = inf when violated, zero otherwise.
//
class Collision : public EnergyTerm {
protected:
    typedef Eigen::Matrix<double,3,1> Vec3;
    typedef Eigen::Matrix<int,3,1> Vec3i;
    typedef Eigen::Matrix<double,4,1> Vec4;
    typedef Eigen::Matrix<int,4,1> Vec4i;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatX;
    typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMat;
    int idx; // constrained vertex
    std::shared_ptr<ConstraintSet> p_;
    bool active;
    double weight;
    double volume;

public:
    int get_dim() const { return 3; }
    double get_weight() const { return weight; }
    double get_volume() const { return volume; }
    //void set_pin( const Vec3 &p ){ pin = p; }
    void set_active( bool a ){ active = a; }

    Collision( int idx_, const std::shared_ptr<ConstraintSet> p) : idx(idx_), p_(p), active(true) {
        // Because we usually use bulk mod of rubber for elastics,
        // We'll make a really strong rubber and use that for pin.
        admm::Lame lame = admm::Lame::soft_rubber();
        weight = std::sqrt(lame.bulk_modulus()*2.0);
        volume = 2.0;
        //std::cout << "weight = " << weight << std::endl;
    }

    void get_reduction( std::vector< Eigen::Triplet<double> > &triplets ){
        const int col = 3*idx;
        triplets.emplace_back( 0, col+0, 1.0 );
        triplets.emplace_back( 1, col+1, 1.0 );
        triplets.emplace_back( 2, col+2, 1.0 );
    }

    void prox( const MatX &W, VecX &zi, const VecX &vi ){
        (void)(W);
        (void)(vi);
        //detect collision
        PassiveCollision::Payload p_payload(idx);
        for( size_t j=0; j<p_->collider->passive_objs.size(); ++j ){
            p_->collider->passive_objs[j]->signed_distance(zi, p_payload);
        }
        if( p_payload.dx < 0 ){
            zi = p_payload.point;
        }

    }

    double prox_for_strain_limiting_energy( VecX &zi ){
        (void)(zi);
        throw std::runtime_error("**Collision Error: Energy not implemented");
    }

    double energy( const VecX &F ){
        (void)(F); // could use hookean spring energy, but not really accurate
        throw std::runtime_error("**Collision Error: Energy not implemented");
    }
    double energyLBFGS( const VecX &F ){
        (void)(F); // could use hookean spring energy, but not really accurate
        throw std::runtime_error("**Collision Error: Energy not implemented");
    }

    void gradient( const VecX &x, VecX &grad ){
        (void)(x); (void)(grad);
        throw std::runtime_error("**Collision Error: No gradient for hard constraint");
    }

    void get_gradient( const VecX &F, VecX &grad ){
        (void)(F); (void)(grad);
        throw std::runtime_error("**Collision Error: No gradient for hard constraint");
    }

}; // end class Collision

} // end namespace admm

#endif




