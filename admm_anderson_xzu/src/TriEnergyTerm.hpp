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

#ifndef ADMM_TRIENERGYTERM_H
#define ADMM_TRIENERGYTERM_H 1

#include "EnergyTerm.hpp"

namespace admm {


//
//	Creates a bunch of energy terms from a mesh
//
template <typename IN_SCALAR, typename TYPE>
inline void create_tris_from_mesh( std::vector< std::shared_ptr<EnergyTerm> > &energyterms,
    const IN_SCALAR *verts, const int *inds, int n_tris, const Lame &lame, const int vertex_offset ){
    typedef Eigen::Matrix<int,3,1> Vec3i;
    typedef Eigen::Matrix<double,3,1> Vec3;
    for( int i=0; i<n_tris; ++i ){
        Vec3i tri(inds[i*3+0], inds[i*3+1], inds[i*3+2]);
        std::vector<Vec3> triverts = {
            Vec3( verts[tri[0]*3+0], verts[tri[0]*3+1], verts[tri[0]*3+2] ),
            Vec3( verts[tri[1]*3+0], verts[tri[1]*3+1], verts[tri[1]*3+2] ),
            Vec3( verts[tri[2]*3+0], verts[tri[2]*3+1], verts[tri[2]*3+2] )
        };
        tri += Vec3i(1,1,1)*vertex_offset;
        energyterms.emplace_back( std::make_shared<TYPE>( tri, triverts, lame ) );
    }
} // end create from mesh


//
//	The tri energy term base class
//
class TriEnergyTerm : public EnergyTerm {
protected:
    typedef Eigen::Matrix<int,3,1> Vec3i;
    typedef Eigen::Matrix<double,3,1> Vec3;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatX;
    Vec3i tri;
    Lame lame;

    double area;
    double weight;
    Eigen::Matrix2d rest_pose;

public:
    int get_dim() const { return 6; }
    double get_weight() const { return weight; }
    double get_volume() const { return area; }

    TriEnergyTerm( const Vec3i &tri_, const std::vector<Vec3> &verts, const Lame &lame_ );

    void get_reduction( std::vector< Eigen::Triplet<double> > &triplets );
    void get_gradient_for_vertices( const VecX &grad, VecX &all_grad);

    // Unless derived from uses linear strain (no area conservation)
    virtual void prox(const MatX &W, VecX &zi, const VecX &vi );
    virtual double energy( const VecX &F );
    virtual double energyLBFGS( const VecX &F );
    virtual void gradient( const VecX &F, VecX &grad );
    virtual void get_gradient( const VecX &F, VecX &grad );

}; // end class TriEnergyTerm


} // ns admm
#endif // trienergyterm
