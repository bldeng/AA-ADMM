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

#ifndef ADMM_TETENERGYTERM_H
#define ADMM_TETENERGYTERM_H

#include "EnergyTerm.hpp"
#include "MCL/Newton.hpp"
#include "MCL/LBFGS.hpp"
#include "MCL/NonLinearCG.hpp"

namespace admm {


//
//	Creates a bunch of tet energy terms from a tet mesh
//
template <typename IN_SCALAR, typename TYPE>
inline void create_tets_from_mesh( std::vector< std::shared_ptr<EnergyTerm> > &energyterms,
    const IN_SCALAR *verts, const int *inds, int n_tets, const Lame &lame, const int vertex_offset ){
    typedef Eigen::Matrix<int,4,1> Vec4i;
    typedef Eigen::Matrix<double,3,1> Vec3;
    for( int i=0; i<n_tets; ++i ){
        Vec4i tet(inds[i*4+0], inds[i*4+1], inds[i*4+2], inds[i*4+3]);
        std::vector<Vec3> tetverts = {
            Vec3( verts[tet[0]*3+0], verts[tet[0]*3+1], verts[tet[0]*3+2] ),
            Vec3( verts[tet[1]*3+0], verts[tet[1]*3+1], verts[tet[1]*3+2] ),
            Vec3( verts[tet[2]*3+0], verts[tet[2]*3+1], verts[tet[2]*3+2] ),
            Vec3( verts[tet[3]*3+0], verts[tet[3]*3+1], verts[tet[3]*3+2] )
        };
        tet += Vec4i(1,1,1,1)*vertex_offset;
        energyterms.emplace_back( std::make_shared<TYPE>(tet, tetverts, lame) );
    }
} // end create from mesh


//
//	The tet energy term base class
//
class TetEnergyTerm : public EnergyTerm {
protected:
    typedef Eigen::Matrix<int,4,1> Vec4i;
    typedef Eigen::Matrix<double,3,1> Vec3;
    typedef Eigen::Matrix<double,9,1> Vec9;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatX;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
    Vec4i tet;
    Lame lame;
    double volume;
    double weight;
    VecX zi_;
    Eigen::Matrix3d edges_inv;

public:
    int get_dim() const { return 9; }
    double get_weight() const { return weight; }
    double get_volume() const { return volume; }

    TetEnergyTerm( const Vec4i &tet_, const std::vector<Vec3> &verts, const Lame &lame_ );

    void get_reduction( std::vector< Eigen::Triplet<double> > &triplets );

    // Unless derived from uses linear tet strain (no volume conservation)
    virtual void prox(const MatX &W, VecX &zi , const VecX &vi);
    virtual double prox_for_strain_limiting_energy( VecX &zi );
    virtual double energy( const VecX &F );
    virtual double energyLBFGS( const VecX &F );
    //virtual double gradient( const VecX &F, VecX &grad );
    virtual void gradient(const VecX &x, VecX &grad);
    virtual void get_gradient( const VecX &F, VecX &grad );

}; // end class TetEnergyTerm


//
//	HyperElastic Tet base class
//
class HyperElasticTet : public TetEnergyTerm {
public:
    typedef Eigen::Matrix<double,3,3> Mat3;
    typedef Eigen::Matrix<double,9,9> Mat9;
    class Prox : public mcl::optlib::Problem<double,9> {
    public:
        virtual void set_x0(const Vec3 x0_)=0; // set x0 (quad penalty)
        virtual bool converged(const Vec9 &x0, const Vec9 &x1, const Vec9 &grad){
            return ( grad.norm() < 1e-10 || (x0-x1).norm() < 1e-10 );
        }
        virtual void set_W(const Mat9 W_)=0;
        virtual void set_v(const Vec9 vi_)=0;
        virtual void set_vol(const double vol_)=0;

        virtual void U_gradient(const Mat3 &F, Mat3 &grad)=0;
    };
    mcl::optlib::LBFGS<double,9> solver;

    // Returns a pointer to the local problem (constitutive model)
    virtual Prox* get_problem()=0;

    // Sets some default values for weight
    HyperElasticTet( const Vec4i &tet_, const std::vector<Vec3> &verts, const Lame &lame_ ) :
        TetEnergyTerm(tet_,verts,lame_) {}

    // Local step on the problem returned from get_problem
    void prox(const MatX &W, VecX &zi, const VecX &vi);
    double prox_for_strain_limiting_energy( VecX &zi );
    double energy( const VecX &F );
    double energyLBFGS( const VecX &F );
    //double gradient( const VecX &F, VecX &grad );
    void gradient(const VecX &x, VecX &grad);
    void get_gradient( const VecX &F, VecX &grad );
};


//
//	NeoHookean Tet
//
class NeoHookeanTet : public HyperElasticTet {
public:
    class NHProx : public HyperElasticTet::Prox {
    public:
        double mu, lambda, k, vol;
        Vec3 x0;
        Mat9 W;
        Vec9 vi;
        void set_lame(const Lame &lame){
            mu = lame.mu;
            lambda = lame.lambda;
            k = lame.bulk_modulus();
        }
        void set_x0(const Vec3 x0_){ x0=x0_; }
        void set_W(const Mat9 W_){ W=W_; }
        void set_v(const Vec9 vi_){ vi=vi_; }
        void set_vol(const double vol_){ vol=vol_; }
        double energy_density(const Vec9 &x) const;
        double value(const Vec9 &x);
        double valueLBFGS(const Vec9 &x);
        double gradient(const Vec9 &x, Vec9 &grad);
        void U_gradient(const Mat3 &F, Mat3 &grad);
    } problem;

    Prox* get_problem(){ return &problem; }
    NeoHookeanTet( const Vec4i &tet_, const std::vector<Vec3> &verts, const Lame &lame_ ) :
        HyperElasticTet(tet_,verts,lame_) { problem.set_lame(lame_); }
};


//
//	St-VK Tet
//
class StVKTet : public HyperElasticTet {
public:
    class StVKProx : public HyperElasticTet::Prox {
    public:
        double mu, lambda, k, vol;
        Vec3 x0;
        Mat9 W;
        Vec9 vi;
        void set_lame(const Lame &lame){
            mu = lame.mu;
            lambda = lame.lambda;
            k = lame.bulk_modulus();
        }
        void set_x0(const Vec3 x0_){ x0=x0_; }
        void set_W(const Mat9 W_){ W=W_; }
        void set_v(const Vec9 vi_){ vi=vi_; }
        void set_vol(const double vol_){ vol=vol_; }
        double energy_density(const Vec9 &x) const;
        double value(const Vec9 &x);
        double valueLBFGS(const Vec9 &x);
        double gradient(const Vec9 &x, Vec9 &grad);
        void U_gradient(const Mat3 &F, Mat3 &grad);
        static inline double ddot( Vec3 &a, Vec3 &b ) { return (a*b.transpose()).trace(); }
        static inline double v3trace( const Vec3 &v ) { return v[0]+v[1]+v[2]; }
    } problem;

    Prox* get_problem(){ return &problem; }
    StVKTet( const Vec4i &tet_, const std::vector<Vec3> &verts, const Lame &lame_ ) :
        HyperElasticTet(tet_,verts,lame_) { problem.set_lame(lame_); }
};

} // end namespace admm

#endif




