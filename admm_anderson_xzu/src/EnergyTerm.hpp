// Original work Copyright (c) 2016, University of Minnesota
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

#ifndef ADMM_FORCE_H
#define ADMM_FORCE_H 1

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <memory>
#include <iostream>

namespace admm {


//
//	Lame constants
//
class Lame {
public:
    static Lame rubber(){ return Lame(10000000,0.499); } // true rubber
    static Lame soft_rubber(){ return Lame(10000000,0.399); } // fun rubber!
    static Lame very_soft_rubber(){ return Lame(1000000,0.299); } // more funner!

    double mu, lambda;
    double bulk_modulus() const { return lambda + (2.0/3.0)*mu; }

    // Hard strain limiting (e.g. [0.95,1.05]), default no limit
    // with  min: -inf to 1, max: 1 to inf.
    // In practice if max>99 it's basically no limiting.
    double limit_min, limit_max;

    // k: Youngs (Pa), measure of stretch
    // v: Poisson, measure of incompressibility
    Lame( double k, double v ) :
        mu(k/(2.0*(1.0+v))),
        lambda(k*v/((1.0+v)*(1.0-2.0*v))),
        limit_min(-100.0),
        limit_max(100.0) {
    }

    // Use custom mu, lambda
    Lame(): limit_min(-100.0),
        limit_max(100.0) {}
};


//
//	Base class tets: Linear (non-corotated) elastic
//
class EnergyTerm {
private:
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatX;
    typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMat;
    int g_index; // global idx (starting row of reduction matrix)

public:

    virtual ~EnergyTerm() {}

    // Called by the solver to create the global reduction and weight matrices
    inline void get_reduction( std::vector< Eigen::Triplet<double> > &triplets, std::vector<double> &weights );

    // Called by the solver for a local step update
    //inline void update( const SparseMat &D, const VecX &x, VecX &z, VecX &u );

    // Called by the solver for a local step update
    inline void update_z(const SparseMat &D, const SparseMat &W_inv, const SparseMat &W, const VecX &x, VecX &z, const VecX &u , const VecX &c);

    // Called by the solver for a local step update
    inline void update_u(const SparseMat &D, const SparseMat &W, const VecX &x, const VecX &z, VecX &u , const VecX &c);

    // Computes energy of the force (used for debugging)
    inline double energy( const SparseMat &D, const VecX &x );

    inline double get_all_energy(const VecX &x );
    // Compute energy and gradient for a first-order opt. solver
//	inline double gradient( const SparseMat &D, const VecX &x, VecX &grad );
    inline void gradient( const SparseMat &D, const VecX &x, VecX &grad );

    inline void get_all_gradient(const VecX &z, VecX &grad );

    inline void compute_x_gradient(const VecX &grad, VecX &all_grad );

    // Dimension of deformation gradient.
    virtual int get_dim() const = 0;

    // Return the scalar weight of the energy term
    virtual double get_weight() const = 0;

    virtual double get_volume() const = 0;

protected:

    virtual void get_gradient_for_vertices( const VecX &grad, VecX &all_grad) = 0;

    // Get a local reduction matrix
    virtual void get_reduction( std::vector< Eigen::Triplet<double> > &triplets ) = 0;

    // Proximal update
    virtual void prox( const MatX &W, VecX &zi, const VecX &vi) = 0;

    // Returns energy of the force
    virtual double energy( const VecX &F ) = 0;

    virtual double energyLBFGS( const VecX &F ) = 0;

    // Computes a first-order update (energy and gradient)
    virtual void gradient(const VecX &x, VecX &grad) = 0;

    virtual void get_gradient( const VecX &F, VecX &grad ) = 0;

}; // end class EnergyTerm

//
//  Implementation
//

inline void EnergyTerm::get_reduction( std::vector< Eigen::Triplet<double> > &triplets, std::vector<double> &weights ){
    std::vector< Eigen::Triplet<double> > temp_triplets;

    get_reduction( temp_triplets );
    int n_trips = temp_triplets.size();
    g_index = weights.size();

    for( int i=0; i<n_trips; ++i ){
        const Eigen::Triplet<double> &trip = temp_triplets[i];
        triplets.emplace_back( trip.row()+g_index, trip.col(), trip.value() );
    }

    int dim = get_dim();
    double w = get_weight();
    if( w <= 0.0 ){
        throw std::runtime_error("**EnergyTerm::get_reduction Error: Some weight leq 0");
    }
    for( int i=0; i<dim; ++i ){ weights.emplace_back( w ); }
}

inline void EnergyTerm::update_u( const SparseMat &D, const SparseMat &W, const VecX &x, const VecX &z, VecX &u, const VecX &c ){
    int dof = x.rows();
    int dim = get_dim();
    VecX Dix = D.block(g_index,0,dim,dof)*x;
    VecX ui = u.segment(g_index,dim);
    VecX zi = W.block(g_index,g_index,dim,dim) * z.segment(g_index,dim);
    VecX ci = c.segment(g_index,dim);
    ui += (Dix - zi - ci);
    u.segment(g_index,dim) = ui;
}

inline void EnergyTerm::update_z( const SparseMat &D, const SparseMat &W_inv, const SparseMat &W, const VecX &x, VecX &z, const VecX &u, const VecX &c ){
    int dof = x.rows();
    const int dim = get_dim();
    VecX Dix = D.block(g_index,0,dim,dof)*x;
    VecX ui = u.segment(g_index,dim);
    VecX ci = c.segment(g_index,dim);
    VecX vi = Dix + ui - ci;
    VecX zi = W_inv.block(g_index,g_index,dim,dim) * vi;
    vi = zi;
    Eigen::Matrix<double,9,9> Wi = W.block(g_index,g_index,dim,dim);
    prox(Wi, zi, vi);
    z.segment(g_index,dim) = zi;
}

inline double EnergyTerm::energy( const SparseMat &D, const VecX &x ){
    (void)(D);
    int dim = get_dim();
    VecX Dix = x.segment(g_index,dim);
    return energy( Dix );
}

inline double EnergyTerm::get_all_energy( const VecX &x ){
    int dim = get_dim();
    VecX Dix = x.segment(g_index,dim);
    return energyLBFGS( Dix );
}

inline void EnergyTerm::gradient( const SparseMat &D, const VecX &x, VecX &grad ){
    int dof = x.rows();
    int dim = get_dim();
    VecX Dix = D.block(g_index,0,dim,dof)*x;
    gradient( Dix, grad );
}

inline void EnergyTerm::get_all_gradient( const VecX &z, VecX &grad ){
    int dim = get_dim();
    VecX zi = z.segment(g_index,dim);
    VecX grad_z = grad.segment(g_index,dim);
    get_gradient( zi, grad_z );
    grad.segment(g_index,dim) = grad_z;
}

inline void EnergyTerm::compute_x_gradient(const VecX &grad, VecX &all_grad ){
    int dim = get_dim();
    VecX gi = grad.segment(g_index,dim);
    get_gradient_for_vertices(gi, all_grad);
}
} // end namespace admm

#endif




