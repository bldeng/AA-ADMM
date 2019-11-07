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

#include "TetEnergyTerm.hpp"
#include <iostream>
#include "FastSVD.hpp"

using namespace admm;
using namespace Eigen;

//
//	TetEnergyTerm
//

TetEnergyTerm::TetEnergyTerm( const Vec4i &tet_, const std::vector<Vec3> &verts, const Lame &lame_ ) :
    tet(tet_), lame(lame_), volume(0.0), weight(0.0) {

    // Compute inv rest pose
    Matrix<double,3,3> edges; // B
    edges.col(0) = verts[1] - verts[0];
    edges.col(1) = verts[2] - verts[0];
    edges.col(2) = verts[3] - verts[0];
    edges_inv = edges.inverse();

    volume = (edges).determinant() / static_cast<double>(6.0);
    if( volume < 0 ){
        throw std::runtime_error("**TetEnergyTerm Error: Inverted initial tet");
    }

    double k = lame.bulk_modulus(); // a little stiffer than we usually want, but a fine approx.
    weight = std::sqrt(k*volume); // admm weight. If you change this, see TetEnergyTerm::prox.
}

void TetEnergyTerm::get_reduction( std::vector< Eigen::Triplet<double> > &triplets ){

    Matrix<double,4,3> S; // Selector
    S.setZero();
    S(0,0) = -1;	S(0,1) = -1;	S(0,2) = -1;
    S(1,0) =  1;
            S(2,1) =  1;
                    S(3,2) =  1;
    Eigen::Matrix<double,4,3> D = S * edges_inv;
    Eigen::Matrix<double,3,4> Dt = D.transpose(); // Reduction

    const int rows[3] = { 0, 3, 6 };
    const int cols[4] = { 3*tet[0], 3*tet[1], 3*tet[2], 3*tet[3] };
    for( int r=0; r<3; ++r ){
        for( int c=0; c<4; ++c ){
            double value = Dt(r,c);
            for( int j=0; j<3; ++j ){
                triplets.emplace_back( rows[r]+j, cols[c]+j, value );
            }
        }
    }
}

void TetEnergyTerm::prox(const MatX &W, VecX &zi, const VecX &vi ){
    (void)(W);
    (void)(vi);
    typedef Matrix<double,9,1> Vector9d;
    zi_ = zi;
    Matrix<double,3,3> F = Map<Matrix<double,3,3> >(zi.data());
    JacobiSVD< Matrix<double,3,3>, FullPivHouseholderQRPreconditioner > svd(F, ComputeFullU | ComputeFullV);
    Vec3 S = Vec3::Ones();
    // Flip last singular value if inverted
    if( F.determinant() < 1e-16 ){ S[2] = -1.0; }
    // Project onto constraint
    Matrix<double,3,3> proj = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    Vector9d p = Map<Vector9d>(proj.data());
    // Update zi
    // Note, the below only works if w^2 = k*volume.
    zi = 0.5 * ( p + zi );

    // If w^2 != k*volume, use this:
//	double k = lame.bulk_modulus();
//	double kv = k * volume;
//	double w2 = weight*weight;
//	zi = (kv*p + w2*zi) / (w2 + kv);
}

double TetEnergyTerm::prox_for_strain_limiting_energy( VecX &zi ){
    (void)(zi);
    throw std::runtime_error("**Collision Error: Energy not implemented");
}

double TetEnergyTerm::energy( const VecX &vecF ) {
    typedef Matrix<double,9,1> Vector9d;
    Vector9d vecF_ = vecF; // data() function is non-const
    Matrix<double,3,3> F = Map<Matrix<double,3,3> >(vecF_.data());
    JacobiSVD< Matrix<double,3,3>, FullPivHouseholderQRPreconditioner > svd(F, ComputeFullU | ComputeFullV);
    double k = lame.bulk_modulus();
    vecF_ = zi_ - vecF;
    return 0.5 * k * volume * (( svd.singularValues() - Vec3::Ones() ).squaredNorm() + vecF_.squaredNorm());
}

double TetEnergyTerm::energyLBFGS( const VecX &vecF ) {
    typedef Matrix<double,9,1> Vector9d;
    Vector9d vecF_ = vecF; // data() function is non-const
    Matrix<double,3,3> F = Map<Matrix<double,3,3> >(vecF_.data());
    JacobiSVD< Matrix<double,3,3>, FullPivHouseholderQRPreconditioner > svd(F, ComputeFullU | ComputeFullV);
    double k = lame.bulk_modulus();
    return 0.5 * k * volume * ( svd.singularValues() - Vec3::Ones() ).squaredNorm();
}

void TetEnergyTerm::gradient( const VecX &F, VecX &grad ) {
    typedef Matrix<double,9,1> Vector9d;
    Vector9d copyF = F; // data() function is non-const
    Matrix<double,3,3> F_ = Map<Matrix<double,3,3> >(copyF.data());
    JacobiSVD< Matrix<double,3,3>, FullPivHouseholderQRPreconditioner > svd(F_, ComputeFullU | ComputeFullV);

    double k = lame.bulk_modulus();
    Matrix<double,3,3> G = k * volume * (F_ - svd.matrixU() * svd.matrixV().transpose());
    grad = Map<Vector9d>(G.data());
    grad = grad + copyF - zi_;

}

void TetEnergyTerm::get_gradient( const VecX &F, VecX &grad ) {
    //std::cout << "tet get gradient." << std::endl;
    typedef Matrix<double,9,1> Vector9d;
    Vector9d copyF = F; // data() function is non-const
    Matrix<double,3,3> F_ = Map<Matrix<double,3,3> >(copyF.data());
    JacobiSVD< Matrix<double,3,3>, FullPivHouseholderQRPreconditioner > svd(F_, ComputeFullU | ComputeFullV);

    double k = lame.bulk_modulus();
    Matrix<double,3,3> G = k * volume * (F_ - svd.matrixU() * svd.matrixV().transpose());
    grad = Map<Vector9d>(G.data());
}

//
//	HyperElasticTet implementation
//

void HyperElasticTet::prox(const MatX &W, VecX &zi, const VecX &vi ){
    typedef Matrix<double,9,1> Vec9;
    Prox *problem = get_problem();
    problem->set_W(W);
    problem->set_v(vi);
    problem->set_vol(volume);

    Vec9 zcopy = zi;
    solver.minimize( *problem, zcopy );

    zi = zcopy;
}

double HyperElasticTet::prox_for_strain_limiting_energy( VecX &zi ){
    (void)(zi);
    throw std::runtime_error("**Collision Error: Energy not implemented");
}

double HyperElasticTet::energy( const VecX &F_ ){
    Prox *problem = get_problem();
    // Excludes ADMM penalty since x = x0
    return problem->value(F_)*volume;
}

double HyperElasticTet::energyLBFGS( const VecX &F_ ){
    Prox *problem = get_problem();
    return problem->valueLBFGS(F_)*volume;
}

void HyperElasticTet::gradient( const VecX &F_, VecX &grad ){
    Prox *problem = get_problem();
    Vec9 g = grad;
    problem->gradient(F_,g);//*volume

    grad = g;
    //return e;
}

void HyperElasticTet::get_gradient( const VecX &F_, VecX &grad ){
    typedef Matrix<double,3,3> Mat3;
    Prox *problem = get_problem();

    Matrix<double,9,1> Fcopy = F_;
    Matrix<double,3,3> F = Map<Matrix<double,3,3> >(Fcopy.data());

    Mat3 G = F;
    problem->U_gradient(F, G);
    F = volume*G;
    grad = Map<Eigen::Matrix<double,9,1>>(F.data());
}

//
//	NeoHookean
//

double NeoHookeanTet::NHProx::energy_density(const Vec9 &x) const {
    Matrix<double,9,1> xcopy = x;
    Mat3 F = Map<Mat3>(xcopy.data());

    double J = F.determinant();//S[0]*S[1]*S[2];

    F = F.transpose()*F;
    double I_1 = F.trace();
    double I_3 = J*J;

    double log_I3 = std::log( I_3 );
    double t1 = 0.5 * mu * ( I_1 - log_I3 - 3.0 );
    double t2 = 0.125 * lambda * log_I3 * log_I3;
    double r = t1 + t2;
    return r;
}

double NeoHookeanTet::NHProx::value(const Vec9 &x){
    double t1 = energy_density(x); // U(Dx)
    Vec9 quad_vec = vi-x;
    double t2 = 0.5 * k * quad_vec.squaredNorm();// quad penalty
    return t1 + t2;
}

double NeoHookeanTet::NHProx::valueLBFGS(const Vec9 &x){
    double t1 = energy_density(x); // U(Dx)
    return t1;
}

double NeoHookeanTet::NHProx::gradient(const Vec9 &x, Vec9 &grad){

    Matrix<double,9,1> xcopy = x;
    Mat3 F = Map<Mat3>(xcopy.data());
    Mat3 G = Map<Mat3>(grad.data());
    U_gradient(F, G);
    grad = Map<Vec9>(G.data());
    grad = vol * (grad + k * (x-vi));//
    return vol*value(x);
}

void NeoHookeanTet::NHProx::U_gradient(const Mat3 &F, Mat3 &grad){
    Mat3 F_inv_t = F.inverse();
    double J = F.determinant();//S(0, 0) * S(1, 1) * S(2, 2);

    grad = mu * (F - F_inv_t.transpose()) + lambda * std::log(J) * F_inv_t.transpose();
}
//
//	St Venant-Kirchhoff
//

double StVKTet::StVKProx::value(const Vec9 &x){
//    if( x[0]<0.0 || x[1]<0.0 || x[2]<0.0 ){
//		// No Mr. Linesearch, you have gone too far!
//		return std::numeric_limits<float>::max();
//	}
    double t1 = energy_density(x); // U(Dx)
    Vec9 quad_vec = vi-x;
    double t2 = 0.5 * k * quad_vec.squaredNorm(); // quad penalty
    return t1 + t2;
}

double StVKTet::StVKProx::valueLBFGS(const Vec9 &x){
    double t1 = energy_density(x); // U(Dx)
    return t1;
}

double StVKTet::StVKProx::energy_density(const Vec9 &x) const {

    Matrix<double,9,1> xcopy = x;
    Mat3 F = Map<Mat3>(xcopy.data());

    Mat3 I;
    I.setIdentity();

    Mat3 E = 0.5 * (F.transpose()*F-I);
    double E_trace = E.trace();
    Mat3 E_E = E.transpose()*E;
    double r = mu * E_E.trace() + 0.5 * lambda * E_trace * E_trace;
    return r;

}

double StVKTet::StVKProx::gradient(const Vec9 &x, Vec9 &grad){

    Matrix<double,9,1> xcopy = x;
    Mat3 F = Map<Mat3>(xcopy.data());

    Mat3 G = Map<Mat3>(grad.data());
    U_gradient(F, G);
    grad = Map<Vec9>(G.data());
    grad = vol * (grad + k * (x-vi));

    return vol*value(x);
}

void StVKTet::StVKProx::U_gradient(const Mat3 &F, Mat3 &grad){
    Mat3 I;
    I.setIdentity();

    Mat3 E = 0.5 * (F.transpose()*F-I);
    grad = F * (2.0 * mu * E + lambda * E.trace() * I);
}

