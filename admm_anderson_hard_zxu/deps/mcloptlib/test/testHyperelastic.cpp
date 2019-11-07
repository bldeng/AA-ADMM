// The MIT License (MIT)
// Copyright (c) 2017 Matt Overby
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <iostream>
#include "MCL/LBFGS.hpp"
#include "MCL/NonLinearCG.hpp"
#include "MCL/Newton.hpp"
#include <memory>

using namespace mcl::optlib;
typedef std::shared_ptr< Minimizer<double,3> > StVKSolver;

//
//	Test on a near-zero volume isotropic NeoHookean tet
//	which usually causes linesearch rate-blocking.
//	Note that this isn't normal NH, it's NH subject to
//	a quadratic that penalizes change from its initial state (x0).
//
class NeoHookean : public Problem<double,3> {
public:
	typedef Eigen::Matrix<double,3,1> Vec3;
	typedef Eigen::Matrix<double,3,3> Mat3;

	double mu, lambda, k;
	Vec3 x0;
	NeoHookean() : mu(3.33556e+06), lambda(1.66444e+09), k(1.66667e+09),
		x0(2.74498,1.18735,0.863) {}

	double energy_density(const Vec3 &x) const {
		double J = x[0]*x[1]*x[2];
		double I_1 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
		double I_3 = J*J;
		double log_I3 = std::log( I_3 );
		double t1 = 0.5 * mu * ( I_1 - log_I3 - 3.0 );
		double t2 = 0.125 * lambda * log_I3 * log_I3;
		double r = t1 + t2;
		return r;
	}

	double value(const Vec3 &x){
		if( x[0]<0.0 || x[1]<0.0 || x[2]<0.0 ){
			// No Mr. Linesearch, you have gone too far!
			return std::numeric_limits<float>::max();
		}
		double t1 = energy_density(x); // U(Dx)
		double t2 = (k*0.5) * (x-x0).squaredNorm(); // quad penalty
		return t1 + t2;
	}

	double gradient(const Vec3 &x, Vec3 &grad){
		double J = x[0]*x[1]*x[2];
		if( J <= 0.0 ){
			throw std::runtime_error("BAD DET IN GRAD");
		} else {
			Eigen::Vector3d invSigma(1.0/x[0],1.0/x[1],1.0/x[2]);
			grad = (mu * (x - invSigma) + lambda * std::log(J) * invSigma) + k*(x-x0);
		}
		return value(x);
	}

	void hessian(const Vec3 &x, Mat3 &hess){
		double det = x[0]*x[1]*x[2]; // J
		if( det <= 0.0 ){
			throw std::runtime_error("BAD DET IN HESSIAN");
		}

		Vec3 inv_x(1.0/x[0],1.0/x[1],1.0/x[2]);
		Mat3 invXSqMat = Mat3::Zero();
		invXSqMat(0,0) = 1.0 / (x[0]*x[0]);
		invXSqMat(1,1) = 1.0 / (x[1]*x[1]);
		invXSqMat(2,2) = 1.0 / (x[2]*x[2]);
		double logJ = std::log(det);
		Mat3 I = Mat3::Identity();
		hess = mu*(I - 2.0*invXSqMat) +
			lambda * logJ * invXSqMat +
			lambda * Mat3(inv_x * inv_x.transpose()) +
			k*I; // quad penalty
	}

	bool converged(const Vec3 &x0, const Vec3 &x1, const Vec3 &grad){
		return ( grad.norm() < 1e-10 || (x1-x0).norm() < 1e-10 );
	}
};


class StVK : public Problem<double,3> {
public:
	typedef Eigen::Matrix<double,3,1> Vec3;

	static inline double ddot( Vec3 &a, Vec3 &b ) { return (a*b.transpose()).trace(); }
	static inline double v3trace( const Vec3 &v ) { return v[0]+v[1]+v[2]; }

	double mu, lambda, k;
	Vec3 x0;
	StVK() : mu(3.33556e+06), lambda(1.66444e+09), k(1.66667e+09),
		x0(3.70093,0.0247722,0.0210168) {}

	double value(const Vec3 &x){
		if( x[0]<0.0 || x[1]<0.0 || x[2]<0.0 ){
			return std::numeric_limits<float>::max(); // No Mr. Linesearch, you have gone too far!
		}
		double t1 = energy_density(x); // U(Dx)
		double t2 = (k*0.5) * (x-x0).squaredNorm(); // quad penalty
		return t1 + t2;
	}

	double energy_density(const Vec3 &x) const {
		Vec3 x2( x[0]*x[0], x[1]*x[1], x[2]*x[2] );
		Vec3 st = 0.5 * ( x2 - Vec3(1,1,1) ); // strain tensor
		double st_tr2 = v3trace(st)*v3trace(st);
		double r = ( mu * ddot( st, st ) + ( lambda * 0.5 * st_tr2 ) );
		return r;
	}


	double gradient(const Vec3 &x, Vec3 &grad){
		Vec3 term1(
			mu * x[0]*(x[0]*x[0] - 1.0),
			mu * x[1]*(x[1]*x[1] - 1.0),
			mu * x[2]*(x[2]*x[2] - 1.0)
		);
		Vec3 term2 = 0.5 * lambda * ( x.dot(x) - 3.0 ) * x;
		grad = term1 + term2 + k*(x-x0);
		return value(x);
	}

	bool converged(const Vec3 &x0, const Vec3 &x1, const Vec3 &grad){
		return ( grad.norm() < 1e-10 || (x1-x0).norm() < 1e-10 );
	}
};


int main(){
	srand(100);

	std::vector< StVKSolver > solver_list;
//	solver_list.emplace_back( std::make_shared< NonLinearCG<double,3> >( NonLinearCG<double,3>() ) );
	solver_list.emplace_back( std::make_shared< LBFGS<double,3> >( LBFGS<double,3>() ) );
//	solver_list.emplace_back( std::make_shared< Newton<double,3> >( Newton<double,3>() ) );

	NeoHookean nh;
	StVK stvk;
	bool use_nh = true; // set to false to use stvk
	bool success = true;

	//
	//	Loop through solvers and test
	//
	int n_solvers = solver_list.size();
	for( int i=0; i<n_solvers; ++i ){

		solver_list[i]->set_max_iters(100);
		solver_list[i]->set_verbose(1);

		Eigen::Vector3d x; // init guess and final solution
		int iters = -1; // <= 0 is an error
		if( use_nh ){
			x = nh.x0;
			iters = solver_list[i]->minimize( nh, x );
		} else {
			x = stvk.x0;
			iters = solver_list[i]->minimize( stvk, x );
		}

		for( int i=0; i<3; ++i ){
			if( std::isnan(x[i]) || std::isinf(x[i]) ){
				std::cerr << "(" << i << ") Bad values in x: " << x[i] << std::endl;
				success = false;
			}
		}

		Eigen::Vector3d grad;
		if( use_nh ){
			nh.gradient(x,grad);
		} else {
			stvk.gradient(x,grad);
		}

		double gn = grad.norm();
		std::cout << "iters: " << iters << ", grad norm: " << gn << std::endl;
		if( gn > 1e-4 ){
			std::cerr << "(" << i << ") Failed to minimize: |grad| = " << gn << std::endl;
			success = false;
		}

		if( iters == Minimizer<double,3>::FAILURE ){
			std::cerr << "(" << i << ") Failed to minimize" << std::endl;
			success = false;
		}

	} // end loop solvers

	if( success ){
		std::cout << "Success" << std::endl;
		return EXIT_SUCCESS;
	}
	return EXIT_FAILURE;
}



