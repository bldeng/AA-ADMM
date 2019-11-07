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

#ifndef MCL_NONLINEARCG_H
#define MCL_NONLINEARCG_H

#include "Armijo.hpp"
#include "Minimizer.hpp"

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM>
class NonLinearCG : public Minimizer<Scalar,DIM> {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VectorX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatrixX;

public:
	int max_iters;

	NonLinearCG() : max_iters(100) {}
	void set_max_iters( int iters ){ max_iters = iters; }
	void set_verbose( int v ){ (void)(v); } // TODO

	int minimize(Problem<Scalar,DIM> &problem, VectorX &x){

		VectorX grad, grad_old, p, x_last;
		if( DIM == Eigen::Dynamic ){
			int dim = x.rows();
			x_last.resize(dim);
			grad.resize(dim);
			grad_old.resize(dim);
			p.resize(dim);
		}

		int iter=0;
		for( ; iter<max_iters; ++iter ){

			problem.gradient(x, grad);

			if( iter==0 ){ p = -grad; }
			else {
				Scalar beta = grad.dot(grad) / (grad_old.dot(grad_old));
				p = -grad + beta*p;
			}

			Scalar alpha = Armijo<Scalar, DIM, decltype(problem)>::linesearch(x, p, problem, 1);

			if( alpha <= 0 ){
				printf("NonLinearCG::minimize: Failure in linesearch\n");
				return Minimizer<Scalar,DIM>::FAILURE;
			}

			x_last = x;
			x += alpha*p;
			grad_old = grad;

			if( problem.converged(x_last,x,grad) ){ break; }
		}
		return iter;
	} // end minimize

};

} // ns optlib
} // ns mcl

#endif
