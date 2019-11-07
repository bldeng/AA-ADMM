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

#ifndef MCL_ARMIJO_H
#define MCL_ARMIJO_H

#include "Problem.hpp"

namespace mcl {
namespace optlib {

// Backtracking-Armijo
template<typename Scalar, int DIM, typename P>
class Armijo {
public:
	typedef Eigen::Matrix<Scalar,DIM,1> VectorX;
	typedef Eigen::Matrix<Scalar,DIM,DIM> MatrixX;

	static Scalar linesearch(const VectorX &x, const VectorX &p, P &problem, Scalar alpha_init) {

		const Scalar tau = 0.7;
		const Scalar beta = 0.3;
		const int max_iter = 1000000;
		Scalar alpha = std::abs(alpha_init);
		VectorX grad;
		if( DIM == Eigen::Dynamic ){ grad = VectorX::Zero(x.rows()); }
		Scalar f_x = problem.gradient(x, grad);
		Scalar bgdp = beta*grad.dot(p);

		int iter = 0;
		for( iter=0; iter < max_iter; ++iter ){
			Scalar f_xap = problem.value(x + alpha*p);
			Scalar f_x_a = f_x + alpha*bgdp; // Armijo condition
			if( f_xap <= f_x_a ){ break; }
			alpha *= tau;
		}

		if( iter == max_iter ){
			printf("Armijo::linesearch Error: Reached max_iters\n");
			return -1;
		}

		return alpha;
	}
};

}
}

#endif
