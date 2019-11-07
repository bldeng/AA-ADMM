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

#ifndef MCL_MINIMIZER_H
#define MCL_MINIMIZER_H

#include "Problem.hpp"

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM>
class Minimizer {
public:
    typedef Eigen::Matrix<Scalar,DIM,1> VectorX;
	static const int FAILURE = -1; // returned by minimize if an error is encountered

	virtual void set_max_iters( int iters ) = 0;
	virtual void set_verbose( int v ) = 0;
	virtual int minimize(Problem<Scalar,DIM> &problem, VectorX &x) = 0;
};

} // ns optlib
} // ns mcl

#endif
