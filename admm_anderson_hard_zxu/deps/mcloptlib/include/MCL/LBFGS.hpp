// The MIT License (MIT)
// Copyright (c) 2017 University of Minnesota
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

#ifndef MCL_LBFGS_H
#define MCL_LBFGS_H

#include "Armijo.hpp"
#include "MoreThuente.hpp"
#include "Minimizer.hpp"

namespace mcl {
namespace optlib {

enum LINE_SEARCH_ALGORITHM
{
    ///
    /// Backtracking method with the Armijo condition.
    /// The backtracking method finds the step length such that it satisfies
    /// the sufficient decrease (Armijo) condition,
    /// \f$f(x + a \cdot d) \le f(x) + \beta' \cdot a \cdot g(x)^T d\f$,
    /// where \f$x\f$ is the current point, \f$d\f$ is the current search direction,
    /// \f$a\f$ is the step length, and \f$\beta'\f$ is the value specified by
    /// \ref LBFGSParam::ftol. \f$f\f$ and \f$g\f$ are the function
    /// and gradient values respectively.
    ///
    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1,

    ///
    /// The backtracking method with the defualt (regular Wolfe) condition.
    /// An alias of `LBFGS_LINESEARCH_BACKTRACKING_WOLFE`.
    ///
    LBFGS_LINESEARCH_BACKTRACKING = 2,

    ///
    /// Backtracking method with regular Wolfe condition.
    /// The backtracking method finds the step length such that it satisfies
    /// both the Armijo condition (`LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`)
    /// and the curvature condition,
    /// \f$g(x + a \cdot d)^T d \ge \beta \cdot g(x)^T d\f$, where \f$\beta\f$
    /// is the value specified by \ref LBFGSParam::wolfe.
    ///
    LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2,

    ///
    /// Backtracking method with strong Wolfe condition.
    /// The backtracking method finds the step length such that it satisfies
    /// both the Armijo condition (`LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`)
    /// and the following condition,
    /// \f$\vert g(x + a \cdot d)^T d\vert \le \beta \cdot \vert g(x)^T d\vert\f$,
    /// where \f$\beta\f$ is the value specified by \ref LBFGSParam::wolfe.
    ///
    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3
};

// L-BFGS implementation based on Nocedal & Wright Numerical Optimization book (Section 7.2)
// DIM = dimension of the problem
// M = history window
//
// Original Author: Ioannis Karamouzas
//
template<typename Scalar, int DIM, int M=8>
class LBFGS : public Minimizer<Scalar,DIM> {
private:
    typedef Eigen::Matrix<Scalar,DIM,1> VecX;
    typedef Eigen::Matrix<Scalar,DIM,M> MatM;
    typedef Eigen::Matrix<Scalar,DIM,1> VecM;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Map<Vector> MapVec;

public:
	int max_iters;
	bool show_denom_warning; // Print out warning for zero denominators

    LBFGS() : max_iters(100), show_denom_warning(false) {}
	void set_max_iters( int iters ){ max_iters = iters; }
	void set_verbose( int v ){ show_denom_warning = v > 0 ? true : false; }

private:

    Matrix                    m_s;      // History of the s vectors
    Matrix                    m_y;      // History of the y vectors
    Vector                    m_ys;     // History of the s'y values
    Vector                    m_alpha;  // History of the step lengths
    Vector                    m_fx;     // History of the objective function values
    Vector                    m_xp;     // Old x
    VecX                    m_grad;   // New gradient
    VecX                    m_gradp;  // Old gradient
    Vector                    m_drt;    // Moving direction

    int m              = 6;
    Scalar epsilon        = Scalar(1e-6);
    int past           = 1;
    Scalar delta          = Scalar(1e-16);
    int linesearch     = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    int max_linesearch = 2000;
    Scalar min_step       = Scalar(1e-20);
    Scalar max_step       = Scalar(1e+20);
    Scalar ftol           = Scalar(1e-4);
    Scalar wolfe          = Scalar(0.9);

    inline void reset(int n)
    {
        m_s.resize(n, m);
        m_y.resize(n, m);
        m_ys.resize(m);
        m_alpha.resize(m);
        m_xp.resize(n);
        //m_grad.resize(n);
        //m_gradp.resize(n);
        m_drt.resize(n);
        if(past > 0) {
            m_fx.resize(past);
        }
    }

public:
    //LineSearch = LineSearchBacktracking
    void LineSearch(Problem<Scalar,DIM> &problem, Scalar& fx, VecX& x, VecX& grad,
                               Scalar& step, const VecX& drt, const VecX& xp)
        {
            // Decreasing and increasing factors
            const Scalar dec = 0.5;
            const Scalar inc = 2.1;

            // Check the value of step
            if(step <= Scalar(0))
                std::invalid_argument("'step' must be positive");

            // Save the function value at the current x
            const Scalar fx_init = fx;
            // Projection of gradient on the search direction
            const Scalar dg_init = grad.dot(drt);
            // Make sure d points to a descent direction
            if(dg_init > 0)
                std::logic_error("the moving direction increases the objective function value");

            const Scalar dg_test = ftol * dg_init;
            Scalar width;

            int iter;
            for(iter = 0; iter < max_linesearch; iter++)
            {
                x.noalias() = xp + step * drt;
                // Evaluate this candidate
                fx = problem.gradient(x, grad);

                if(fx > fx_init + step * dg_test)
                {
                    width = dec;
                } else {
                    // Armijo condition is met
                    if(linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO) {
                        break;}

                    const Scalar dg = grad.dot(drt);
                    if(dg < wolfe * dg_init)
                    {
                        width = inc;
                    } else {
                        // Regular Wolfe condition is met
                        if(linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE)
                            break;

                        if(dg > -wolfe * dg_init)
                        {
                            width = dec;
                        } else {
                            // Strong Wolfe condition is met
                            break;
                        }
                    }
                }

                if(iter >= max_linesearch)
                    throw std::runtime_error("the line search routine reached the maximum number of iterations");

                if(step < min_step)
                    throw std::runtime_error("the line search step became smaller than the minimum value allowed");

                if(step > max_step)
                    throw std::runtime_error("the line search step became larger than the maximum value allowed");

                step *= width;
            }
            //std::cout << "line search iter = " << iter << std::endl;
        }

    int minimize(Problem<Scalar,DIM> &problem, VecX &x){
        const int n = x.size();
        const int fpast = past;
        reset(n);

        // Evaluate function and compute gradient
        double fx = problem.gradient(x, m_grad);
        Scalar xnorm = x.norm();
        Scalar gnorm = m_grad.norm();

        if(fpast > 0)
            m_fx(0) = fx;

        // Early exit if the initial x is already a minimizer
        if(gnorm <= epsilon * std::max(xnorm, Scalar(1.0)))
        {
            return 1;
        }

        // Initial direction
        m_drt.noalias() = -m_grad;
        // Initial step
        Scalar step = Scalar(1.0) / m_drt.norm();

        int k = 1;
        int end = 0;
        for( ; ; )
        {
            // Save the curent x and gradient
            m_xp.noalias() = x;
            m_gradp.noalias() = m_grad;

            // Line search to update x, fx and gradient
            LineSearch(problem, fx, x, m_grad, step, m_drt, m_xp);

            // New x norm and gradient norm
            xnorm = x.norm();
            gnorm = m_grad.norm();

            // Convergence test -- gradient
            if(gnorm <= epsilon * std::max(xnorm, Scalar(1.0)))
            {
                return k;
            }
            // Convergence test -- objective function value
            if(fpast > 0)
            {
                //if(k >= fpast && std::abs((m_fx[k % fpast] - fx) / fx) < delta)
                if(k >= fpast && std::abs(m_fx(k % fpast) - fx) < delta)
                    return k;

                m_fx(k % fpast) = fx;
            }
            // Maximum number of iterations
            if(max_iters != 0 && k >= max_iters)
            {
                return k;
            }

            // Update s and y
            // s_{k+1} = x_{k+1} - x_k
            // y_{k+1} = g_{k+1} - g_k
            MapVec svec(&m_s(0, end), n);
            MapVec yvec(&m_y(0, end), n);
            svec.noalias() = x - m_xp;
            yvec.noalias() = m_grad - m_gradp;

            // ys = y's = 1/rho
            // yy = y'y
            Scalar ys = yvec.dot(svec);
            Scalar yy = yvec.squaredNorm();
            m_ys(end) = ys;

            // Recursive formula to compute d = -H * g
            m_drt.noalias() = -m_grad;
            int bound = std::min(m, k);
            end = (end + 1) % m;
            int j = end;
            for(int i = 0; i < bound; i++)
            {
                j = (j + m - 1) % m;
                MapVec sj(&m_s(0, j), n);
                MapVec yj(&m_y(0, j), n);
                m_alpha(j) = sj.dot(m_drt) / m_ys(j);
                m_drt.noalias() -= m_alpha(j) * yj;
            }

            m_drt *= (ys / yy);

            for(int i = 0; i < bound; i++)
            {
                MapVec sj(&m_s(0, j), n);
                MapVec yj(&m_y(0, j), n);
                Scalar beta = yj.dot(m_drt) / m_ys(j);
                m_drt.noalias() += (m_alpha(j) - beta) * sj;
                j = (j + 1) % m;
            }

            // step = 1.0 as initial guess
            step = Scalar(1.0);
            k++;
    }
    }
};

}
}

#endif
