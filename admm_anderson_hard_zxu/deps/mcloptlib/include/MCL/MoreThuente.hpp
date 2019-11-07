// The MIT License (MIT)
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
//
// From https://github.com/PatWie/CppNumericalSolvers
//

#ifndef MCL_MORETHUENTE_H
#define MCL_MORETHUENTE_H

#include "Problem.hpp"

namespace mcl {
namespace optlib {

template<typename Scalar, int DIM, typename P>
class MoreThuente {
private:
	typedef Eigen::Matrix<Scalar,DIM,1> VectorX;

public:

	static Scalar linesearch(const VectorX &x, const VectorX &p, P &problem, Scalar alpha_init){
		Scalar alpha = std::abs(alpha_init);
		cvsrch(problem, x, alpha, p);
		return alpha;
	}

	static void cvsrch(P &problem, const VectorX &x0, Scalar &stp, const VectorX &s) {
		int info           = 0;
		int infoc          = 1;
		const Scalar xtol   = 1e-15;
		const Scalar ftol   = 1e-4;
		const Scalar gtol   = 1e-2;
		const Scalar stpmin = 1e-15;
		const Scalar stpmax = 1e15;
		const Scalar xtrapf = 4;
		const int maxfev   = 20;
		int nfev           = 0;
		int dim = x0.rows();

		VectorX g;
		if( DIM == Eigen::Dynamic ){ g = VectorX::Zero(dim); }
		else{ g.setZero(); }

		Scalar f = problem.gradient(x0, g);
		Scalar dginit = g.dot(s);
		if (dginit >= 0.0) {
			// no descent direction
			return;
		}

		bool brackt      = false;
		bool stage1      = true;

		Scalar finit      = f;
		Scalar dgtest     = ftol * dginit;
		Scalar width      = stpmax - stpmin;
		Scalar width1     = 2 * width;
		VectorX x = x0.eval();

		Scalar stx        = 0.0;
		Scalar fx         = finit;
		Scalar dgx        = dginit;
		Scalar sty        = 0.0;
		Scalar fy         = finit;
		Scalar dgy        = dginit;

		Scalar stmin = 0.0;
		Scalar stmax = 0.0;

		const int max_iters = 100000;
		int iter = 0;
		for( ; iter<max_iters; ++iter ){

			// make sure we stay in the interval when setting min/max-step-width
			if (brackt) {
				stmin = std::min(stx, sty);
				stmax = std::max(stx, sty);
			} else {
				stmin = stx;
				stmax = stp + xtrapf * (stp - stx);
			}

			// Force the step to be within the bounds stpmax and stpmin.
			stp = std::max(stp, stpmin);
			stp = std::min(stp, stpmax);

			// Oops, let us return the last reliable values
			if (
			(brackt && ((stp <= stmin) || (stp >= stmax)))
			|| (nfev >= maxfev - 1 ) || (infoc == 0)
			|| (brackt & (stmax - stmin <= xtol * stmax))) {
				stp = stx;
			}

			// test new point
			x = x0 + stp * s;
			f = problem.gradient(x, g);
			nfev++;
			Scalar dg = g.dot(s);
			Scalar ftest1 = finit + stp * dgtest;

			// all possible convergence tests
			if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0))
				info = 6;

			if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest))
				info = 5;
	
			if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest)))
				info = 4;
		
			if (nfev >= maxfev)
				info = 3;
	
			if (brackt & (stmax - stmin <= xtol * stmax))
				info = 2;
	
			if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))
				info = 1;

			// terminate when convergence reached
			if (info != 0)
				return;

			if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol)*dginit))
				stage1 = false;

			if (stage1 & (f <= fx) & (f > ftest1)) {
				Scalar fm = f - stp * dgtest;
				Scalar fxm = fx - stx * dgtest;
				Scalar fym = fy - sty * dgtest;
				Scalar dgm = dg - dgtest;
				Scalar dgxm = dgx - dgtest;
				Scalar dgym = dgy - dgtest;

				cstep( stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

				fx = fxm + stx * dgtest;
				fy = fym + sty * dgtest;
				dgx = dgxm + dgtest;
				dgy = dgym + dgtest;
			} else {
				// this is ugly and some variables should be moved to the class scope
				cstep( stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
			}

			if (brackt) {
				if (fabs(sty - stx) >= 0.66 * width1)
					stp = stx + 0.5 * (sty - stx);

				width1 = width;
				width = fabs(sty - stx);
			}

		} // end while true

		if( iter == max_iters ){
			throw std::runtime_error("MoreThuente::linesearch Error: Reached max_iter");
		}

		return;
	}

	static void cstep(Scalar& stx, Scalar& fx, Scalar& dx, Scalar& sty, Scalar& fy, Scalar& dy, Scalar& stp,
	Scalar& fp, Scalar& dp, bool& brackt, Scalar& stpmin, Scalar& stpmax, int& info) {
		info = 0;
		bool bound = false;

		// Check the input parameters for errors.
		if ((brackt & ((stp <= std::min(stx, sty) ) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
		| (stpmax < stpmin)) {
			return;
		}

		Scalar sgnd = dp * (dx / fabs(dx));

		Scalar stpf = 0;
		Scalar stpc = 0;
		Scalar stpq = 0;

		if (fp > fx) {
			info = 1;
			bound = true;
			Scalar theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
			Scalar s = std::max(theta, std::max(dx, dp));
			Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
			if (stp < stx)
				gamma = -gamma;

			Scalar p = (gamma - dx) + theta;
			Scalar q = ((gamma - dx) + gamma) + dp;
			Scalar r = p / q;
			stpc = stx + r * (stp - stx);
			stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
			if (fabs(stpc - stx) < fabs(stpq - stx))
				stpf = stpc;
			else
				stpf = stpc + (stpq - stpc) / 2;

			brackt = true;
		} else if (sgnd < 0.0) {
			info = 2;
			bound = false;
			Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
			Scalar s = std::max(theta, std::max(dx, dp));
			Scalar gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
			if (stp > stx)
				gamma = -gamma;

			Scalar p = (gamma - dp) + theta;
			Scalar q = ((gamma - dp) + gamma) + dx;
			Scalar r = p / q;
			stpc = stp + r * (stx - stp);
			stpq = stp + (dp / (dp - dx)) * (stx - stp);
			if (fabs(stpc - stp) > fabs(stpq - stp))
				stpf = stpc;
			else
				stpf = stpq;

			brackt = true;
		} else if (fabs(dp) < fabs(dx)) {
			info = 3;
			bound = 1;
			Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
			Scalar s = std::max(theta, std::max( dx, dp));
			Scalar gamma = s * sqrt(std::max(static_cast<Scalar>(0.), (theta / s) * (theta / s) - (dx / s) * (dp / s)));
			if (stp > stx)
				gamma = -gamma;

			Scalar p = (gamma - dp) + theta;
			Scalar q = (gamma + (dx - dp)) + gamma;
			Scalar r = p / q;
			if ((r < 0.0) & (gamma != 0.0)) {
				stpc = stp + r * (stx - stp);
			} else if (stp > stx) {
				stpc = stpmax;
			} else {
				stpc = stpmin;
			}
			stpq = stp + (dp / (dp - dx)) * (stx - stp);
			if (brackt) {
				if (fabs(stp - stpc) < fabs(stp - stpq)) {
					stpf = stpc;
				} else {
					stpf = stpq;
				}
			} else {
				if (fabs(stp - stpc) > fabs(stp - stpq)) {
					stpf = stpc;
				} else {
					stpf = stpq;
				}
			}
		} else {
			info = 4;
			bound = false;
			if (brackt) {
				Scalar theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
				Scalar s = std::max(theta, std::max(dy, dp));
				Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
				if (stp > sty)
					gamma = -gamma;

				Scalar p = (gamma - dp) + theta;
				Scalar q = ((gamma - dp) + gamma) + dy;
				Scalar r = p / q;
				stpc = stp + r * (sty - stp);
				stpf = stpc;
			} else if (stp > stx)
				stpf = stpmax;
			else {
				stpf = stpmin;
			}
		}

		if (fp > fx) {
			sty = stp;
			fy = fp;
			dy = dp;
		} else {
			if (sgnd < 0.0) {
				sty = stx;
				fy = fx;
				dy = dx;
			}
			stx = stp;
			fx = fp;
			dx = dp;
		}

		stpf = std::min(stpmax, stpf);
		stpf = std::max(stpmin, stpf);
		stp = stpf;

		if (brackt & bound) {
			if (sty > stx) {
				stp = std::min(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
			} else {
				stp = std::max(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
			}
		}

		return;

	} // end cstep

};

}
}

#endif
