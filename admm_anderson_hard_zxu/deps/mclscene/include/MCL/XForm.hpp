// Copyright (c) 2017 University of Minnesota
// 
// MCLSCENE Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
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
//
// By Matt Overby (http://www.mattoverby.net)

#ifndef MCL_XFORM_H
#define MCL_XFORM_H 1

#include "Vec.hpp"

// TODO this class with my own XForm implementation since I don't like Eigen's

namespace mcl {

template <typename T, int dim=3> using XForm = Eigen::Transform<T,dim,Eigen::Affine>;

namespace xform {

	// Makes an identity matrix
	template <typename T> static XForm<T> identity(){
		XForm<T> r;
		r.setIdentity();
		return r;
	}

	// Makes a scale matrix
	// Usage: XForm<float> s = xform::make_scale(1.f, 2.f, 3.f);
	template <typename T> static XForm<T> make_scale(T x, T y, T z){
		XForm<T> r;
		r.setIdentity();
		r.data()[0] = x; r.data()[5] = y; r.data()[10] = z;
		return r;
	}

	// Makes a translation matrix
	// Usage: XForm<float> t = xform::make_trans(1.f, 2.f, 3.f);
	template <typename T> static XForm<T> make_trans(T x, T y, T z){
		XForm<T> r;
		r.setIdentity();
		r.data()[12] = x; r.data()[13] = y; r.data()[14] = z;
		return r;
	}

	template <typename T> static XForm<T> make_trans( const Vec3<T> &t ){
		return make_trans<T>(t[0],t[1],t[2]);
	}

	// Makes a rotation matrix
	// Usage: Xform<float> t = xform::make_rot(45.f, Vec3f(0,1,0));
	template <typename T> static XForm<T> make_rot(T angle_deg, const Vec3<T> &axis){
		T rx = axis[0]; T ry = axis[1]; T rz = axis[2];
		T l = sqrt(rx*rx + ry*ry + rz*rz);
		T angle = angle_deg * M_PI / 180.f;
		T l1 = 1.f/l, x = rx * l1, y = ry * l1, z = rz * l1;
		T s = std::sin(angle), c = std::cos(angle);
		T xs = x*s, ys = y*s, zs = z*s, c1 = 1.f-c;
		T xx = c1*x*x, yy = c1*y*y, zz = c1*z*z;
		T xy = c1*x*y, xz = c1*x*z, yz = c1*y*z;
		T mat[16] = {xx+c,  xy+zs, xz-ys, 0,
				xy-zs, yy+c,  yz+xs, 0,
				xz+ys, yz-xs, zz+c,  0,
				0, 0, 0, 1 };
		XForm<T> r;
		std::memcpy(r.data(), mat, 16*sizeof(T));
		return r;
	}

	// Makes a view matrix
	// Usage: XForm<float> v = xform::make_view(eye, viewdir, Vec3f(0,1,0));
	template <typename T> static XForm<T> make_view(const Vec3<T> &eye, const Vec3<T> &dir, const Vec3<T> &up){
		Vec3<T> w = dir*-1.f; w.normalize();
		Vec3<T> u = up.cross(w);
		Vec3<T> v = w.cross(u);
		XForm<T> r;
		r.setIdentity();
		for(size_t i=0; i<3; ++i){
			r.data()[4*i] = u[i];
			r.data()[4*i+1] = v[i];
			r.data()[4*i+2] = w[i];
		}
		r.data()[12] = -eye.dot(u);
		r.data()[13] = -eye.dot(v);
		r.data()[14] = -eye.dot(w);
		return r;
	}

	// Makes a view matrix (from a lookat point)
	// Usage: XForm<float> v = xform::make_lookat(eye, Vec3f(0,0,0), Vec3f(0,1,0));
	template <typename T> static XForm<T> make_lookat(
		const Vec3<T> &eye, const Vec3<T> &point, const Vec3<T> &up){
		Vec3<T> dir = point-eye;
		return xform::make_view(eye,dir,up);
	}

	// Makes a perspective matrix
	// Usage: XForm<float> p = xform::make_persp(45, width/height, 1e-3f, 1e6f);
	template <typename T> static XForm<T> make_persp(T fov_deg, T aspect, T near, T far){
		T fov = fov_deg * M_PI / 180.f;
		T cossinf = std::cos(fov/2.f) / std::sin(fov/2.f);
		XForm<T> r;
		r.setIdentity();
		r.data()[0] = cossinf/aspect;
		r.data()[5] = cossinf;
		r.data()[10] = -(near+far)/(far-near);
		r.data()[14] = -(2.f*near*far)/(far-near);
		r.data()[11] = -1.f;
		r.data()[15] = 0.f;
		return r;
	}

	template <typename T> static inline std::string to_string( const XForm<T> &xf ){
		std::stringstream ss;
		ss << xf.data()[0];
		for( int i=1; i<16; ++i ){ ss << ' ' << xf.data()[i]; }
		return ss.str();
	}

	template <typename T> static inline XForm<T> from_string( const std::string &s ){
		std::stringstream ss; ss << s;
		XForm<T> result;
		// TODO some testing to make sure there are tokens left
		for( int i=0; i<16; ++i ){ ss >> result.data()[i]; }
		return result;
	}

} // ns xform

} // ns mcl

#endif
