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

#ifndef MCL_VEC3_H
#define MCL_VEC3_H 1

// Helper to hush api warnings
#define MCL_UNUSED(x) (void)(x)

#include <Eigen/Geometry>

namespace mcl {

	// Common types that I don't feel like typing out all the time:
	template <typename T> using Vec4 = Eigen::Matrix<T,4,1>;
	template <typename T> using Vec3 = Eigen::Matrix<T,3,1>;
	template <typename T> using Vec2 = Eigen::Matrix<T,2,1>;
	typedef Vec4<float> Vec4f;
	typedef Vec3<float> Vec3f;
	typedef Vec2<float> Vec2f;
	typedef Vec4<double> Vec4d;
	typedef Vec3<double> Vec3d;
	typedef Vec2<double> Vec2d;
	typedef Vec4<int> Vec4i;
	typedef Vec3<int> Vec3i;
	typedef Vec2<int> Vec2i;

namespace vec {

	template <typename T> // Returns a normalized vector
	static inline Vec3<T> normalized(const Vec3<T> &v){ Vec3<T> t=v; t.normalize(); return t; }

	template <typename T, size_t D> // Vec output as transpose with spaces instead of tabs
	static inline std::string to_str(const Eigen::Matrix<T,D,1> &v){
		std::stringstream ss; ss << v[0];
		for( size_t i=1; i<D; ++i ){ ss << ' ' << v[i]; }
		return ss.str();
	}

	template <typename T> // Compute barycentric coords for a point on a triangle
	static inline Vec3<T> barycoords(const Vec3<T> &p, const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2){
		Vec3<T> v0 = p1 - p0, v1 = p2 - p0, v2 = p - p0;
		T d00 = v0.dot(v0);
		T d01 = v0.dot(v1);
		T d11 = v1.dot(v1);
		T d20 = v2.dot(v0);
		T d21 = v2.dot(v1);
		T invDenom = 1.0 / (d00 * d11 - d01 * d01);
		Vec3<T> r;
		r[1] = (d11 * d20 - d01 * d21) * invDenom;
		r[2] = (d00 * d21 - d01 * d20) * invDenom;
		r[0] = 1.0 - r[1] - r[2];
		return r;
	}

	template <typename T> // scalar triple product
	static inline T scalar_triple_product( const Vec3<T> &u, const Vec3<T> &v, const Vec3<T> &w ){ return u.dot(v.cross(w)); }
	
	template <typename T> // Compute barycentric coords for a point in a tet
	static inline Vec4<T> barycoords(const Vec3<T> &p, const Vec3<T> &a, const Vec3<T> &b, const Vec3<T> &c, const Vec3<T> &d){
		Vec3<T> vap = p - a;
		Vec3<T> vbp = p - b;
		Vec3<T> vab = b - a;
		Vec3<T> vac = c - a;
		Vec3<T> vad = d - a;
		Vec3<T> vbc = c - b;
		Vec3<T> vbd = d - b;
		T va6 = scalar_triple_product(vbp, vbd, vbc);
		T vb6 = scalar_triple_product(vap, vac, vad);
		T vc6 = scalar_triple_product(vap, vad, vab);
		T vd6 = scalar_triple_product(vap, vab, vac);
		T v6 = 1.0 / scalar_triple_product(vab, vac, vad);
		return Vec4<T>(va6*v6, vb6*v6, vc6*v6, vd6*v6);
	}
	
	template <typename T> // Spherical coords to cartesian
	static inline Vec3<T> spherical_to_cartesian(T theta, T phi){
		T sin_t = std::sin(theta); T cos_t = std::cos(theta);
		T sin_p = std::sin(phi); T cos_p = std::cos(phi);
		return Vec3<T>( sin_t * sin_p, sin_t * cos_p, cos_t );
	}

	template <typename T>  // Cartesian coords to spherical
	static inline Vec2<T> cartesian_to_spherical(const Vec3<T> &v){
		Vec2<T> r( std::acos(v[2]), std::atan2(v[1], v[0]) );
		if(r[1] < 0){ r[1] += 2*M_PI; }
		return r;
	}

} // end namespace vec

// Randoms (u1, u2, etc...): 0 to 1
// Putting it here until I find a better spot
namespace sample {

	template<typename T> // Uniformly samples a cone (e.g. spotlight)
	static inline Vec3<T> uniform_cone( T u1, T u2, T max_theta ){
		T cos_theta = (1 - u1) + u1 * std::cos(max_theta);
		T sin_theta = std::sqrt(1 - cos_theta*cos_theta);
		T phi = u2 * 2 * M_PI;
		return Vec3<T>( std::cos(phi)*sin_theta, std::sin(phi)*sin_theta, cos_theta );
	}

	template<typename T> // Cosine weighted hemisphere sampling (e.g. diffuse reflection)
	static inline Vec3<T> cosine_hemisphere( T u1, T u2 ){
		T r = std::sqrt( u1 );
		T theta = 2 * M_PI * u2;
		return Vec3<T>( r * std::cos(theta), r * std::sin(theta), std::sqrt( std::max(T(0), T(1.f-u1)) ) );
	}

}; // end namespace sample

} // end namespace mcl

/*
//
//	trimesh and mcl::Vec xforms:
//
namespace trimesh {

	template <typename T, typename U>
	static inline mcl::Vec3<T> operator*(const trimesh::XForm<U> &m, const mcl::Vec3<T> &v){
		mcl::Vec3<T> r;
		r[0] = m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12];
		r[1] = m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13];
		r[2] = m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14];
		return r;
	}

	template <typename T, typename U>
	static inline mcl::Vec4<T> operator*(const trimesh::XForm<U> &m, const mcl::Vec4<T> &v){
		mcl::Vec4<T> r;
		r[0] = m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12]*v[3];
		r[1] = m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13]*v[3];
		r[2] = m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14]*v[3];
		r[3] = m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]*v[3];
		return r;
	}

} // end namespace trimesh
*/

#endif
