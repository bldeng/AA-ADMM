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

//
// Several static functions for projecting a point onto a geometric surface.
// Will return a point on the surface that is nearest to the one given./
//


#ifndef MCL_PROJECTION_H
#define MCL_PROJECTION_H 1

#include <math.h>
#include "Vec.hpp"

namespace mcl {
namespace projection {

	//	Projection on Triangle
	template <typename T> static Vec3<T> point_on_triangle( const Vec3<T> &point, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3 );

	//	Projection on Sphere
	template <typename T> static Vec3<T> point_on_sphere( const Vec3<T> &point, const Vec3<T> &center, const T &rad );

	//	Projection on a Box
	template <typename T> static Vec3<T> point_on_box( const Vec3<T> &point, const Vec3<T> &bmin, const Vec3<T> &bmax );

	//	Point in tet
	template <typename T> static bool point_in_tet( const Vec3<T> &point, const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3 );

	//	Helper functions
	template <typename T> static T myclamp( const T &val ){ return val < 0 ? 0 : (val > 1 ? 1 : val); }

}; // end namespace Projection

//
//	Implementation
//

template <typename T>
Vec3<T> projection::point_on_triangle( const Vec3<T> &point, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3 ){

	Vec3<T> edge0 = p2 - p1;
	Vec3<T> edge1 = p3 - p1;
	Vec3<T> v0 = p1 - point;

	T a = edge0.dot( edge0 );
	T b = edge0.dot( edge1 );
	T c = edge1.dot( edge1 );
	T d = edge0.dot( v0 );
	T e = edge1.dot( v0 );
	T det = a*c - b*b;
	T s = b*e - c*d;
	T t = b*d - a*e;

	const T zero(0);
	const T one(1);

	if ( s + t < det ) {
		if ( s < zero ) {
		    if ( t < zero ) {
			if ( d < zero ) {
			    s = myclamp( -d/a );
			    t = zero;
			}
			else {
			    s = zero;
			    t = myclamp( -e/c );
			}
		    }
		    else {
			s = zero;
			t = myclamp( -e/c );
		    }
		}
		else if ( t < zero ) {
		    s = myclamp( -d/a );
		    t = zero;
		}
		else {
		    T invDet = one / det;
		    s *= invDet;
		    t *= invDet;
		}
	}
	else {
		if ( s < zero ) {
		    T tmp0 = b+d;
		    T tmp1 = c+e;
		    if ( tmp1 > tmp0 ) {
			T numer = tmp1 - tmp0;
			T denom = a-T(2)*b+c;
			s = myclamp( numer/denom );
			t = one-s;
		    }
		    else {
			t = myclamp( -e/c );
			s = zero;
		    }
		}
		else if ( t < zero ) {
		    if ( a+d > b+e ) {
			T numer = c+e-b-d;
			T denom = a-T(2)*b+c;
			s = myclamp( numer/denom );
			t = one-s;
		    }
		    else {
			s = myclamp( -e/c );
			t = zero;
		    }
		}
		else {
		    T numer = c+e-b-d;
		    T denom = a-T(2)*b+c;
		    s = myclamp( numer/denom );
		    t = one - s;
		}
	}

	return ( p1 + edge0*s + edge1*t );

} // end project triangle


template <typename T>
Vec3<T> projection::point_on_sphere( const Vec3<T> &point, const Vec3<T> &center, const T &rad ){
	Vec3<T> dir = point-center;
	dir.normalize();
	return ( center + dir*rad );
} // end project sphere


template <typename T>
Vec3<T> projection::point_on_box( const Vec3<T> &point, const Vec3<T> &bmin, const Vec3<T> &bmax ){
	// Loops through axes and moves point to nearest surface
	Vec3<T> x = point;
	T dx = std::numeric_limits<T>::max();
	for( int i=0; i<3; ++i ){
		T dx_max = std::abs(bmax[i]-point[i]);
		T dx_min = std::abs(bmin[i]-point[i]);
		if( dx_max < dx ){
			x = point;
			x[i] = bmax[i];
			dx = dx_max;
		}
		if( dx_min < dx ){
			x = point;
			x[i] = bmin[i];
			dx = dx_min;
		}
	}
	return x;
} // end project box


template <typename T>
bool check_norm( const Vec3<T> &point,
	const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3 ){
	const Vec3<T> n = (p1 - p0).cross(p2 - p0);
	const T dp3 = n.dot(p3 - p0);
	const T dp = n.dot(point - p0);
	return (dp3*dp>0);
}

template <typename T>
bool projection::point_in_tet( const Vec3<T> &point,
	const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3 ){
	return check_norm<T>(point, p0, p1, p2, p3) && check_norm<T>(point, p1, p2, p3, p0) &&
		check_norm<T>(point, p2, p3, p0, p1) && check_norm<T>(point, p3, p0, p1, p2);
}

} // end namespace mcl

#endif

