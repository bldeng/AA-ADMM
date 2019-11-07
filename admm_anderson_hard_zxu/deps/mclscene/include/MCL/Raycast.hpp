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
// Several common classes and static functions related to ray tracing
//

#ifndef MCL_RAYINTERSECT_H
#define MCL_RAYINTERSECT_H 1

#include <memory>
#include "Vec.hpp"

namespace mcl {
namespace raycast {

	// Ray
	template<typename T> class Ray {
	public:
		Ray() : origin(0,0,0), direction(0,0,-1), eps(1e-5) {}
		Ray( Vec3<T> o, Vec3<T> d, T e=T(1e-5) ) : origin(o), direction(d), eps(e) {}
		Vec3<T> origin, direction;
		T eps;
	};

	// Payload
	template<typename T> class Payload {
	public:
		Payload() : t_min(1e-5), t_max(std::numeric_limits<T>::max()),
			bary(0,0,0), n(0,0,0), hit_point(0,0,0), material(-1) {}
		mutable T t_min, t_max;
		mutable Vec3<T> bary, n, hit_point;
		mutable int material; // usually an index into an array
	};

	//
	//	Intersection functions and utility
	//

	// Ideal specular reflection
	template<typename T> static Vec3<T> reflect( const Vec3<T> &incident, const Vec3<T> &norm ){
		return ( incident - T(2) * norm * norm.dot( incident ) );
	}

	// ray -> axis aligned bounding box
	// Returns true/false only and does not set the payload.
	template<typename T> static bool ray_aabb( const Ray<T> *ray,
		const Vec3<T> &min, const Vec3<T> &max );

	// ray -> triangle without smoothed normals
	template<typename T> static bool ray_triangle( const Ray<T> *ray,
		const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2,
		Payload<T> *payload );

	// ray -> triangle with smoothed normals
	template<typename T> static bool ray_triangle( const Ray<T> *ray,
		const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2,
		const Vec3<T> &n0, const Vec3<T> &n1, const Vec3<T> &n2,
		Payload<T> *payload );

	// ray -> sphere
	template<typename T> static bool ray_sphere( const Ray<T> *ray,
		const Vec3<T> &center, const T &radius,
		Payload<T> *payload );

} // end namespace raycast

//
//	Implementation below
//

// ray -> axis aligned bounding box
template<typename T> static bool mcl::raycast::ray_aabb( const Ray<T> *ray,
	const Vec3<T> &min, const Vec3<T> &max ){

	// First check if origin is inside AABB
	bool in_box = true;
	for( int i=0; i<3; ++i ){
		if( ray->origin[i] > max[i] ){ in_box = false; break; }
		if( ray->origin[i] < min[i] ){ in_box = false; break; }
	}
	if( in_box ){ return true; }

	// Compute edge t values
	T txmin=0.f, txmax=0.f;
	T dirX = 1.f / ray->direction[0];
	if( dirX >= 0.0 ){
		txmin = dirX * ( min[0] - ray->origin[0] );
		txmax = dirX * ( max[0] - ray->origin[0] );
	}
	else{
		txmax = dirX * ( min[0] - ray->origin[0] );
		txmin = dirX * ( max[0] - ray->origin[0] );
	}

	T tymin=0.f, tymax=0.f;
	T dirY = 1.f / ray->direction[1];
	if( ray->direction[1] >= 0.0 ){
		tymin = dirY * ( min[1] - ray->origin[1] );
		tymax = dirY * ( max[1] - ray->origin[1] );
	}
	else{
		tymax = dirY * ( min[1] - ray->origin[1] );
		tymin = dirY * ( max[1] - ray->origin[1] );
	}

	// Now check x/y axis
	if( txmin > tymax || tymin > txmax ){ return false; }

	T tzmin=0.f, tzmax=0.f;
	T dirZ = 1.f / ray->direction[2];
	if( ray->direction[2] >= 0.0 ){
		tzmin = dirZ * ( min[2] - ray->origin[2] );
		tzmax = dirZ * ( max[2] - ray->origin[2] );
	}
	else{
		tzmax = dirZ * ( min[2] - ray->origin[2] );
		tzmin = dirZ * ( max[2] - ray->origin[2] );
	}

	// Finally, check z axis
	if( txmin > tzmax || tzmin > txmax ){ return false; }
	if( tymin > tzmax || tzmin > tymax ){ return false; }

	return true;

} // end ray box intersection


// ray -> triangle without smoothed normals
template<typename T> static bool mcl::raycast::ray_triangle( const Ray<T> *ray,
	const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2, Payload<T> *payload ){

	// Compute hit point
	const Vec3<T> e0 = p1 - p0;
	const Vec3<T> e1 = p0 - p2;
	Vec3<T> n = e1.cross( e0 );
	const Vec3<T> e2 = ( 1.0 / n.dot( ray->direction ) ) * ( p0 - ray->origin );
	T t = n.dot( e2 );
	bool hit = (t<payload->t_max) && (t>payload->t_min);
	if( !hit ){ return false; }

	// Compute bary coords
	const Vec3<T> i  = ray->direction.cross( e2 );
	T beta  = i.dot( e1 );
	T gamma = i.dot( e0 );
	T alpha = 1.0 - beta - gamma;
	bool bary_test = alpha>0 && beta>0 && gamma>0 && (alpha+beta+gamma)<=1;
	if( !bary_test ){ return false; }

	// Compute payload variables
	n.normalize();
	payload->n = n;
	payload->t_max = t;
	payload->hit_point = ray->origin + ray->direction*t;
	payload->bary = Vec3<T>(alpha,beta,gamma);
	return true;

} // end  ray -> triangle


// ray -> triangle with smoothed normals
template<typename T> static bool mcl::raycast::ray_triangle( const Ray<T> *ray,
	const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2,
	const Vec3<T> &n0, const Vec3<T> &n1, const Vec3<T> &n2, Payload<T> *payload ){

	if( ray_triangle<T>( ray, p0, p1, p2, payload ) ){
		payload->n = payload->bary[0]*n0 + payload->bary[1]*n1 + payload->bary[2]*n2;
		return true;
	}

	return false;

} // end  ray -> triangle


// ray -> sphere
template<typename T> static bool mcl::raycast::ray_sphere( const Ray<T> *ray,
	const Vec3<T> &center, const T &radius, Payload<T> *payload ){

	const Vec3<T> s = ray->origin - center;
	T a = ray->direction.dot( ray->direction );
	T b = ( (2*s).dot( ray->direction ) );
	T c = s.dot(s)-( radius*radius );
	T disc = b*b - 4*a*c;
	if( disc < 0 ){ return false; }

	disc = std::sqrt( disc );
	T q = 0;
	if( b < 0 ){ q = ( -b - disc )/2; }
	else{ q = ( -b + disc )/2; }
	T t = ( -0.5 )*( b+disc )/a;
	bool hit = (t<payload->t_max) && (t>payload->t_min);
	if( !hit ){ return false; }

	// Set payload
	payload->t_max = t;
	payload->hit_point = ray->origin + ray->direction*t;
	payload->n = ( payload->hit_point - center );
	payload->n.normalize();
	return true;

} // end ray sphere intersection


} // end namespace mcl

#endif
