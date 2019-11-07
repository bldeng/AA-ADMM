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

#ifndef MCL_VISITOR_H
#define MCL_VISITOR_H 1

#include "Vec.hpp"
#include "Projection.hpp"
#include "Raycast.hpp"

namespace mcl {
namespace bvh {

// Visitor class for traversing the tree
// PDIM is the dimension of the primitive,
// i.e. 1=verts, 2=edges, 3=tris, 4=tets
template <typename T, short PDIM>
class Visitor {
typedef Eigen::AlignedBox<T,3> AABB;
public:
	// See if we hit the current node. Return true to check node's children.
	virtual bool hit_aabb( const AABB &aabb ) = 0;

	// See if we hit the current primitive. Return true to stop traversing.
	virtual bool hit_prim( int prim ) = 0;

	// Return true if the left node should be checked before the right one.
	virtual bool check_left_first( const AABB &left, const AABB &right ) = 0;
};


// Point in tet
template <typename T>
class PointInTet : public Visitor<T,4> {
typedef Eigen::AlignedBox<T,3> AABB;
public:
	Vec3<T> point; // query point
	int hit_tet; // intersected tet
	std::vector<int> skip_vert_idx;  // vert index to skip (for self collision)
	const T *verts;
	const int *inds;
	PointInTet( Vec3<T> point_, const T *verts_, const int *inds_ );
	bool hit_aabb( const AABB &aabb );
	bool hit_prim( int prim );
	bool check_left_first( const AABB &left, const AABB &right );
};


// Nearest point on surface
template <typename T>
class NearestTriangle : public Visitor<T,3> {
typedef Eigen::AlignedBox<T,3> AABB;
public:
	Vec3<T> point; // query point
	Vec3<T> proj; // nearest point on hit_tri
	int hit_tri; // triangle idx
	std::vector<int> skip_vert_idx; // vert index to skip (for self collision)
	T curr_nearest; // current nearest distance to tri
	const T *verts;
	const int *inds;
	NearestTriangle( Vec3<T> point_, const T *verts_, const int *inds_ );
	bool hit_aabb( const AABB &aabb );
	bool hit_prim( int prim );
	bool check_left_first( const AABB &left, const AABB &right );
};

/*
// Raycast with multi-hit (counter)
template <typename T>
class RayMultiHit : public Visitor<T,3> {
typedef Eigen::AlignedBox<T,3> AABB;
public:
	raycast::Ray<T> ray;
	int skip_vert_idx; // vert index to skip (for self collision)
	int hit_count; // number of intersections
	const T *verts;
	const int *inds;
	RayMultiHit( Vec3<T> point_, const T *verts_, const int *inds_ );
	bool hit_aabb( const AABB &aabb );
	bool hit_prim( int prim );
	bool check_left_first( const AABB &left, const AABB &right );
};
*/

//
//	Implementation
//

//
// PointInTet
//

template <typename T> 
PointInTet<T>::PointInTet( Vec3<T> point_, const T *verts_, const int *inds_ ) :
	point(point_), verts(verts_), inds(inds_), hit_tet(-1) {}

template <typename T> 
bool PointInTet<T>::hit_aabb( const AABB &aabb ){
	if( aabb.isEmpty() ){ throw std::runtime_error("PointInTet Error: Empty AABB"); }
	return aabb.contains(point);
}

template <typename T> 
bool PointInTet<T>::hit_prim( int prim ){
	Vec4i tet( inds[prim*4+0], inds[prim*4+1], inds[prim*4+2], inds[prim*4+3] );
	int n_skip = skip_vert_idx.size();
	for( int i=0; i<n_skip; ++i ){
		for( int j=0; j<4; ++j ){
			if( skip_vert_idx[i]==tet[j] ){ return false; }
		}
	}
	Vec3<T> v0( verts[tet[0]*3+0], verts[tet[0]*3+1], verts[tet[0]*3+2] );
	Vec3<T> v1( verts[tet[1]*3+0], verts[tet[1]*3+1], verts[tet[1]*3+2] );
	Vec3<T> v2( verts[tet[2]*3+0], verts[tet[2]*3+1], verts[tet[2]*3+2] );
	Vec3<T> v3( verts[tet[3]*3+0], verts[tet[3]*3+1], verts[tet[3]*3+2] );
	if( projection::point_in_tet<T>( point, v0, v1, v2, v3 ) ){
		hit_tet = prim;
		return true;
	}
	return false;
}

template <typename T> 
bool PointInTet<T>::check_left_first( const AABB &left, const AABB &right ){
	T left_ed = left.squaredExteriorDistance( point );
	return left_ed <= right.squaredExteriorDistance( point );
}

//
// NearestTriangle
//

template <typename T> 
NearestTriangle<T>::NearestTriangle( Vec3<T> point_, const T *verts_, const int *inds_ ) :
	point(point_), verts(verts_), inds(inds_), hit_tri(-1), proj(-1,-1,-1),
	curr_nearest(std::numeric_limits<T>::max()) {}

template <typename T> 
bool NearestTriangle<T>::hit_aabb( const AABB &aabb ){
	return aabb.squaredExteriorDistance(point) < curr_nearest;
}

template <typename T> 
bool NearestTriangle<T>::hit_prim( int prim ){
	Vec3i tri( inds[prim*3+0], inds[prim*3+1], inds[prim*3+2] );
	int n_skip = skip_vert_idx.size();
	for( int i=0; i<n_skip; ++i ){
		for( int j=0; j<3; ++j ){
			if( skip_vert_idx[i]==tri[j] ){ return false; }
		}
	}
	Vec3<T> v0( verts[tri[0]*3+0], verts[tri[0]*3+1], verts[tri[0]*3+2] );
	Vec3<T> v1( verts[tri[1]*3+0], verts[tri[1]*3+1], verts[tri[1]*3+2] );
	Vec3<T> v2( verts[tri[2]*3+0], verts[tri[2]*3+1], verts[tri[2]*3+2] );

	Vec3<T> p = projection::point_on_triangle( point, v0, v1, v2 );
	T dist = (p-point).squaredNorm();
	if( dist > curr_nearest ){ return false; }

	curr_nearest = dist;
	hit_tri = prim;
	proj = p;
	return false; // return false to keep checking other tris
}

template <typename T> 
bool NearestTriangle<T>::check_left_first( const AABB &left, const AABB &right ){
	T left_ed = left.squaredExteriorDistance( point );
	return left_ed < right.squaredExteriorDistance( point );
}

//
// RayMultiHit
//
/*
template <typename T> 
RayMultiHit<T>::RayMultiHit( Vec3<T> point_, const T *verts_, const int *inds_ ) :
	ray(point_, Vec3<T>(0,1,0)), verts(verts_), inds(inds_), skip_vert_idx(-1), hit_count(0) {}

template <typename T> 
bool RayMultiHit<T>::hit_aabb( const AABB &aabb ){
	return raycast::ray_aabb( &ray, aabb.min(), aabb.max() );
}

template <typename T> 
bool RayMultiHit<T>::hit_prim( int prim ){

	Vec3i tri( inds[prim*3+0], inds[prim*3+1], inds[prim*3+2] );
	for( int i=0; i<3; ++i ){ if( tri[i]==skip_vert_idx ){ return false; } }
	Vec3<T> v0( verts[tri[0]*3+0], verts[tri[0]*3+1], verts[tri[0]*3+2] );
	Vec3<T> v1( verts[tri[1]*3+0], verts[tri[1]*3+1], verts[tri[1]*3+2] );
	Vec3<T> v2( verts[tri[2]*3+0], verts[tri[2]*3+1], verts[tri[2]*3+2] );

	// We'll use a fresh payload each time so we don't only find
	// nearest intersections. This lets us count how many times we
	// hit a surface, and sort them from nearest to furthest.
	raycast::Payload<T> tmpPayload;
	bool hit = raycast::ray_triangle( &ray, v0, v1, v2, &tmpPayload );
	if( hit ){ hit_count++; }
	return false; // return false to keep checking other tris
}

template <typename T> 
bool RayMultiHit<T>::check_left_first( const AABB &left, const AABB &right ){
	return raycast::ray_aabb( &ray, left.min(), left.max() );
}
*/
//
// End
//

} // end ns bvh
} // end ns mcl

#endif
