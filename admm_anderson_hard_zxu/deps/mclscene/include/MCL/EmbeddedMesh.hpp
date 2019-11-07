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

#ifndef MCL_EMBEDDEDMESH_H
#define MCL_EMBEDDEDMESH_H

#include "TetMesh.hpp"
#include "TriangleMesh.hpp"
#include "Projection.hpp"

namespace mcl {

class EmbeddedMesh {
public:
	typedef std::shared_ptr<EmbeddedMesh> Ptr;
	static std::shared_ptr<EmbeddedMesh> create(){
		return std::make_shared<EmbeddedMesh>();
	}

	EmbeddedMesh(){
		embedded = TriangleMesh::create();
		lattice = TetMesh::create();
	}

	// Copy constructor makes copies (and new pointers) of inner meshes
	EmbeddedMesh( const EmbeddedMesh &mesh );

	int flags;
	std::shared_ptr<TriangleMesh> embedded; // the embedded surface-only mesh
	std::shared_ptr<TetMesh> lattice; // tet mesh that embeds the surface
	std::vector<Vec4f> barycoords; // per emb vert bary coords
	std::vector<int> vert_to_tet; // mapping from emb vert to tet idx

	// Updates the positions/normals of the embedded vertices
	// after a change in the lattice.
	inline void update_embedded();

	// Computes barycoords and vert_to_tet by mapping embedded
	// vertices into the lattice. Assumes a lattice has already been created.
	// Unlike update_embedded, you don't need to call this more than once.
	inline bool update_lattice();

	// Generates a lattice around the embedded triangle mesh.
	// If tets/verts already exist in the lattice, they are removed and re-generated.
	// To simply compute barycoords/vert_to_tet, use update_lattice().
	// Tess is the approx number of cubes to generate on the largest face.
	// TODO: Better lattice gen, I do it poorly at the moment
	inline bool gen_lattice( int tess=7 );

	template<typename T> void apply_xform( const XForm<T,3> &xf );

	// Returns aabb
	inline Eigen::AlignedBox<float,3> bounds();

	// Computes volume-weighted masses for each lattice vertex.
	// density_kgm3 is the unit-volume density (e.g. soft rubber: 1100)
	// TODO Exact mass computation from
	// "Fast Viscoelastic Behavior with Thin Features" (2008) by Wojtan and Turk.
	inline void weighted_masses( std::vector<float> &m, float density_kgm3=1100.0 );

private:
	static inline void gen_tets( Vec3f min, Vec3f max, std::vector<Vec3f> &verts, std::vector<Vec4i> &tets );

	static inline float baryweight( short i, const Vec4f &bary ){ return bary[i] / ( bary.dot(bary) ); }

}; // end class EmbeddedMesh

EmbeddedMesh::EmbeddedMesh( const EmbeddedMesh &mesh ) : flags(0) {
	embedded = std::make_shared<TriangleMesh>( *(mesh.embedded) );
	lattice = std::make_shared<TetMesh>( *(mesh.lattice) );
	barycoords = mesh.barycoords;
	vert_to_tet = mesh.vert_to_tet;
}

inline void EmbeddedMesh::update_embedded(){

	const size_t nv = vert_to_tet.size();
	const int nt = lattice->tets.size();
	if( nv != embedded->vertices.size() || nv != barycoords.size() ){
		std::cerr << "**EmbeddedMesh::update_embedded Error: was there a topology change?" << std::endl;
		return;
	}

	#pragma omp parallel for
	for( size_t i=0; i<nv; ++i ){
		int t_idx = vert_to_tet[i];
		if( t_idx < 0 || t_idx >= nt ){ continue; }
		const Vec4i &tet = lattice->tets[t_idx];
		const Vec4f &barys = barycoords[i];
		embedded->vertices[i] = barys[0]*lattice->vertices[tet[0]] +
			barys[1]*lattice->vertices[tet[1]] +
			barys[2]*lattice->vertices[tet[2]] +
			barys[3]*lattice->vertices[tet[3]];
	} // end loop verts

	embedded->need_normals(true);

} // end update embedded


inline bool EmbeddedMesh::gen_lattice( int tess ){
	typedef Eigen::AlignedBox<float,3> AABB;
	lattice->clear();

	// TODO THIS WHOLE THING

	AABB aabb;
	const int nv = embedded->vertices.size();
	for( int i=0; i<nv; ++i ){ aabb.extend( embedded->vertices[i] ); }


	float step_scalar = 1.f/float(tess);
	Vec3f min = aabb.min();
	Vec3f max = aabb.max();
	Vec3f step = (max-min)*step_scalar;
	float step_coeff = step.maxCoeff();
	step[0] = step_coeff;
	step[1] = step_coeff;
	step[2] = step_coeff;
	min -= step_scalar*step;
	max += step_scalar*step;

	std::vector<Vec3f> verts;
	std::vector<Vec4i> tets;

	// Generate the lattice
	for( float x=min[0]-step[0]; x<max[0]; x += step[0] )
	for( float y=min[1]-step[1]; y<max[1]; y += step[1] )
	for( float z=min[2]-step[2]; z<max[2]; z += step[2] ){
		Vec3f lower(x,y,z);
		Vec3f upper = lower+step;
		gen_tets( lower, upper, verts, tets );
	} // end loop grid

	// Compute bary coords and remove any tet that doesn't contain a vertex
	barycoords.resize( nv, Vec4f(-1,-1,-1,-1) ); // barys within the tet
	vert_to_tet.resize( nv, -1 ); // vertex to tet mapping
	std::vector<int> num_v_in_t( tets.size(), 0 ); // num verts in a tet

	omp_lock_t writelock;
	omp_init_lock(&writelock);

	#pragma omp parallel for
	for( int i=0; i<nv; ++i ){
		Vec3f point = embedded->vertices[i];
		const int nt = tets.size();
		for( int j=0; j<nt; ++j ){ // TODO accel structure
			Vec4i tet = tets[j];
			if( !mcl::projection::point_in_tet( point,
				verts[tet[0]], verts[tet[1]], verts[tet[2]], verts[tet[3]] ) ){ continue; }
			barycoords[i] = mcl::vec::barycoords( point,
				verts[tet[0]], verts[tet[1]], verts[tet[2]], verts[tet[3]] );
			vert_to_tet[i] = j;
		}
		if( vert_to_tet[i] < 0 ){ continue; }

		omp_set_lock(&writelock);
		num_v_in_t[ vert_to_tet[i] ]++;
		omp_unset_lock(&writelock);
	}

	omp_destroy_lock(&writelock);

	// Remove any tets that do not contain vertices
	int n_tets = tets.size();
	lattice->clear();
	lattice->vertices = verts;
	lattice->tets.reserve(n_tets);
	for( int i=0; i<n_tets; ++i ){
		if( num_v_in_t[i] > 0 ){
			lattice->tets.emplace_back( tets[i] );
		}
	}

	verts.clear();
	tets.clear();

	// remove vertices that aren't a part of a tet
	// and combine ones that are close
	lattice->refine();

	// Update embedded verts again
	#pragma omp parallel for
	for( int i=0; i<nv; ++i ){
		Vec3f point = embedded->vertices[i];
		const int nt = lattice->tets.size();
		for( int j=0; j<nt; ++j ){ // TODO accel structure
			Vec4i tet = lattice->tets[j];
			if( !mcl::projection::point_in_tet( point,
				lattice->vertices[tet[0]],
				lattice->vertices[tet[1]],
				lattice->vertices[tet[2]],
				lattice->vertices[tet[3]] ) ){ continue; }
			barycoords[i] = mcl::vec::barycoords( point,
				lattice->vertices[tet[0]],
				lattice->vertices[tet[1]],
				lattice->vertices[tet[2]],
				lattice->vertices[tet[3]] );
			if( barycoords[i].minCoeff() < 0 || barycoords[i].sum()-1e-3 > 1 ){
				std::stringstream ss;
				ss << "BAD BARYS: tet(" << vert_to_tet[i] << "), sum(" << barycoords[i].sum() <<
					"), coords(" << barycoords[i].transpose() << ")" << std::endl;
				printf("%s", ss.str().c_str() );
				vert_to_tet[i] = -1;
				continue;
			}
			vert_to_tet[i] = j;
		}
	}

	// Double check values
	n_tets = lattice->tets.size();
	for( int i=0; i<nv; ++i ){
		int v2t = vert_to_tet[i];
		if( v2t < 0 || v2t >= n_tets ){
			std::cerr << "Emb vert " << i << " has tet index " << v2t << ", with " << n_tets << " tets" << std::endl;
			return false;
		}
	}

	return true;

} // end gen lattice

template<typename T>
void EmbeddedMesh::apply_xform( const XForm<T,3> &xf ){
	if( lattice->vertices.size() == 0 ){
		std::cerr << "EmbeddedMesh::apply_xform Error: XForm is applied to lattice, which hasn't been set" << std::endl;
		return;
	}
	lattice->apply_xform( xf );
	update_embedded();
}


inline Eigen::AlignedBox<float,3> EmbeddedMesh::bounds() {
	size_t n_lat_verts = lattice->vertices.size();
	size_t n_emb_verts = embedded->vertices.size();
	Eigen::AlignedBox<float,3> aabb;
	if( n_lat_verts == 0 ){
		for( size_t i=0; i<n_emb_verts; ++i ){ aabb.extend( embedded->vertices[i] ); }
	}
	else {
		for( size_t i=0; i<n_lat_verts; ++i ){ aabb.extend( lattice->vertices[i] ); }
	}
	return aabb;
}


inline void EmbeddedMesh::weighted_masses( std::vector<float> &m, float density_kgm3 ){

	if( !lattice ){ return; }
	m.resize( lattice->vertices.size(), 0.f );
	int n_tets = lattice->tets.size();
	for( int t=0; t<n_tets; ++t ){
		Vec4i tet = lattice->tets[t];
		Eigen::Matrix<float,3,3> edges;
		edges.col(0) = lattice->vertices[tet[1]] - lattice->vertices[tet[0]];
		edges.col(1) = lattice->vertices[tet[2]] - lattice->vertices[tet[0]];
		edges.col(2) = lattice->vertices[tet[3]] - lattice->vertices[tet[0]];
		float v = std::abs( (edges).determinant()/6.f );
		float tet_mass = density_kgm3 * v;
		m[ tet[0] ] += tet_mass / 4.f;
		m[ tet[1] ] += tet_mass / 4.f;
		m[ tet[2] ] += tet_mass / 4.f;
		m[ tet[3] ] += tet_mass / 4.f;
	}

} // end weighted masses


inline void EmbeddedMesh::gen_tets( Vec3f min, Vec3f max, std::vector<Vec3f> &verts, std::vector<Vec4i> &tets ){

	// Top plane, clockwise looking down
	Vec3f a = max;
	Vec3f b( min[0], max[1], max[2] );
	Vec3f c( min[0], max[1], min[2] );
	Vec3f d( max[0], max[1], min[2] );

	// Bottom plan, clockwise looking down
	Vec3f e( max[0], min[1], max[2] );
	Vec3f f( min[0], min[1], max[2] );
	Vec3f g( min[0], min[1], min[2] );
	Vec3f h( max[0], min[1], min[2] );

	// Add the verts
	int nv = verts.size();
	verts.emplace_back( a ); // 0
	verts.emplace_back( b ); // 1
	verts.emplace_back( c ); // 2
	verts.emplace_back( d ); // 3
	verts.emplace_back( e ); // 4
	verts.emplace_back( f ); // 5
	verts.emplace_back( g ); // 6
	verts.emplace_back( h ); // 7

	// Pack 5 tets into the cube
	Vec4i t0( 0, 5, 7, 4 );
	Vec4i t1( 5, 7, 2, 0 );
	Vec4i t2( 5, 0, 2, 1 );
	Vec4i t3( 7, 2, 0, 3 );
	Vec4i t4( 5, 2, 7, 6 );
	Vec4i offset(nv,nv,nv,nv);

	// Add the tets
	tets.emplace_back( t0+offset );
	tets.emplace_back( t1+offset );
	tets.emplace_back( t2+offset );
	tets.emplace_back( t3+offset );
	tets.emplace_back( t4+offset );

} // end gen tets for a cube



inline bool EmbeddedMesh::update_lattice(){
throw std::runtime_error("**EmbeddedMesh TODO: update_lattice");
/*
	const int nv = embedded->vertices.size();
	vert_to_tet.resize( nv, -1 );
	barycoords.resize( nv, Vec4f(0,0,0,0) );
	std::vector<int> num_v_in_t( lattice->tets.size(), 0 ); // num verts in a tet

	omp_lock_t writelock;
	omp_init_lock(&writelock);

	#pragma omp parallel for
	for( int i=0; i<nv; ++i ){
		const Vec3f &point = embedded->vertices[i];
		const int nt = lattice->tets.size();
		for( int j=0; j<nt; ++j ){
			const Vec4i &tet = lattice->tets[j];
			if( !mcl::projection::point_in_tet( point,
				lattice->vertices[tet[0]],
				lattice->vertices[tet[1]],
				lattice->vertices[tet[2]],
				lattice->vertices[tet[3]] ) ){ continue; }
			barycoords[i] = mcl::vec::barycoords( point,
				lattice->vertices[tet[0]],
				lattice->vertices[tet[1]],
				lattice->vertices[tet[2]],
				lattice->vertices[tet[3]] );
			vert_to_tet[i] = j;
		}
		if( vert_to_tet[i] < 0 ){ continue; }

		omp_set_lock(&writelock);
		num_v_in_t[ vert_to_tet[i] ]++;
		omp_unset_lock(&writelock);
	}

	omp_destroy_lock(&writelock);
*/
	return true;

} // end compute barys and vert to tet



} // namespace mcl

#endif
