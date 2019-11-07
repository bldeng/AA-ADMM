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

#ifndef MCL_TETMESH_H
#define MCL_TETMESH_H

#include <vector>
#include <memory>
#include "XForm.hpp"
#include "HashKeys.hpp"
#include <iostream>

namespace mcl {

class TetMesh {
public:
	typedef std::shared_ptr<TetMesh> Ptr;
	static std::shared_ptr<TetMesh> create(){
		return std::make_shared<TetMesh>();
	}

	TetMesh() : flags(0) {}

	// Data
	int flags;
	std::vector< Vec4i > tets; // all elements
	std::vector< Vec3f > vertices; // all vertices in the mesh
	std::vector< Vec3f > normals; // zero length for all non-surface normals
	std::vector< Vec3i > faces; // surface triangles
	std::vector< Vec2f > texcoords; // per vertex uv coords
	std::vector< Vec2i > edges; // unique tet edges

	// Get per-vertex data.
	// If normals have not been set, they are computed.
	inline void get_vertex_data(
		float* &vertices, int &num_vertices,
		float* &normals, int &num_normals,
		float* &texcoords, int &num_texcoords
	);

	// Get primitive data.
	// If edges/faces are requested but have not been set, they are computed.
	// Dimension describes the prim type, i.e. 2 = edges, 3 = triangles, 4 = tets, etc...
	inline void get_primitive_data( short dimension, int* &prims, int &num_prims );

	template<typename T> void apply_xform( const XForm<T,3> &xf );

	// Returns AABB
	inline Eigen::AlignedBox<float,3> bounds();

	// Finds and stores the surface trianges
	inline void need_faces( bool recompute=false );

	// Computes per-vertex normals
	inline void need_normals( bool recompute=false );

	// Creates unique edges of the tets or faces
	inline void need_edges( bool recompute=false, bool surface_only=true );

	// Removes vertices not indexed by a tet, and combines vertices
	// that are within eps distance (lowest index is kept).
	inline void refine( float eps=1e-6f );

	// Computes volume-weighted masses for each vertex
	// density_kgm3 is the unit-volume density (e.g. soft rubber: 1100)
	// See: https://www.engineeringtoolbox.com/density-solids-d_1265.html
	inline void weighted_masses( std::vector<float> &m, float density_kgm3=1100.0 );

	// Returns a list of vertex indices that are on the surface
	inline void surface_inds( std::vector<int> &surf_inds );

	// Clear all mesh data
	inline void clear();

}; // end class TetMesh


//
//	Implementation
//


inline void TetMesh::get_vertex_data(
	float* &verts, int &num_vertices,
	float* &norms, int &num_normals,
	float* &tex, int &num_texcoords ){
	if( normals.size() != vertices.size() ){ need_normals(); }
	num_vertices = vertices.size();
	num_normals = normals.size();
	num_texcoords = texcoords.size();
	if( num_vertices > 0 ){ verts = &vertices[0][0]; }
	if( num_normals > 0 ){ norms = &normals[0][0]; }
	if( num_texcoords > 0 ){ tex = &texcoords[0][0]; }

} // end get vertex data


inline void TetMesh::get_primitive_data( short dim, int* &prims, int &num_prims ){
	if( dim == 2 ){
		need_edges(); // compute edges if we don't have them
		num_prims = edges.size();
		if( num_prims > 0 ){ prims = &edges[0][0]; }
	}
	else if( dim == 3 && faces.size() > 0 ){
		num_prims = faces.size();
		prims = &faces[0][0];
	}
	else if( dim == 4 && tets.size() > 0 ){
		num_prims = tets.size();
		prims = &tets[0][0];
	}
} // end get prim data


template<typename T>
void TetMesh::apply_xform( const XForm<T,3> &xf_ ){
	Eigen::Transform<float,3,Eigen::Affine> xf = xf_.template cast<float>();
	int nv = vertices.size();
	for(int i=0; i<nv; ++i){ vertices[i] = xf * vertices[i]; }
} // end apply xform


inline Eigen::AlignedBox<float,3> TetMesh::bounds(){
	size_t n_verts = vertices.size();
	Eigen::AlignedBox<float,3> aabb;
	for( size_t i=0; i<n_verts; ++i ){ aabb.extend( vertices[i] ); }
	return aabb;
}


// Loop through all of the tets and counts the number of times a face is indexed,
// with the indices of the face sorted from low to high. If a face only exists
// on a tet once, it's an outer facing face.
inline void TetMesh::need_faces( bool recompute ){
	if( faces.size()>0 && !recompute ){ return; }
	std::unordered_map< hashkey::sint3, int > face_ids;
	size_t n_tets = tets.size();
	for( size_t t=0; t<n_tets; ++t ){
		int p0 = tets[t][0];
		int p1 = tets[t][1];
		int p2 = tets[t][2];
		int p3 = tets[t][3];
		hashkey::sint3 curr_faces[4];
		curr_faces[0] = hashkey::sint3( p0, p1, p3 );
		curr_faces[1] = hashkey::sint3( p0, p2, p1 );
		curr_faces[2] = hashkey::sint3( p0, p3, p2 );
		curr_faces[3] = hashkey::sint3( p1, p2, p3 );
		for( int f=0; f<4; ++f ){
			if( face_ids.count(curr_faces[f]) == 0 ){ face_ids[ curr_faces[f] ] = 1; }
			else{ face_ids[ curr_faces[f] ] += 1; }
		}
	}
	faces.clear();
	std::unordered_map< hashkey::sint3, int >::iterator faceIt = face_ids.begin();
	for( ; faceIt != face_ids.end(); ++faceIt ){
		if( faceIt->second == 1 ){
			hashkey::sint3 f = faceIt->first;
			faces.emplace_back( Vec3i( f.orig_v[0], f.orig_v[1], f.orig_v[2] ) );
		}
	}
} // end need faces


inline void TetMesh::need_normals( bool recompute ){
	const size_t nv = vertices.size();
	if( nv == normals.size() && !recompute ){ return; }
	if( nv != normals.size() ){ normals.resize( vertices.size() ); }
	std::fill( normals.begin(), normals.end(), Vec3f(0,0,0) );
	size_t nf = faces.size();
	if( nf == 0 ){ need_faces(); }
	for( size_t i = 0; i < nf; ++i ){
		const Vec3f &p0 = vertices[faces[i][0]];
		const Vec3f &p1 = vertices[faces[i][1]];
		const Vec3f &p2 = vertices[faces[i][2]];
		Vec3f a = p0-p1, b = p1-p2, c = p2-p0;
		float l2a = a.squaredNorm(), l2b = b.squaredNorm(), l2c = c.squaredNorm();
		if (!l2a || !l2b || !l2c){ continue; }
		Vec3f facenormal = a.cross( b );
		normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
		normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
		normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
	}
	for(size_t i=0; i<nv; ++i){
		if( normals[i].squaredNorm() > 0 ){ normals[i].normalize(); }
	}
} // end compute normals


inline void TetMesh::need_edges( bool recompute, bool surface_only ){
	if( edges.size()>0 && !recompute ){ return; }
	std::unordered_map< hashkey::sint2, int > edge_ids;
	if( surface_only ){
		need_faces();
		int n_faces = faces.size();
		for( int f=0; f<n_faces; ++f ){
			edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][0],faces[f][1]), 1) );
			edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][0],faces[f][2]), 1) );
			edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][1],faces[f][2]), 1) );
		}
	} else {
		int n_tets = tets.size();
		for( int f=0; f<n_tets; ++f )
		for( int i=0; i<4; ++i )
		for( int j=0; j<4; ++j ){
			if( i==j ){ continue; }
			edge_ids.emplace( std::make_pair( hashkey::sint2(tets[f][i],tets[f][j]), 1) );
		}
	}
	edges.clear();
	std::unordered_map< hashkey::sint2, int >::iterator it = edge_ids.begin();
	for( ; it != edge_ids.end(); ++it ){
		edges.emplace_back( Vec2i(it->first[0],it->first[1]) );
	}
} // end compute edges


inline void TetMesh::refine( float eps ){

	double eps2 = double(eps)*double(eps); // so we can use squaredNorm
	int n_tets = tets.size();
	int n_verts_0 = vertices.size();
	std::vector<int> vert_refs( n_verts_0, 0 );

	for( int i=0; i<n_tets; ++i ){
		Vec4i tet = tets[i];
		for( int j=0; j<4; ++j ){
			int idx = tet[j];
			int new_idx = idx;
			Vec3d curr_x = vertices[idx].cast<double>();

			// Loop through the vertices.
			// Keep the lowest index if two vertices are within eps
			for( int k=0; k<n_verts_0; ++k ){
				double d2 = (vertices[k].cast<double>()-curr_x).squaredNorm();
				if( d2 <= eps2 && k < new_idx ){
					new_idx = k;
					break;
				}
			}

			// Update tet index
			tets[i][j] = new_idx;
			vert_refs[ new_idx ] += 1;
		}
	} // end loop tets

	// Now make a list of new vertices
	std::unordered_map<int,int> old_to_new;
	std::vector<Vec3f> old_vertices = vertices;
	vertices.clear();
	int num_ref_verts = 0;
	for( int i=0; i<n_verts_0; ++i ){
		if( vert_refs[i] > 0 ){
			old_to_new.insert({i,num_ref_verts});
			vertices.emplace_back( old_vertices[i] );
			num_ref_verts++;
		}
	}

	// Update tet indices
	for( int i=0; i<n_tets; ++i ){
		Vec4i tet = tets[i];
		for( int j=0; j<4; ++j ){
			int idx = tet[j];
			if( old_to_new.count(idx)==0 ){
				throw std::runtime_error("TetMesh::refine Error: Something went wrong.");
			}
			tets[i][j] = old_to_new[idx];
		}
	}

	// Remake other data if needed
	if( faces.size() ){ need_faces(true); }
	if( edges.size() ){ need_edges(true); }
	if( normals.size() ){ need_normals(true); }

} // end refine

inline void TetMesh::weighted_masses( std::vector<float> &m, float density_kgm3 ){

	m.resize( vertices.size(), 0.f );
	int n_tets = tets.size();
	for( int t=0; t<n_tets; ++t ){
		Vec4i tet = tets[t];
		Eigen::Matrix<float,3,3> edges;
		edges.col(0) = vertices[tet[1]] - vertices[tet[0]];
		edges.col(1) = vertices[tet[2]] - vertices[tet[0]];
		edges.col(2) = vertices[tet[3]] - vertices[tet[0]];
		float v = std::abs( (edges).determinant()/6.f );
		float tet_mass = density_kgm3 * v;
		m[ tet[0] ] += tet_mass / 4.f;
		m[ tet[1] ] += tet_mass / 4.f;
		m[ tet[2] ] += tet_mass / 4.f;
		m[ tet[3] ] += tet_mass / 4.f;
	}

} // end weighted masses

inline void TetMesh::surface_inds( std::vector<int> &surf_inds ){
	bool had_faces = true;
	if( faces.size()==0 ){
		had_faces = false;
		need_faces();
	}

	// Get a list of indices (unique)
	std::unordered_map<int,int> ind_map;
	int n_faces = faces.size();
	for( int i=0; i<n_faces; ++i ){
		ind_map[ faces[i][0] ] = 1;
		ind_map[ faces[i][1] ] = 1;
		ind_map[ faces[i][2] ] = 1;
	}

	// Copy map to vector
	std::unordered_map<int,int>::iterator it = ind_map.begin();
	for( ; it != ind_map.end(); ++it ){
		surf_inds.emplace_back( it->first );
	}

	if( !had_faces ){
		faces.clear();
	}
}

inline void TetMesh::clear(){
	tets.clear();
	vertices.clear();
	normals.clear();
	faces.clear();
	texcoords.clear();
	edges.clear();
} // end clear all data


} // end namespace mcl

#endif
