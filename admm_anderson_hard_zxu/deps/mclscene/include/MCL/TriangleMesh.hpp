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

#ifndef MCL_TRIANGLEMESH_H
#define MCL_TRIANGLEMESH_H 1

#include <vector>
#include <memory>
#include "Vec.hpp"
#include "XForm.hpp"
#include "HashKeys.hpp"

namespace mcl {

class TriangleMesh {
public:
	typedef std::shared_ptr<TriangleMesh> Ptr;
	static std::shared_ptr<TriangleMesh> create(){
		return std::make_shared<TriangleMesh>();
	}

	TriangleMesh() : flags(0) {}

	// Data
	int flags;
	std::vector< Vec3f > vertices; // all vertices in the mesh
	std::vector< Vec3f > normals; // zero length for all non-surface normals
	std::vector< Vec3i > faces; // surface triangles
	std::vector< Vec2f > texcoords; // per vertex uv coords
	std::vector< Vec2i > edges; // unique face edges
	std::vector< Vec2i > exterior_edges; // edges on the boundary only

	// Get per-vertex data.
	// If normals have not been set, they are computed.
	inline void get_vertex_data(
		float* &vertices, int &num_vertices,
		float* &normals, int &num_normals,
		float* &texcoords, int &num_texcoords
	);

	// Get primitive data.
	// If edges are requested but have not been set, they are computed.
	// Dimension describes the prim type, i.e. 2 = edges, 3 = triangles, etc...
	inline void get_primitive_data( short dimension, int* &prims, int &num_prims );

	template<typename T> void apply_xform( const XForm<T,3> &xf );

	// Returns aabb
	inline Eigen::AlignedBox<float,3> bounds();

	// Computes per-vertex normals
	inline void need_normals( bool recompute=false );

	// Creates unique edges of the triangle faces
	inline void need_edges( bool recompute=false );

	// Creates edges as above, but only on the exterior surface
	inline void need_exterior_edges( bool recompute=false );

	// Removes vertices not indexed by a triangle, and combines vertices
	// that are within eps distance (lowest index is kept).
	inline void refine( float eps=1e-6f );

	// Computes area-weighted masses for each vertex.
	// density_kgm2 is the density per unit area.
	// Most cloth, for instance, is like 0.1 to 0.6.
	inline void weighted_masses( std::vector<float> &m, float density_kgm2=0.4f );

	// Clear all mesh data
	inline void clear();

}; // end class TriangleMesh


//
//	Implementation
//


inline void TriangleMesh::get_vertex_data(
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


inline void TriangleMesh::get_primitive_data( short dim, int* &prims, int &num_prims ){
	if( dim == 2 ){
		need_edges(); // compute edges if we don't have them
		num_prims = edges.size();
		if( num_prims > 0 ){ prims = &edges[0][0]; }
	}
	else if( dim == 3 && faces.size() > 0 ){
		num_prims = faces.size();
		prims = &faces[0][0];
	}
} // end get prim data


template<typename T>
void TriangleMesh::apply_xform( const XForm<T,3> &xf_ ){
	Eigen::Transform<float,3,Eigen::Affine> xf = xf_.template cast<float>();
	int nv = vertices.size();
	for(int i=0; i<nv; ++i){ vertices[i] = xf * vertices[i]; }
} // end apply xform


inline Eigen::AlignedBox<float,3> TriangleMesh::bounds(){
	size_t n_verts = vertices.size();
	Eigen::AlignedBox<float,3> aabb;
	for( size_t i=0; i<n_verts; ++i ){ aabb.extend( vertices[i] ); }
	return aabb;
}

inline void TriangleMesh::need_normals( bool recompute ){
	const size_t nv = vertices.size();
	if( nv == normals.size() && !recompute ){ return; }
	if( nv != normals.size() ){ normals.resize( nv ); }
	std::fill( normals.begin(), normals.end(), Vec3f(0,0,0) );
	int nf = faces.size();
	for( int i = 0; i < nf; ++i ){
		const Vec3f &p0 = vertices[faces[i][0]];
		const Vec3f &p1 = vertices[faces[i][1]];
		const Vec3f &p2 = vertices[faces[i][2]];
		Vec3f a = p0-p1, b = p1-p2, c = p2-p0;
		float l2a = a.squaredNorm(), l2b = b.squaredNorm(), l2c = c.squaredNorm();
		if(!l2a || !l2b || !l2c){ continue; }
		Vec3f facenormal = a.cross( b );
		normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
		normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
		normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
	}
	for(size_t i = 0; i < nv; ++i){
		if( normals[i].squaredNorm() > 0 ){ normals[i].normalize(); }
	}
} // end compute normals


inline void TriangleMesh::need_edges( bool recompute ){

	if( edges.size()>0 && !recompute ){ return; }

	// vertex ids -> number of faces using these indices
	std::unordered_map< hashkey::sint2, int > edge_ids;
	int n_faces = faces.size();
	for( int f=0; f<n_faces; ++f ){
		edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][0],faces[f][1]), 1) );
		edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][0],faces[f][2]), 1) );
		edge_ids.emplace( std::make_pair( hashkey::sint2(faces[f][1],faces[f][2]), 1) );
	}

	// Now copy the map into edges
	edges.clear();
	std::unordered_map< hashkey::sint2, int >::iterator it = edge_ids.begin();
	for( ; it != edge_ids.end(); ++it ){
		edges.emplace_back( Vec2i(it->first[0],it->first[1]) );
	}

} // end compute edges

// Loops through all of the triangles and counts the number of
// times an edge was indexed. If once, it's a surface edge.
inline void TriangleMesh::need_exterior_edges( bool recompute ){
	if( edges.size()>0 && !recompute ){ return; }

	std::unordered_map< hashkey::sint2, int > edge_ids;
	size_t n_faces = faces.size();
	for( size_t t=0; t<n_faces; ++t ){
		int p0 = faces[t][0];
		int p1 = faces[t][1];
		int p2 = faces[t][2];
		hashkey::sint2 curr_edges[3];
		curr_edges[0] = hashkey::sint2( p0, p1 );
		curr_edges[1] = hashkey::sint2( p0, p2 );
		curr_edges[2] = hashkey::sint2( p1, p2 );
		for( int f=0; f<3; ++f ){
			if( edge_ids.count(curr_edges[f]) == 0 ){ edge_ids[ curr_edges[f] ] = 1; }
			else{ edge_ids[ curr_edges[f] ] += 1; }
		}
	}
	exterior_edges.clear();
	std::unordered_map< hashkey::sint2, int >::iterator edgeIt = edge_ids.begin();
	for( ; edgeIt != edge_ids.end(); ++edgeIt ){
		if( edgeIt->second == 1 ){
			hashkey::sint2 f = edgeIt->first;
			exterior_edges.emplace_back( Vec2i( f.orig_v[0], f.orig_v[1] ) );
		}
	}

} // end compute edges



inline void TriangleMesh::refine( float eps ){

	double eps2 = double(eps)*double(eps); // so we can use squaredNorm
	int n_faces = faces.size();
	int n_verts_0 = vertices.size();
	std::vector<int> vert_refs( n_verts_0, 0 );

	for( int i=0; i<n_faces; ++i ){
		Vec3i face = faces[i];
		for( int j=0; j<3; ++j ){
			int idx = face[j];
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

			// Update face index
			faces[i][j] = new_idx;
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

	// Update face indices
	for( int i=0; i<n_faces; ++i ){
		Vec3i face = faces[i];
		for( int j=0; j<3; ++j ){
			int idx = face[j];
			if( old_to_new.count(idx)==0 ){
				throw std::runtime_error("TetMesh::refine Error: Something went wrong.");
			}
			faces[i][j] = old_to_new[idx];
		}
	}

	// Remake other data if needed
	if( edges.size() ){ need_edges(true); }
	if( normals.size() ){ need_normals(true); }

} // end refine

inline void TriangleMesh::weighted_masses( std::vector<float> &m, float density_kgm2 ){

	m.resize( vertices.size(), 0.f );
	int n_faces = faces.size();
	for( int f=0; f<n_faces; ++f ){
		Vec3i face = faces[f];
		Vec3f edge1 = vertices[ face[1] ] - vertices[ face[0] ];
		Vec3f edge2 = vertices[ face[2] ] - vertices[ face[0] ];
		float area = 0.5f * (edge1.cross(edge2)).norm();
		float tri_mass = density_kgm2 * area;
		m[ face[0] ] += tri_mass / 3.f;
		m[ face[1] ] += tri_mass / 3.f;
		m[ face[2] ] += tri_mass / 3.f;
	}

} // end weighted masses


inline void TriangleMesh::clear(){
	vertices.clear();
	normals.clear();
	faces.clear();
	texcoords.clear();
	edges.clear();
} // end clear all data

} // end namespace mcl

#endif
