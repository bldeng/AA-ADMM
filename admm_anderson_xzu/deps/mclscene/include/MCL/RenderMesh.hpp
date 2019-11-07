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
//	A render mesh helps with interop between mclscene and OpenGL.
//	TODO:
//	- double precision base mesh
//	- colors/materials
//	- subdivide
//	- make flat
//

#ifndef MCL_RENDERMESH_H
#define MCL_RENDERMESH_H 1

#include "Texture.hpp"
#include "TriangleMesh.hpp"
#include "TetMesh.hpp"
#include "Material.hpp"
#include <memory>

namespace mcl {

class RenderMesh {
public:
	typedef std::shared_ptr<RenderMesh> Ptr;

	static inline std::shared_ptr<RenderMesh> create(
		std::shared_ptr<TriangleMesh> mesh, int options=DEFAULT ){
		return std::make_shared<RenderMesh>( mesh, options );
	}

	static inline std::shared_ptr<RenderMesh> create(
		std::shared_ptr<TetMesh> mesh, int options=DEFAULT ){
		return std::make_shared<RenderMesh>( mesh, options );
	}

	enum {
		DEFAULT = 1 << 0, // regular triangle mesh
		WIREFRAME = 1 << 1, // draw as wireframe
		INVISIBLE = 1 << 2, // don't draw mesh
		DYNAMIC = 1 << 3, // The verts they are a-changin'
	};
	int flags;

	RenderMesh( std::shared_ptr<TriangleMesh> mesh, int options=DEFAULT );
	RenderMesh( std::shared_ptr<TetMesh> mesh, int options=DEFAULT );

	// Copy vertex data to GPU. If the pointer'd mesh changes,
	// it's up to you to call this function.
	// Reload can be any number of the RELOAD enum if we are updating.
	enum { // load flags:
		ALL = 1 << 0, // reload all
		VERTICES = 1 << 1, // reload vertices
		NORMALS = 1 << 2, // reload norms
		COLORS = 1 << 3, // reload colors
		PRIMS = 1 << 4, // reload primitives
	};
	inline void load_buffers( int load=ALL );

	// Get the model matrix.
	inline const mcl::XForm<float> &get_model() const { return model; }

	// Returns aabb for the render mesh
	inline Eigen::AlignedBox<float,3> bounds();

	// Draws the mesh with current settings.
	inline void draw();

	// Vertex data. These are pointers to the actual
	// data stored in a tet/tri mesh.
	float *vertices, *normals, *colors, *texcoords;
	int num_vertices, num_normals, num_colors, num_texcoords;

	// Primitive data
	int *prims;
	int num_prims;

	// OpenGL handles
	unsigned int tex_id;
	unsigned int verts_vbo, normals_vbo, colors_vbo, texcoords_vbo, prims_ibo, vao;

	// Model matrix
	XForm<float> model;

	// Material info
	material::Phong phong;

private:

	// I could use a base class for meshes, but I don't really want to.
	std::shared_ptr<TriangleMesh> trimeshPtr;
	std::shared_ptr<TetMesh> tetmeshPtr;
	std::unique_ptr<mcl::Texture> texture;

	// Data that is allocated when not found in the wrapped mesh
	std::vector<Vec2f> texcoords_data;
	std::vector<Vec3f> colors_data;

	int last_prim_size;
	inline void init(); // called by constructors
	inline void get_data(); // gets data from the mesh ptr
	inline void subdivide_mesh();
	inline void make_flat();
};


//
//	Implementation
//


RenderMesh::RenderMesh( std::shared_ptr<TriangleMesh> mesh, int opt ){ flags=opt; trimeshPtr=mesh; init(); }
RenderMesh::RenderMesh( std::shared_ptr<TetMesh> mesh, int opt ){ flags=opt; tetmeshPtr=mesh; init(); }

inline void RenderMesh::init(){
	vertices = nullptr;
	num_vertices = 0;
	normals = nullptr;
	num_normals = 0;
	colors = nullptr;
	num_colors = 0;
	texcoords = nullptr;
	num_texcoords = 0;
	prims = nullptr;
	num_prims = 0;
	tex_id = 0;
	verts_vbo = 0;
	normals_vbo = 0;
	colors_vbo = 0;
	texcoords_vbo = 0;
	prims_ibo = 0;
	vao = 0;
	last_prim_size = 0;
	model.setIdentity();
//	if( flags & WIREFRAME ){
//		phong.diff.setZero();
//		phong.amb = mcl::Vec3f(0.5,0,0);
//		phong.spec.setZero();
//	}
//	else { phong = material::autoPhong(); }
	phong = material::autoPhong();
}

inline void RenderMesh::get_data(){

	// Update vertex pointers and double check we have valid data
	if( trimeshPtr ){
		trimeshPtr->need_normals();
		trimeshPtr->get_vertex_data( vertices, num_vertices, normals, num_normals, texcoords, num_texcoords );
		if( flags & WIREFRAME ){
			trimeshPtr->get_primitive_data( 2, prims, num_prims );
			last_prim_size = 2;
		}
		else {
			trimeshPtr->get_primitive_data( 3, prims, num_prims );
			last_prim_size = 3;
		}
	}
	else if( tetmeshPtr ){
		tetmeshPtr->need_normals();
		tetmeshPtr->get_vertex_data( vertices, num_vertices, normals, num_normals, texcoords, num_texcoords );
		if( flags & WIREFRAME ){
			tetmeshPtr->get_primitive_data( 2, prims, num_prims );
			last_prim_size = 2;
		}
		else {
			tetmeshPtr->get_primitive_data( 3, prims, num_prims );
			last_prim_size = 3;
		}
	}

	// Fill colors if none exist
	if( num_colors != num_vertices ){
		if( flags & WIREFRAME ){ colors_data.resize( num_vertices, Vec3f(0,0,0) ); }
		else{ colors_data.resize( num_vertices, Vec3f(1,0,0) ); }
		colors = &colors_data[0][0];
		num_colors = colors_data.size();
	}

	// Fill texcoords if non exist
	if( num_texcoords != num_vertices ){
		texcoords_data.resize( num_vertices, Vec2f(0,0) );
		texcoords = &texcoords_data[0][0];
		num_texcoords = texcoords_data.size();
	}

} // end get data


inline void RenderMesh::load_buffers( int load ){

	get_data();

	// Check to make sure we have data
	if( num_vertices<=0 || num_normals<=0 || num_prims<=0 ){
		std::cerr << "**RenderMesh::update Error: No data for " <<
			(trimeshPtr ? "trimesh " : "tetmesh ") <<
			"(" << num_vertices << ", " << num_normals << ", " <<  num_prims << ")" << std::endl;
		return;
	}

	// Create the buffer for indices
	if( !prims_ibo || (load & (ALL|PRIMS)) ){
		int dim = flags & WIREFRAME ? 2 : 3;
		if( dim != last_prim_size ){
			get_data(); // reload primitive data
		}
		if( !prims_ibo ){
			glGenBuffers(1, &prims_ibo);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, prims_ibo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_prims*sizeof(int)*dim, prims, GL_STATIC_DRAW);
		} else {
			glBindBuffer(GL_ARRAY_BUFFER, prims_ibo);
			glBufferSubData(GL_ARRAY_BUFFER, 0, num_prims*sizeof(int)*dim, prims);
		}
	}

	// Now copy vertex and face data to GPU
	GLenum draw_mode = GL_STATIC_DRAW;
	if( flags & DYNAMIC ){ draw_mode = GL_DYNAMIC_DRAW; }
	const int stride = sizeof(float)*3;

	if( !verts_vbo ){ // Create the buffer for vertices
		glGenBuffers(1, &verts_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, verts_vbo);
		glBufferData(GL_ARRAY_BUFFER, num_vertices*stride, vertices, draw_mode);
	} else if( load & (ALL|VERTICES) ){ // Otherwise update
		glBindBuffer(GL_ARRAY_BUFFER, verts_vbo);
		glBufferSubData( GL_ARRAY_BUFFER, 0, num_vertices*stride, vertices );
	}

	if( !normals_vbo ){ // Create the buffer for normals
		glGenBuffers(1, &normals_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
		glBufferData(GL_ARRAY_BUFFER, num_normals*stride, normals, draw_mode);
	} else if( load & (ALL|NORMALS) ){ // Otherwise update
		glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
		glBufferSubData( GL_ARRAY_BUFFER, 0, num_normals*stride, normals );
	}

	if( !colors_vbo ){ // Create the buffer for colors
		glGenBuffers(1, &colors_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
		glBufferData(GL_ARRAY_BUFFER, num_colors*stride, colors, draw_mode);
	} else if( load & (ALL|COLORS) ){ // Otherwise update
		glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
		glBufferSubData( GL_ARRAY_BUFFER, 0, num_colors*stride, colors );
	}

	 // Create the buffer for tex coords, these won't change
	if( !texcoords_vbo ){
		glGenBuffers(1, &texcoords_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, texcoords_vbo);
		glBufferData(GL_ARRAY_BUFFER, num_texcoords*sizeof(float)*2, texcoords, GL_STATIC_DRAW);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// Create the VAO
	// NOTE: This has to match the shader.
	if( !vao ){

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		// location=0 is the vertex
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, verts_vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, 0);

		// location=1 is the color
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, 0);

		// location=2 is the normal
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, 0);

		// location=3 is the tex coord (eventually)
//		glEnableVertexAttribArray(2);
//		glBindBuffer(GL_ARRAY_BUFFER, texcoords_vbo);
//		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, 0);

		// Done setting data for the vao
		glBindVertexArray(0);
	}

} // end copy data to GPU


inline Eigen::AlignedBox<float,3> RenderMesh::bounds(){

	Eigen::AlignedBox<float,3> aabb;

	if( trimeshPtr ){
		int n_verts = trimeshPtr->vertices.size();
		for( int i=0; i<n_verts; ++i ){
			aabb.extend( trimeshPtr->vertices[i] );
		}
	}
	else if( tetmeshPtr ){
		int n_verts = tetmeshPtr->vertices.size();
		for( int i=0; i<n_verts; ++i ){
			aabb.extend( tetmeshPtr->vertices[i] );
		}
	}

	return aabb;
}

inline void RenderMesh::draw(){
	if( !prims_ibo ){ load_buffers(); }
	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, prims_ibo);
	if( flags & WIREFRAME ){
		glDrawElements(GL_LINES, num_prims*2, GL_UNSIGNED_INT, 0);
	} else {
		glDrawElements(GL_TRIANGLES, num_prims*3, GL_UNSIGNED_INT, 0);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}


inline void RenderMesh::subdivide_mesh(){
throw std::runtime_error( "**TODO: RenderMesh::subdivide_mesh" );
/*
	// Copy vertex data to tempmesh
	tempmesh.vertices.clear(); tempmesh.vertices.reserve( num_vertices );
	tempmesh.faces.clear(); tempmesh.faces.reserve( num_faces );
	tempmesh.texcoords.clear(); tempmesh.texcoords.reserve( num_texcoords );
	for( int i=0; i<num_vertices; ++i ){ tempmesh.vertices.push_back( trimesh::vec( vertices[i*3+0], vertices[i*3+1], vertices[i*3+2] ) ); }
	for( int i=0; i<num_faces; ++i ){ tempmesh.faces.push_back( trimesh::TriMesh::Face( faces[i*3+0], faces[i*3+1], faces[i*3+2] ) ); }
	for( int i=0; i<num_texcoords; ++i ){ tempmesh.texcoords.push_back( trimesh::vec2( texcoords[i*2+0], texcoords[i*2+1] ) ); }

	// Do subdivision with trimesh
	trimesh::subdiv( &tempmesh ); // creates faces
	tempmesh.need_normals(true);

	// We will use the app data stored with that object.
	vertices = &tempmesh.vertices[0][0];
	normals = &tempmesh.normals[0][0];
	faces = &tempmesh.faces[0][0];
	texcoords = &tempmesh.texcoords[0][0];
	num_vertices = tempmesh.vertices.size();
	num_normals = tempmesh.normals.size();
	num_faces = tempmesh.faces.size();
	num_texcoords = tempmesh.texcoords.size();
*/
} // end subdivide mesh


inline void RenderMesh::make_flat(){
throw std::runtime_error( "**TODO: RenderMesh::make_flat" );
/*
	tempmesh.vertices.reserve( num_vertices );
	tempmesh.faces.reserve( num_faces );
	tempmesh.texcoords.reserve( num_texcoords );
	for( int i=0; i<num_faces; ++i ){
		Vec3i f( faces[i*3], faces[i*3+1], faces[i*3+2] );
		int v_idx = tempmesh.vertices.size();
		tempmesh.vertices.push_back( Vec3f( vertices[f[0]*3], vertices[f[0]*3+1], vertices[f[0]*3+2] ) );
		tempmesh.vertices.push_back( Vec3f( vertices[f[1]*3], vertices[f[1]*3+1], vertices[f[1]*3+2] ) );
		tempmesh.vertices.push_back( Vec3f( vertices[f[2]*3], vertices[f[2]*3+1], vertices[f[2]*3+2] ) );
		tempmesh.faces.push_back( Vec3i(v_idx,v_idx+1,v_idx+2) );
		if( num_texcoords ){
			tempmesh.texcoords.push_back( Vec2f( texcoords[f[0]*2], texcoords[f[0]*2+1] ) );
			tempmesh.texcoords.push_back( Vec2f( texcoords[f[1]*2], texcoords[f[1]*2+1] ) );
			tempmesh.texcoords.push_back( Vec2f( texcoords[f[2]*2], texcoords[f[2]*2+1] ) );
		}
	}
	tempmesh.need_normals(true);
	tempmesh.need_edges(true);

	// We will use the app data stored with that object.
	vertices = &tempmesh.vertices[0][0];
	normals = &tempmesh.normals[0][0];
	faces = &tempmesh.faces[0][0];
	texcoords = &tempmesh.texcoords[0][0];
	edges = &tempmesh.edges[0][0];
	num_vertices = tempmesh.vertices.size();
	num_normals = tempmesh.normals.size();
	num_faces = tempmesh.faces.size();
	num_texcoords = tempmesh.texcoords.size();
	num_edges = tempmesh.edges.size();
*/
} // end make flat shading

} // end namespace mcl

#endif
