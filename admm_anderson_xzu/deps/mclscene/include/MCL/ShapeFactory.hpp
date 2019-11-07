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
// Much of the code in this file comes from trimesh2
// by Szymon Rusinkiewicz (gfx.cs.princeton.edu/proj/trimesh2)
//

#ifndef MCL_FACTORY_H
#define MCL_FACTORY_H 1

#include "TriangleMesh.hpp"
#include "TetMesh.hpp"

namespace mcl {
namespace factory {

	// A basic sphere
	static inline std::shared_ptr<TriangleMesh> make_sphere( Vec3f center, float radius, int tess );

	// A 2x2x2 non-symmetric cube centered at origin
	static inline std::shared_ptr<TriangleMesh> make_tri_cube( int tess );

	// Creates connected 1x1x1 cubes, each cube has 5 densly packed tets.
	// Starts at 0,0,0 and builds in the positive axis directions.
	static inline std::shared_ptr<TetMesh> make_tet_blocks( int cubes_x, int cubes_y, int cubes_z );

	// Creates connected 1x1 squares, each cube has 4 triangles (uniform build).
	// CCW winding with normals facing +z
	// Starts at 0,0,0 and builds in the positive x and y directions.
	static inline std::shared_ptr<TriangleMesh> make_tri_blocks( int cubes_x, int cubes_y );

	// A box with specified boxmin/boxmax
	static inline std::shared_ptr<TriangleMesh> make_box( int tess, Vec3f boxmin, Vec3f boxmax );

	// A 2x2 plane along x/y axis
	static inline std::shared_ptr<TriangleMesh> make_plane( int tess_x, int tess_y );

	// Cylinder centered on z axis. tess_c is the triangle resultion at the base.
	static inline std::shared_ptr<TriangleMesh> make_cyl( int tess_c, int tess_l, float r );

	// Create a torus
//	static inline std::shared_ptr<TriangleMesh> make_torus( int tess, float inner_rad, float outer_rad );

	// Helpers from trimesh2
	static inline void mkpoint(std::shared_ptr<TriangleMesh> mesh, float x, float y, float z);
	static inline void mkface(std::shared_ptr<TriangleMesh> mesh, int v1, int v2, int v3);
	static inline void mkquad(std::shared_ptr<TriangleMesh> mesh, int ll, int lr, int ul, int ur);
	static inline void mkquad_sym(std::shared_ptr<TriangleMesh> mesh, int ll, int lr, int ul, int ur, int cent);
	static inline void remove_vertices(TriangleMesh *mesh, const std::vector<bool> &toremove);
	static inline void remove_unused_vertices(TriangleMesh *mesh);
	static inline void remove_faces(TriangleMesh *mesh, const std::vector<bool> &toremove);
	static inline void remap_verts(TriangleMesh *mesh, const std::vector<int> &remap_table);
	static inline void reorder_verts(TriangleMesh *mesh);
	template <typename T> static inline T sqr( T val ){ return val*val; }

}; // end namespace factory

//
//	Implementation
//

static inline std::shared_ptr<TriangleMesh> factory::make_sphere( Vec3f center, float radius, int tess ){

	float m_2pi = M_PI*2.f;
	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	mesh->vertices.reserve(2+tess*(tess-1));
	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j < tess; j++) {
		float th = M_PI * j / tess;
		float z = -cos(th);
		float r = sin(th);
		for (int i = 0; i < tess; i++) {
			float ph = m_2pi * i / tess;
			mkpoint(mesh, r*cos(ph), r*sin(ph), z);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(2*tess*tess - 2*tess);

	for (int i = 0; i < tess; i++)
		mkface(mesh, 0, ((i+1)%tess)+1, i+1);

	for (int j = 0; j < tess-2; j++) {
		int base = 1 + j * tess;
		for (int i = 0; i < tess; i++) {
			int i1 = (i+1)%tess;
			mkquad(mesh, base + i, base + i1,
				base+tess+i, base+tess+i1);
		}
	}

	int base = 1 + (tess-2)*tess;
	for (int i = 0; i < tess; i++){ mkface(mesh, base+i, base+((i+1)%tess), base+tess); }

	// Now scale it by the radius
	XForm<float> s_xf = xform::make_scale(radius,radius,radius);
	mesh->apply_xform( s_xf );

	// Translate so center is correct
	XForm<float> t_xf = xform::make_trans<float>(center[0],center[1],center[2]);
	mesh->apply_xform( t_xf );

	return mesh;

} // end make sphere

/*
static inline std::shared_ptr<TriangleMesh> factory::make_beam( int chunks, int tess ){

	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	for( int b=0; b<chunks; ++b ){

		std::shared_ptr<TriangleMesh> box = make_tri_cube( tess );

		// Now translate box1 and box2
		XForm<float> xf = xform::make_trans<float>( b*2.f, 0, 0 );
		box->apply_xform( xf );
		box->need_normals();

		// Remove faces on the -x and +x
		std::vector<bool> toremove( box->faces.size(), false );
		for( size_t f=0; f<box->faces.size(); ++f ){
			if( b > 0 ){ // remove -x
				if( box->trinorm(f).dot( Vec3f(-1.f,0.f,0.f) ) > 0.f ){ toremove[f]=true; }
			}
			if( b < chunks-1 ){ // remove +x
				if( box->trinorm(f).dot( Vec3f(1.f,0.f,0.f) ) > 0.f ){ toremove[f]=true; }
			}
		}

		remove_faces( box.get(), toremove );
		remove_unused_vertices( box.get() );

		int prev_verts = mesh->vertices.size();

		for( size_t i=0; i<box->vertices.size(); ++i ){ mesh->vertices.push_back( box->vertices[i] ); }
		for( size_t i=0; i<box->faces.size(); ++i ){
			for( int j=0; j<3; ++j ){ box->faces[i][j] += prev_verts; }
			mesh->faces.push_back( box->faces[i] );
		}

	}

	mesh->collapse_points(0.1f);

	return mesh;

} // end make beam
*/

static inline std::shared_ptr<TriangleMesh> factory::make_tri_cube( int tess ){

	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	if (tess < 1)
		tess = 1;

	mesh->vertices.reserve(6*sqr(tess)+2);
	for (int j = 0; j < tess+1; j++) {
		float y = 1.0f - 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, y, -1);
		}
	}
	for (int j = 1; j < tess; j++) {
		float z = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, -1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, 1, y, z);
		}
		for (int i = 0; i < tess; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, 1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, -1, y, z);
		}
	}
	for (int j = 0; j < tess+1; j++) {
		float y = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, y, 1);
		}
	}

	mesh->faces.reserve(12*sqr(tess));
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = i + j * (tess+1);
			mkquad(mesh, ind, ind+tess+1, ind+1, ind+tess+2);
		}
	}

	int topstart = sqr(tess+1) + 4*tess*(tess-1);
	for (int j = 0; j < tess; j++) {
		int next = sqr(tess+1) + 4*tess*(j-1);
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = sqr(tess+1)-1 - i;
				lr = ll - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess*(tess+1) - i*(tess+1);
				lr = ll - (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess + i*(tess+1);
				ur = ul + (tess+1);
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = i;
				lr = i + 1;
			}
			if (j == tess-1) {
				ul = topstart + sqr(tess+1)-1 - i;
				ur = ul - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess + i*(tess+1);
				lr = ll + (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess*(tess+1) - i*(tess+1);
				ur = ul - (tess+1);
			}
			if (i == tess-1) {
				if (j != 0)
					lr -= 4*tess;
				if (j != tess-1)
					ur -= 4*tess;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
	}
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = topstart + i + j * (tess+1);
			mkquad(mesh, ind, ind+1, ind+tess+1, ind+tess+2);
		}
	}

	return mesh;
}

static inline std::shared_ptr<TetMesh> factory::make_tet_blocks( int cubes_x, int cubes_y, int cubes_z ){

	std::shared_ptr<TetMesh> mesh = std::make_shared<TetMesh>();

	cubes_x = std::max(1,cubes_x);
	cubes_y = std::max(1,cubes_y);
	cubes_z = std::max(1,cubes_z);
	Vec3f start(0,0,0);

	for( int x=0; x<cubes_x; ++x ){
	for( int y=0; y<cubes_y; ++y ){
	for( int z=0; z<cubes_z; ++z ){

		Vec3f min = start + Vec3f(x,y,z);
		Vec3f max = min + Vec3f(1,1,1);

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
		int nv = mesh->vertices.size();
		mesh->vertices.emplace_back( a ); // 0
		mesh->vertices.emplace_back( b ); // 1
		mesh->vertices.emplace_back( c ); // 2
		mesh->vertices.emplace_back( d ); // 3
		mesh->vertices.emplace_back( e ); // 4
		mesh->vertices.emplace_back( f ); // 5
		mesh->vertices.emplace_back( g ); // 6
		mesh->vertices.emplace_back( h ); // 7

		// Pack 5 tets into the cube
		Vec4i t0( 0, 5, 7, 4 );
		Vec4i t1( 5, 7, 2, 0 );
		Vec4i t2( 5, 0, 2, 1 );
		Vec4i t3( 7, 2, 0, 3 );
		Vec4i t4( 5, 2, 7, 6 );
		Vec4i offset(nv,nv,nv,nv);

		// Add the tets
		mesh->tets.emplace_back( t0+offset );
		mesh->tets.emplace_back( t1+offset );
		mesh->tets.emplace_back( t2+offset );
		mesh->tets.emplace_back( t3+offset );
		mesh->tets.emplace_back( t4+offset );

	}}} // x, y, z

	// Join colocated vertices and update indices.
	mesh->refine();

	return mesh;

} // end make tetmesh cube

static inline std::shared_ptr<TriangleMesh> factory::make_tri_blocks( int cubes_x, int cubes_y ){

	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	cubes_x = std::max(1,cubes_x);
	cubes_y = std::max(1,cubes_y);
	Vec3f start(0,0,0);
	int z = 0;

	// Adds redundant vertices than removes them with the refine method.
	for( int x=0; x<cubes_x; ++x ){
	for( int y=0; y<cubes_y; ++y ){

		Vec3f min = start + Vec3f(x,y,z);
		Vec3f max = min + Vec3f(1,1,z);
		Vec3f mid = ( min + max ) * 0.5f;


		Vec3f a = min;
		Vec3f b( max[0], min[1], z );
		Vec3f c = max;
		Vec3f d( min[0], max[1], z );
		Vec3f e = mid;

		// Add the verts
		int nv = mesh->vertices.size();
		mesh->vertices.emplace_back( a ); // 0
		mesh->vertices.emplace_back( b ); // 1
		mesh->vertices.emplace_back( c ); // 2
		mesh->vertices.emplace_back( d ); // 3
		mesh->vertices.emplace_back( e ); // 4

		// Set triangel and add to mesh
		Vec3i t0( 3, 0, 4 );
		Vec3i t1( 0, 1, 4 );
		Vec3i t2( 1, 2, 4 );
		Vec3i t3( 2, 3, 4 );
		Vec3i offset(nv,nv,nv);
		mesh->faces.emplace_back( t0+offset );
		mesh->faces.emplace_back( t1+offset );
		mesh->faces.emplace_back( t2+offset );
		mesh->faces.emplace_back( t3+offset );

	}} // x, y

	// Join colocated vertices and update indices.
	mesh->refine();

	return mesh;
}

static inline std::shared_ptr<TriangleMesh> factory::make_box( int tess, Vec3f boxmin, Vec3f boxmax ){

	// Just make a cube and translate/scale it
	std::shared_ptr<TriangleMesh> box = factory::make_tri_cube( tess );
	XForm<float> t1 = xform::make_trans<float>( 1, 1, 1 ); // translate boxmin to origin
	XForm<float> scale1 = xform::make_scale<float>( 0.5, 0.5, 0.5 ); // max box 1x1x1
	XForm<float> scale2 = xform::make_scale<float>( boxmax[0]-boxmin[0], boxmax[1]-boxmin[1], boxmax[2]-boxmin[2] );
	XForm<float> t2 = xform::make_trans<float>( boxmin[0], boxmin[1], boxmin[2] ); // translate to boxmin
	XForm<float> xf = t2 * scale2 * scale1 * t1;
	box->apply_xform( xf );
	return box;

}

static inline std::shared_ptr<TriangleMesh> factory::make_plane( int tess_x, int tess_y ){

	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	if (tess_x < 1)
		tess_x = 1;
	if (tess_y < 1)
		tess_y = 1;

	int n_verts = (tess_x+1)*(tess_y+1)+(tess_x*tess_y);
	mesh->vertices.reserve(n_verts);
//	double y_step = 1.0 / tess_y;
//	double x_step = 1.0 / tess_x;

	// Make grid of vertices
	for( int x=0; x<(tess_x+1); ++x ){
		for( int y=0; y<(tess_y+1); ++y ){
			float xp = -1.0f + 2.0f * x / tess_x;
			float yp = -1.0f + 2.0f * y / tess_y;
			mkpoint(mesh, xp, yp, 0);
		}
	}

	// Make center points
	for( int x=0; x<(tess_x); ++x ){
		for( int y=0; y<(tess_y); ++y ){
			float xp = -1.0f + 2.0f * x / tess_x;
			float yp = -1.0f + 2.0f * y / tess_y;
			xp += 1.f/tess_x;
			yp += 1.f/tess_y;
			mkpoint(mesh, xp, yp, 0);
		}
	}

	assert( n_verts == (int)mesh->vertices.size() );

	// Make faces
	mesh->faces.reserve(tess_x*tess_y*4);
	for( int x=0; x<tess_x; ++x ){
		for( int y=0; y<tess_y; ++y ){
			int ll=y+x*(tess_y+1);
			int lr=y+(x+1)*(tess_y+1);
			int ul=ll+1;
			int ur=lr+1;
			int cent=(tess_x+1)*(tess_y+1) + x*tess_y + y;
			mkquad_sym(mesh, ll, lr, ul, ur, cent);
		}
	}

	// Make texture coordinates
	mesh->texcoords.reserve( mesh->vertices.size() );
	for( size_t i=0; i<mesh->vertices.size(); ++i ){
		Vec3f p = mesh->vertices[i];
		float u = (p[0]+1.f)/2.f;
		float v = 1.f-(p[1]+1.f)/2.f;
		mesh->texcoords.push_back( Vec2f(u,v) );
	}

	return mesh;

}

static inline std::shared_ptr<TriangleMesh> factory::make_cyl(int tess_c, int tess_l, float r){

	float m_2pi = M_PI*2.f;
	std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();

	if (tess_c < 3)
		tess_c = 3;
	if (tess_l < 1)
		tess_l = 1;

	mesh->vertices.reserve(2+3*tess_c*tess_l-tess_c);

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j <= tess_l; j++) {
		float rr = r * j / tess_l;
		for (int i = 0; i < tess_c; i++) {
			float th = m_2pi * i / tess_c;
			mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
		}
	}
	int side_start = mesh->vertices.size();
	for (int j = 1; j < tess_l; j++) {
		float z = -1.0f + 2.0f * j / tess_l;
		for (int i = 0; i < tess_c; i++) {
			float th = m_2pi * i / tess_c;
			mkpoint(mesh, r*cos(th), r*sin(th), z);
		}
	}
	int top_start = mesh->vertices.size();
	for (int j = tess_l; j > 0; j--) {
		float rr = r * j / tess_l;
		for (int i = 0; i < tess_c; i++) {
			float th = m_2pi * i / tess_c;
			mkpoint(mesh, rr*cos(th), rr*sin(th), 1);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(6*tess_c*tess_l - 2*tess_c);

	for (int i = 0; i < tess_c; i++)
		mkface(mesh, 0, ((i+1)%tess_c)+1, i+1);
	for (int j = 1; j < tess_l; j++) {
		int base = 1 + (j-1) * tess_c;
		for (int i = 0; i < tess_c; i++) {
			int i1 = (i+1)%tess_c;
			mkquad(mesh, base+tess_c+i1, base+tess_c+i,
				base+i1, base+i);
		}
	}

	for (int j = 0; j < tess_l; j++) {
		int base = side_start - tess_c + j * tess_c;
		for (int i = 0; i < tess_c; i++) {
			int i1 = (i+1)%tess_c;
			mkquad(mesh, base + i, base + i1,
				base+tess_c+i, base+tess_c+i1);
		}
	}

	for (int j = 0; j < tess_l-1; j++) {
		int base = top_start + j * tess_c;
		for (int i = 0; i < tess_c; i++) {
			int i1 = (i+1)%tess_c;
			mkquad(mesh, base+tess_c+i1, base+tess_c+i,
				base+i1, base+i);
		}
	}
	int base = top_start + (tess_l-1)*tess_c;
	for (int i = 0; i < tess_c; i++)
		mkface(mesh, base+i, base+((i+1)%tess_c), base+tess_c);
	return mesh;
}

/*
static inline std::shared_ptr<TriangleMesh> factory::make_torus( int tess, float inner_rad, float outer_rad ){

	if (tess < 3) tess = 3;

	make_ccyl( mesh, tess, tess, outer_rad );

	for (int i = 0; i < tess; i++)
		mesh->vertices.pop_back();
	for (size_t i = 0; i < mesh->faces.size(); i++) {
		mesh->faces[i][0] %= mesh->vertices.size();
		mesh->faces[i][1] %= mesh->vertices.size();
		mesh->faces[i][2] %= mesh->vertices.size();
	}

	float r = inner_rad;

	for (int j = 0; j < tess; j++) {
		float th = m_2pi * j / tess;
		Vec3f circlepos(cos(th), sin(th), 0);
		for (int i = 0; i < tess; i++) {
			float ph = m_2pi * i / tess;
			mesh->vertices[i+j*tess] = circlepos +
				cosf(ph)*r*circlepos +
				sinf(ph)*r*vec(0,0,-1);
		}
	}


*/

//
//	Helper functions
//

static inline void factory::mkpoint(std::shared_ptr<TriangleMesh> mesh, float x, float y, float z){
	mesh->vertices.push_back( Vec3f(x,y,z) );
}

static inline void factory::mkface(std::shared_ptr<TriangleMesh> mesh, int v1, int v2, int v3){
	mesh->faces.push_back( Vec3i(v1, v2, v3) );
}

// non symmetric quad
//	*---*
//	|  /|
//	| / |
//	|/  |
//	*---*
static inline void factory::mkquad(std::shared_ptr<TriangleMesh> mesh, int ll, int lr, int ul, int ur){
	mkface(mesh, ll, lr, ur);
	mkface(mesh, ll, ur, ul);
}

// Symmetric quad
//	*---*
//	|\ /|
//	| * |
//	|/ \|
//	*---*
static inline void factory::mkquad_sym(std::shared_ptr<TriangleMesh> mesh, int ll, int lr, int ul, int ur, int cent){
	// Counter clockwise
	mkface(mesh, ll, lr, cent);
	mkface(mesh, lr, ur, cent);
	mkface(mesh, cent, ur, ul);
	mkface(mesh, ll, cent, ul);
}


// Remove the indicated vertices from the TriMesh.
static inline void factory::remove_vertices(TriangleMesh *mesh, const std::vector<bool> &toremove){

	int nv = mesh->vertices.size();

	// Build a table that tells how the vertices will be remapped
	if (!nv)
		return;

	std::vector<int> remap_table(nv);
	int next = 0;
	for (int i = 0; i < nv; i++) {
		if (toremove[i])
			remap_table[i] = -1;
		else
			remap_table[i] = next++;
	}

	// Nothing to delete?
	if (next == nv) {
		return;
	}

	remap_verts(mesh, remap_table);

}

static inline void factory::remove_unused_vertices(TriangleMesh *mesh){

	int nv = mesh->vertices.size();
	if (!nv) return;

	int nf = mesh->faces.size();
	std::vector<bool> unused(nv, true);
	for (int i = 0; i < nf; i++) {
		unused[mesh->faces[i][0]] = false;
		unused[mesh->faces[i][1]] = false;
		unused[mesh->faces[i][2]] = false;
	}

	remove_vertices(mesh, unused);
}

static inline void factory::remove_faces(TriangleMesh *mesh, const std::vector<bool> &toremove){

	int numfaces = mesh->faces.size();
	if (!numfaces) return;

	int next = 0;
	for (int i = 0; i < numfaces; i++) {
		if (toremove[i])
			continue;
		mesh->faces[next++] = mesh->faces[i];
	}
	if (next == numfaces) {
		return;
	}

	mesh->faces.erase(mesh->faces.begin() + next, mesh->faces.end());
}


// Remap vertices according to the given table
//
// Faces are renumbered to reflect the new numbering of vertices, and any
// faces that included a vertex that went away will also be removed.
//
// Any per-vertex properties are renumbered along with the vertices.
static inline void factory::remap_verts(TriangleMesh *mesh, const std::vector<int> &remap_table){

	if (remap_table.size() != mesh->vertices.size()) {
		printf("remap_verts called with wrong table size!\n");
		return;
	}

	// Check what we're doing
//	bool removing_verts = false;
	bool any_left = false;
	int last = -1;
	int nv = mesh->vertices.size();
	for (int i = 0; i < nv; i++) {
		if (remap_table[i] < 0) {
//			removing_verts = true;
		} else {
			any_left = true;
			if (remap_table[i] > last)
				last = remap_table[i];
		}
	}

	if (!any_left) {
		mesh->clear();
		return;
	}

	// Figure out what we have sitting around, so we can remap/recompute
//	bool have_faces = !mesh->faces.empty();
	bool have_normals = !mesh->normals.empty();

	// Remap the vertices and per-vertex properties
	TriangleMesh *oldmesh = new TriangleMesh;
	*oldmesh = *mesh;

#define REMAP(property) mesh->property[remap_table[i]] = oldmesh->property[i]

	for (int i = 0; i < nv; i++) {
		if (remap_table[i] < 0 || remap_table[i] == i)
			continue;
		REMAP(vertices);
		if (have_normals) REMAP(normals);
	}

#define ERASE(property) mesh->property.erase(mesh->property.begin()+last+1, \
					     mesh->property.end())
	ERASE(vertices);
	if (have_normals) ERASE(normals);

	// Renumber faces
	int nf = mesh->faces.size(), nextface = 0;
	for (int i = 0; i < nf; i++) {
		int n0 = (mesh->faces[nextface][0] = remap_table[oldmesh->faces[i][0]]);
		int n1 = (mesh->faces[nextface][1] = remap_table[oldmesh->faces[i][1]]);
		int n2 = (mesh->faces[nextface][2] = remap_table[oldmesh->faces[i][2]]);
		if ((n0 >= 0) && (n1 >= 0) && (n2 >= 0))
			nextface++;
	}
	mesh->faces.erase(mesh->faces.begin() + nextface, mesh->faces.end());

	delete oldmesh;

}


// Reorder vertices in a mesh according to the order in which
// they are referenced by the grid, tstrips, or faces.
static inline void factory::reorder_verts(TriangleMesh *mesh){

	int nv = mesh->vertices.size();
	std::vector<int> remap(nv, -1);
	int next = 0;

	for (size_t i = 0; i < mesh->faces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int v = mesh->faces[i][j];
			if (remap[v] == -1)
				remap[v] = next++;
		}
	}

	if (next != nv) {
		// Unreferenced vertices...  Just stick them at the end.
		for (int i = 0; i < nv; i++)
			if (remap[i] == -1)
				remap[i] = next++;
	}

	remap_verts(mesh, remap);
}


} // end namespace mcl

#endif

