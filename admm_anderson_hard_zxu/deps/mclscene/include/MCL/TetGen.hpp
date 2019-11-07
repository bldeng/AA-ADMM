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
//	NOTE: Requires Tetgen 1.5.0
//	This file is just a wrapper for Tetgen (wias-berlin.de/software/tetgen)
//	by Hang Si. It's assumed you've already linked the appropriate libraries.
//
//	A copy of tetgen version 1.5.0 is in the src/tetgen directory. For an
//	example on how to build and link it, see the CMakeLists.txt in mclscene
//	root directory and the src/tests/test_tetgen test file.
//

#ifndef MCL_TETGEN_H
#define MCL_TETGEN_H

#include "Vec.hpp"
#include "HashKeys.hpp"
#include <iostream>

#define TETLIBRARY
#include "tetgen.h"

namespace mcl {

namespace tetgen {

	struct Settings {
		bool verbose;
		float quality;
		float maxvol;
		float maxvol_percent; // compute max vol as a fraction (0-1) of the total
		Settings() : verbose(false), quality(-1), maxvol(-1), maxvol_percent(-1) {}
		void print() const;
	};

	// Tetrahedralizes a triangle mesh, returns true on success.
	static inline bool make_tetmesh( std::vector<Vec4i> &tets, std::vector<Vec3f> &tet_verts,
		const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts,
		const Settings &settings = Settings() );

	// Makes a tetgenio object from a triangle mesh
	static inline void make_tetgenio( tetgenio &tgio, const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts );

	// Verifies a triangle mesh is closed, returns false if not.
	// Note that a mesh must be "clean" for this to work, i.e. all vertices are a part of a triangle.
	static inline bool verify_closed( const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts );

} // end ns tetgen

//
//  Implementation
//

static bool tetgen::make_tetmesh( std::vector<Vec4i> &tets, std::vector<Vec3f> &tet_verts,
		const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts, const Settings &settings ){

	if( !verify_closed( tris, tri_verts ) ){
		std::cerr << "**TetGen Error: Triangle mesh is not closed" << std::endl;
		return false;
	}

	// Make the tetgen data type
	tetgenio in;
	make_tetgenio( in, tris, tri_verts );

	// Compute maxvol if we want it
	float maxvol = settings.maxvol;
	if( settings.maxvol_percent > 0.f ){

		// Note a totally correct volume, but I'm too lazy to be precise.
		Eigen::AlignedBox<float,3> aabb;
		int n_triverts = tri_verts.size();
		for( int i=0; i<n_triverts; ++i ){
			aabb.extend( tri_verts[i] );
		}
		float mvp = settings.maxvol_percent;
		if( mvp <= 0.f || mvp > 1.f ){
			std::cerr << "**TetGen Error: maxvol_percent should be between 0 and 1" << std::endl;
			return false;
		}
		maxvol = mvp * aabb.volume();
	}

	// Set up the switches
	std::stringstream switches;
	if( settings.verbose ){ switches << "V"; }
	else{ switches << "Q"; }
	if( settings.quality > 0 ){ switches << "q" << settings.quality; }
	if( maxvol > 0 ){ switches << "a" << maxvol; }
	std::cout << "Tetgen Switches: " << switches.str() << std::endl;

	char *c_switches = new char[switches.str().length()+1];
	std::strcpy(c_switches,switches.str().c_str());

	tetgenio out;
	tetrahedralize(c_switches, &in, &out);
	delete c_switches; // don't need switches anymore

	// Make sure we had success
	if( out.numberoftetrahedra == 0 || out.numberofpoints == 0 ){
		std::cerr << "**TetGen Error: Failed to tetrahedralize with settings:" << std::flush;
		settings.print();
		return false;
	}

	tets.clear(); // Copy tets from tetgenio
	tets.reserve( out.numberoftetrahedra );
	for( int i=0; i<out.numberoftetrahedra; ++i ){
		tets.emplace_back(
			out.tetrahedronlist[i*4+0],
			out.tetrahedronlist[i*4+1],
			out.tetrahedronlist[i*4+2],
			out.tetrahedronlist[i*4+3]
		);
	}

	tet_verts.clear(); // Copy verts
	tet_verts.reserve( out.numberofpoints );
	for( int i=0; i<out.numberofpoints; ++i ){
		tet_verts.emplace_back(
			out.pointlist[i*3+0],
			out.pointlist[i*3+1],
			out.pointlist[i*3+2]
		);
	}

	return true;

} // end tetrahedralize

static void tetgen::make_tetgenio( tetgenio &tgio, const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts ){

	// We create copies of vertices/triangle faces since TetGen likes
	// to deallocate on its own (tetgenio destructor). We could set the
	// ptrs to null so it doesn't attempt to deallocate, but I'm not sure
	// if TetGen likes to do other things to the vertices/faces and if
	// that would cause problems. I've tested it and it seems to work fine,
	// but for the sake of robustness I'll just make copies.

	tgio.firstnumber = 0;
	tgio.mesh_dim = 3;
	tgio.numberofpoints = tri_verts.size();
	tgio.pointlist = new double[tgio.numberofpoints * 3];
	for( int i=0; i < tgio.numberofpoints; ++i ){ // Copy verts
		tgio.pointlist[i*3+0] = tri_verts[i][0];
		tgio.pointlist[i*3+1] = tri_verts[i][1];
		tgio.pointlist[i*3+2] = tri_verts[i][2];
	}

//	Doesn't work with max volume setting?
//	tgio.numberoftrifaces = tris.size();
//	tgio.trifacelist = new int[tgio.numberoftrifaces * 3];
//	for( int i=0; i < tgio.numberoftrifaces; ++i ){ // Copy verts
//		mcl::Vec3i f = tris[i];
//		tgio.trifacelist[i*3+0] = f[0];
//		tgio.trifacelist[i*3+1] = f[1];
//		tgio.trifacelist[i*3+2] = f[2];
//	}

	tgio.numberoffacets = tris.size();
	tgio.facetlist = new tetgenio::facet[tgio.numberoffacets];
	tgio.facetmarkerlist = new int[tgio.numberoffacets];

	int n_faces = tris.size();
	for( int i=0; i < n_faces; ++i ){ // Copy tris
		tgio.facetmarkerlist[i] = i;
		tetgenio::facet &f = tgio.facetlist[i];

		// Assuming we don't have holes...
		f.numberofholes = 0;
		f.holelist = NULL;

		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[1];
		tetgenio::polygon &p = f.polygonlist[0];

		p.numberofvertices = 3;
		p.vertexlist = new int[3];
		p.vertexlist[0] = tris[i][0];
		p.vertexlist[1] = tris[i][1];
		p.vertexlist[2] = tris[i][2];
	}
}


static bool tetgen::verify_closed( const std::vector<Vec3i> &tris, const std::vector<Vec3f> &tri_verts ){

	// We need to compute the number of UNIQUE edges.
	// vertex ids -> number of faces using these indices
	std::unordered_map< hashkey::sint2, int > edge_ids;
	int n_faces = tris.size();
	for( int f=0; f<n_faces; ++f ){
		edge_ids.emplace( std::make_pair( hashkey::sint2(tris[f][0],tris[f][1]), 1) );
		edge_ids.emplace( std::make_pair( hashkey::sint2(tris[f][0],tris[f][2]), 1) );
		edge_ids.emplace( std::make_pair( hashkey::sint2(tris[f][1],tris[f][2]), 1) );
	}
	int n_edges = edge_ids.size();
	int n_verts = tri_verts.size();
	if( n_verts + n_faces - n_edges != 2 ){ return false; }
	return true;
}

void tetgen::Settings::print() const {
	std::cout <<
		"\n\t verbose: " << verbose <<
		"\n\t quality: " << quality <<
		"\n\t maxvol: " << maxvol << ( maxvol <= 0 ? " (no constraint)" : " " ) <<
		"\n\t maxvol_percent: " << maxvol_percent << ( maxvol_percent <= 0 ? " (no constraint)" : " " ) <<
	std::endl;
} // end print settings


} // ns mcl

#endif
