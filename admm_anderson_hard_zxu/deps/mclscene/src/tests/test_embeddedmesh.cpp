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

#include <iostream>
#include "MCL/EmbeddedMesh.hpp"
#include "MCL/MeshIO.hpp"

using namespace mcl;

static float volume( const std::vector<Vec3f> &verts, const mcl::Vec4i &tet ){
	Eigen::Matrix<float,3,3> edges;
	edges.col(0) = verts[tet[1]] - verts[tet[0]];
	edges.col(1) = verts[tet[2]] - verts[tet[0]];
	edges.col(2) = verts[tet[3]] - verts[tet[0]];
	return (edges).determinant()/6.f;
}

int main(void){

	EmbeddedMesh mesh;

	// Load the mesh
	std::stringstream bunnyfile;
	bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny.obj";
	mcl::meshio::load_obj( mesh.embedded.get(), bunnyfile.str() );
	XForm<float> xf = xform::make_scale<float>(10,10,10);
	mesh.embedded->apply_xform(xf);

	// Generate the lattice
	if( !mesh.gen_lattice() ){
		std::cerr << "Failed to generate lattice" << std::endl;
		return EXIT_FAILURE;
	}


	int n_tets = mesh.lattice->tets.size();
	int n_verts = mesh.lattice->vertices.size();
	std::vector<int> vert_refs( n_verts, 0 );
	for( int i=0; i<n_tets; ++i ){
		Vec4i tet = mesh.lattice->tets[i];

		// Make sure no tets are inverted
		float v = volume( mesh.lattice->vertices, tet );
		if( v <= 0 ){
			std::cerr << "Inverted tets:\n\t volume: " << v <<
				"\n\t tet (" << i << "): " << tet.transpose() <<
				"\n\t v0: " << mesh.lattice->vertices[tet[0]].transpose() <<
				"\n\t v1: " << mesh.lattice->vertices[tet[1]].transpose() <<
				"\n\t v2: " << mesh.lattice->vertices[tet[2]].transpose() <<
				"\n\t v3: " << mesh.lattice->vertices[tet[3]].transpose() <<
			std::endl;
			return EXIT_FAILURE;
		}

		for( int j=0; j<4; ++j ){ vert_refs[ tet[j] ] += 1; }

	} // end loop tets

	// Compute masses
	std::vector<float> masses;
	mesh.weighted_masses( masses, 1000.f );

	for( int i=0; i<n_verts; ++i ){
		if( vert_refs[i] == 0 ){
			std::cerr << "Had unreferenced vertices!" << std::endl;
			return EXIT_FAILURE;
		}
		if( masses[i] <= 0.f ){
			std::cerr << "Zero mass for idx " << i << std::endl;
			return EXIT_FAILURE;
		}
	}
		

	return EXIT_SUCCESS;
}

