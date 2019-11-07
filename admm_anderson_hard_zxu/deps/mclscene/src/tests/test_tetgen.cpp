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
#include "MCL/TetMesh.hpp"
#include "MCL/MeshIO.hpp"
#include "MCL/TetGen.hpp"

using namespace mcl;

bool test_bunny();

int main(void){
	if( !test_bunny() ){ return EXIT_FAILURE; }
	return EXIT_SUCCESS;
}

bool test_bunny(){

	std::cout << "\n\nTesting closed bunny\n\n" << std::endl;
	{
		mcl::TriangleMesh bunny;
		std::stringstream bunnyfile;
		bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny_closed.obj";
		mcl::meshio::load_obj( &bunny, bunnyfile.str() );

		std::cout << "Tri Bunny has " << bunny.vertices.size() << " verts" << std::endl;

		mcl::TetMesh tetbunny;
		tetgen::Settings settings;
		settings.verbose = true;
		settings.maxvol_percent = 0.1;
		bool s = tetgen::make_tetmesh( tetbunny.tets, tetbunny.vertices, bunny.faces, bunny.vertices, settings );
		if( !s ){ return false; }

		std::cout << "Tet bunny has " << tetbunny.vertices.size() << " verts, and " <<
			tetbunny.tets.size() << " tets" << std::endl;
	}

	// Test on an open bunny
	{
		mcl::TriangleMesh bunny;
		std::stringstream bunnyfile;
		bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny.obj";
		mcl::meshio::load_obj( &bunny, bunnyfile.str() );
		bool s = tetgen::verify_closed( bunny.faces, bunny.vertices );
		if( s ){
			std::cout << "Failed to recognize open mesh" << std::endl;
			return false;
		}
	}

	return true;
}

