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

#include "MCL/TetGen.hpp"
#include "MCL/TriangleMesh.hpp"
#include "MCL/TetMesh.hpp"
#include "MCL/ArgParser.hpp"
#include "MCL/MeshIO.hpp"
#include "SliceViewer.hpp"

using namespace mcl;

void help(){
	std::stringstream ss;
	ss << "\n==========================================\nArgs:\n" <<
		"\t-f: .obj file of the closed mesh\n" <<
		"\t-fvol: max tet volume (0,1] as a fraction of the total (recommended)\n" <<
		"\t-vol: max tet volume (optional)\n" <<
		"\t-vis: visualize the tet mesh (optional)\n" <<
		"\t-v: verbose output (optional)\n" <<
	"==========================================\n";
	printf( "%s", ss.str().c_str() );
}

int main(int argc, char *argv[]){

	ArgParser parser(argc,argv);
	if(	parser.exists("-h") || parser.exists("-help") ||
		parser.exists("--h") || parser.exists("--help") ||
		!parser.exists("-f") ){
		help();
		return 0;
	}

	// Parse some args
	std::string filename = parser.get<std::string>("-f");
	bool verbose_output = false;
	parser.get<bool>("-v", &verbose_output);
	float maxvol_percent = -1.f;
	parser.get<float>("-fvol", &maxvol_percent);
	float maxvol = -1.f;
	parser.get<float>("-vol", &maxvol);
	bool vis = parser.exists( "-vis" );

	// Load the triangle mesh
	TriangleMesh trimesh;
	bool success = meshio::load_obj( &trimesh, filename );
	if( !success ){ return 1; }

	// Generate the tet mesh
	std::shared_ptr<mcl::TetMesh> tetmesh = mcl::TetMesh::create();
	tetgen::Settings settings;
	settings.verbose = verbose_output;
	if( maxvol_percent > 0.f ){ settings.maxvol_percent = maxvol_percent; }
	if( maxvol > 0.f ){ settings.maxvol = maxvol; }
	success = tetgen::make_tetmesh( tetmesh->tets, tetmesh->vertices, trimesh.faces, trimesh.vertices, settings );
	if( !success ){ return 1; }

	// Save the tetmesh to file. Should end in .obj otherwise load_obj would have failed.
	std::string fileout = filename.substr(0, filename.size()-4); // remove .obj
	success = meshio::save_elenode( tetmesh.get(), fileout );
	if( !success ){ return 1; }
	std::cout << "Created files:\n\t" << fileout << ".ele\n\t" << fileout << ".node" << std::endl;
	std::cout << "Tetmesh stats:\n\tvertices: " << tetmesh->vertices.size() << "\n\ttets: " << tetmesh->tets.size() << std::endl;

	// Open a render window to visualize the mesh if desired
	if( vis ){
		SliceViewer viewer;
		viewer.set_mesh( tetmesh );
		if( !viewer.display() ){
			return 1;
		}
	}

	return 0;
}

