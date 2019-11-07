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

#include "MCL/TetMesh.hpp"
#include "MCL/ArgParser.hpp"
#include "MCL/MeshIO.hpp"
#include "SliceViewer.hpp"

using namespace mcl;

void help(){
	std::stringstream ss;
	ss << "\n==========================================\nArgs:\n" <<
		"\t-f: filename file of the tetmesh (without extension)\n" <<
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

	std::string filename = parser.get<std::string>("-f");

	// Generate the tet mesh
	std::shared_ptr<mcl::TetMesh> tetmesh = mcl::TetMesh::create();
	bool success = meshio::load_elenode( tetmesh.get(), filename );
	if( !success ){ return 1; }

	// Open a render window to visualize the mesh
	SliceViewer viewer;
	viewer.set_mesh( tetmesh );
	if( !viewer.display() ){
		return 1;
	}

	return 0;
}

