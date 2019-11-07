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

#ifndef MCL_MESHIO_H
#define MCL_MESHIO_H

#include "TriangleMesh.hpp"
#include "TetMesh.hpp"
#include <iostream>
#include <fstream>

namespace mcl {
namespace meshio {

	static inline bool load_obj( TriangleMesh *mesh, std::string file, bool normals=false, bool texcoords=false, bool colors=false );

	static inline bool save_obj( const TriangleMesh *mesh, std::string file );

	static inline bool load_elenode( TetMesh *mesh, std::string file );

	static inline bool save_elenode( const TetMesh *mesh, std::string file );

}; // namespace io


//
// Implementation
//

namespace io_helper {
	static std::string to_lower( std::string s ){
		std::transform( s.begin(), s.end(), s.begin(), ::tolower );
		return s;
	}
}

static inline bool meshio::load_obj( TriangleMesh *mesh, std::string file,
	bool normals, bool texcoords, bool colors ){

	mesh->clear();
	std::ifstream infile( file.c_str() );

	if( infile.is_open() ){

		std::string line;
		while( std::getline( infile, line ) ){

			std::stringstream ss(line);
			std::string tok; ss >> tok;
			tok = io_helper::to_lower(tok);

			if( tok == "v" ){ // Vertex
				float x, y, z; ss >> x >> y >> z; // vertices
				mesh->vertices.emplace_back( Vec3f(x,y,z) );
				if( colors ){
//					float cx, cy, cz; // colors
//					if( ss >> cx >> cy >> cz ){ mesh->colors.emplace_back( Vec3f(cx,cy,cz) ); }
				}
			}

			else if( tok == "vt" && texcoords ){ // Tex coord
				float u, v; ss >> u >> v;
				mesh->texcoords.emplace_back( Vec2f(u,v) );
			}

			else if( tok == "vn" && normals ){ // Normal
				float x, y, z; ss >> x >> y >> z; // vertices
				mesh->normals.emplace_back( Vec3f(x,y,z) );
			}

			else if( tok == "f" ){ // face
				Vec3i face;
				for( size_t i=0; i<3; ++i ){ // Get the three vertices
					std::vector<int> f_vals;
					{ // Split the string with the / delim
						std::string f_str; ss >> f_str;
						std::stringstream ss2(f_str); std::string s2;
						while( std::getline(ss2, s2, '/') ){
							if( s2.size()==0 ){ continue; }
							f_vals.emplace_back( std::stoi(s2)-1 );
						}
					}
					face[i] = f_vals.size() > 0 ? f_vals[0] : -1;
				}
				if( face[0]>=0 && face[1]>=0 && face[2]>=0 ){ mesh->faces.emplace_back(face); }
			} // end parse face

		} // end loop lines

	} // end load obj
	else { std::cerr << "\n**mcl::meshio::load_obj Error: Could not open file " << file << std::endl; return false; }

	// Double check our file
	if( mesh->texcoords.size() != mesh->vertices.size() && mesh->texcoords.size()>0 && texcoords ){
		std::cerr << "\n**mcl::meshio::load_obj Error: Failed to load texture coordinates." << std::endl;
		mesh->texcoords.clear();
		return false;
	}

	// Double check our file
	if( mesh->normals.size() != mesh->vertices.size() && mesh->normals.size()>0 && normals ){
		std::cerr << "\n**mcl::meshio::load_obj Warning: Normals should be per-vertex (removing them)." << std::endl;
		mesh->normals.clear();
		return false;
	}

	return true;

} // end load obj


static inline bool meshio::save_obj( const TriangleMesh *mesh, std::string filename ){

	bool suppress_tex = true;
	bool suppress_normals = false;

	int fsize = filename.size();
	if( fsize<4 ){ printf("\n**TriangleMesh::save Error: Filetype must be .obj\n"); return false; }
	std::string ftype = filename.substr( fsize-4,4 );
	if( ftype != ".obj" ){ printf("\n**TriangleMesh::save Error: Filetype must be .obj\n"); return false; }

	std::ofstream fs;
	fs.open( filename.c_str() );
	fs << "# Generated with mclscene by Matt Overby (www.mattoverby.net)";

	int nv = mesh->vertices.size();
	int nn = mesh->normals.size();
	int nt = mesh->texcoords.size();
	if( suppress_tex ){ nt=0; }
	if( suppress_normals ){ nn=0; }

	for( int i=0; i<nv; ++i ){
		fs << "\nv " << mesh->vertices[i][0] << ' ' << mesh->vertices[i][1] << ' ' << mesh->vertices[i][2];
	}
	for( int i=0; i<nn; ++i ){
		fs << "\nvn " << mesh->normals[i][0] << ' ' << mesh->normals[i][1] << ' ' << mesh->normals[i][2];
	}
	for( int i=0; i<nt; ++i ){
		fs << "\nvt " << mesh->texcoords[i][0] << ' ' << mesh->texcoords[i][1];
	}
	int nf = mesh->faces.size();
	for( int i=0; i<nf; ++i ){
		Vec3i f = mesh->faces[i]; f[0]+=1; f[1]+=1; f[2]+=1;
		if( nn==0 && nt==0 ){ // no normals or texcoords
			fs << "\nf " << f[0] << ' ' << f[1] << ' ' << f[2];
		} else if( nn==0 && nt>0 ){ // no normals, has texcoords
			fs << "\nf " << f[0]<<'/'<<f[0]<< ' ' << f[1]<<'/'<<f[1]<< ' ' << f[2]<<'/'<<f[2];
		} else if( nn>0 && nt==0 ){ // has normals, no texcoords
			fs << "\nf " << f[0]<<"//"<<f[0]<< ' ' << f[1]<<"//"<<f[1]<< ' ' << f[2]<<"//"<<f[2];
		} else if( nn>0 && nt>0 ){ // has normals and texcoords
			fs << "\nf " << f[0]<<'/'<<f[0]<<'/'<<f[0]<< ' ' << f[1]<<'/'<<f[1]<<'/'<<f[1]<< ' ' << f[2]<<'/'<<f[2]<<'/'<<f[2];
		}
	} // end loop faces

	fs << "\n";
	fs.close();
	return true;

} // end save obj


static inline bool meshio::load_elenode( TetMesh *mesh, std::string file ){

	mesh->clear();

	{ // Load ele

		// Load the vertices of the tetmesh
		std::stringstream ele_file; ele_file << file << ".ele";
		std::ifstream filestream;
		filestream.open( ele_file.str().c_str() );
		if( !filestream ){
			std::cerr << "\n**TetMesh Error: Could not load " << ele_file.str() << std::endl;
			return false;
		}

		std::string header;
		std::getline( filestream, header );
		std::stringstream headerSS(header);
		int n_tets = 0; headerSS >> n_tets;

		mesh->tets.resize( n_tets );
		std::vector< int > tet_set( n_tets, 0 );
		bool starts_with_one = false;

		for( int i=0; i<n_tets; ++i ){
			std::string line;
			std::getline( filestream, line );

			std::stringstream lineSS(line);
			size_t idx;
			int node_ids[4];
			lineSS >> idx >> node_ids[0] >> node_ids[1] >> node_ids[2] >> node_ids[3];

			// Check for 1-indexed
			if( i==0 && idx==1 ){ starts_with_one = true; }
			if( starts_with_one ){
				idx -= 1;
				for( int j=0; j<4; ++j ){ node_ids[j]-=1; }
			}

			if( idx > mesh->tets.size() ){
				std::cerr << "\n**TetMesh Error: Your indices are bad for file " << ele_file.str() << std::endl; return false;
			}

			mesh->tets[idx] = Vec4i( node_ids[0], node_ids[1], node_ids[2], node_ids[3] );
			tet_set[idx] = 1;
		}
		filestream.close();

		size_t tet_set_size = tet_set.size();
		for( size_t i=0; i<tet_set_size; ++i ){
			if( tet_set[i] == 0 ){
				std::cerr << "\n**TetMesh Error: Your indices are bad for file " << ele_file.str() << std::endl;
				return false;
			}
		}
	}

	{ // Load node

		// Load the vertices of the tetmesh
		std::stringstream node_file; node_file << file << ".node";
		std::ifstream filestream;
		filestream.open( node_file.str().c_str() );
		if( !filestream ){
			std::cerr << "\n**TetMesh Error: Could not load " << node_file.str() << std::endl;
			return false;
		}

		std::string header;
		getline( filestream, header );
		std::stringstream headerSS(header);
		int n_nodes = 0; headerSS >> n_nodes;

		mesh->vertices.resize( n_nodes );
		std::vector< int > vertex_set( n_nodes, 0 );
		bool starts_with_one = false;

		for( int i=0; i<n_nodes; ++i ){
			std::string line;
			std::getline( filestream, line );

			std::stringstream lineSS(line);
			double x, y, z;
			size_t idx;
			lineSS >> idx >> x >> y >> z;

			// Check for 1-indexed
			if( i==0 && idx==1 ){ starts_with_one = true; }
			if( starts_with_one ){ idx -= 1; }

			if( idx > mesh->vertices.size() ){
				std::cerr << "\n**TetMesh Error: Your indices are bad for file " << node_file.str() << std::endl;
				return false;
			}

			mesh->vertices[idx] = Vec3f( x, y, z );
			vertex_set[idx] = 1;
		}
		filestream.close();

		size_t vert_set = vertex_set.size();
		for( size_t i=0; i<vert_set; ++i ){
			if( vertex_set[i] == 0 ){
				std::cerr << "\n**TetMesh Error: Your indices are bad for file " << node_file.str() << std::endl;
				return false;
			}
		}

	}

	// Check for inverted tets, and reorder if needed
	int n_tets = mesh->tets.size();
	for( int i=0; i<n_tets; ++i ){
		Vec4i tet = mesh->tets[i];
		Vec3f a = mesh->vertices[tet[0]];
		float V = (mesh->vertices[tet[1]]-a).dot(
			(mesh->vertices[tet[2]]-a).cross(mesh->vertices[tet[3]]-a)
		) / 6.f;
		if( V < 0 ){ // reorder base
			mesh->tets[i][1] = tet[2];
			mesh->tets[i][2] = tet[1];
		}
	}

	if( mesh->vertices.size() == 0 || n_tets == 0 ){
		throw std::runtime_error("\n**TetMesh Error: Problem loading files" );
	}

	return true;

} // end load elenode


static inline bool meshio::save_elenode( const TetMesh *mesh, std::string file ){

	//	ele file
	// http://wias-berlin.de/software/tetgen/fformats.ele.html
	// First line: <# of tetrahedra> <nodes per tetrahedron> <# of attributes>
	// Remaining lines list of # of tetrahedra:
	// <tetrahedron #> <node> <node> <node> <node> ... [attributes]
	{
		std::stringstream elefn;
		elefn << file << ".ele";
		std::ofstream filestream;
		filestream.open( elefn.str().c_str() );
		int n_tets = mesh->tets.size();
		filestream << n_tets << " 4 0\n";
		for( int i=0; i<n_tets; ++i ){
			const Vec4i &tet = mesh->tets[i];
			filestream << "\t" << i << ' ' << tet[0] << ' ' << tet[1] << ' ' << tet[2] << ' ' << tet[3] << "\n";
		}
		filestream << "# Generated by mclscene (www.mattoverby.net)";
		filestream.close();
	} // end ele file

	//	node file
	// http://wias-berlin.de/software/tetgen/fformats.node.html
	// First line: <# of points> <dimension (must be 3)> <# of attributes> <# of boundary markers (0 or 1)>
	// Remaining lines list # of points:
	// <point #> <x> <y> <z> [attributes] [boundary marker]
	{
		std::stringstream nodefn;
		nodefn << file << ".node";
		std::ofstream filestream;
		filestream.open( nodefn.str().c_str() );
		int n_verts = mesh->vertices.size();
		filestream << n_verts << " 3 0 0\n";
		for( int i=0; i<n_verts; ++i ){
			const Vec3f &vert = mesh->vertices[i];
			filestream << "\t" << i << ' ' << vert[0] << ' ' << vert[1] << ' ' << vert[2] << "\n";
		}
		filestream << "# Generated by mclscene (www.mattoverby.net)";
		filestream.close();		
	} // end node file

	return true;
}

}; // namespace mcl

#endif
