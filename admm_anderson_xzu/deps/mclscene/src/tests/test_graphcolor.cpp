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
#include "MCL/GraphColor.hpp"
#include "MCL/TetMesh.hpp"
#include "MCL/MeshIO.hpp"
#include <sys/time.h>

using namespace mcl;
typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMat;

bool test_springs();
bool test_bunny();
bool test_dillo();
int make_tri_constraints( mcl::TriangleMesh *mesh, std::vector< Eigen::Triplet<float> > &trips );
int make_tet_constraints( mcl::TetMesh *mesh, std::vector< Eigen::Triplet<float> > &trips );

int main(void){
	srand (time(NULL));
	if( !test_springs() ){ return EXIT_FAILURE; }
	if( !test_dillo() ){ return EXIT_FAILURE; }
	if( !test_bunny() ){ return EXIT_FAILURE; }
	return EXIT_SUCCESS;
}

bool test_bunny(){

	std::cout << "Testing bunny" << std::endl;

	mcl::TriangleMesh bunny;
	std::stringstream bunnyfile;
	bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny.obj";
	mcl::meshio::load_obj( &bunny, bunnyfile.str() );

	std::cout << "Tri Bunny has " << bunny.vertices.size() << " verts" << std::endl;

	std::vector< Eigen::Triplet<float> > triplets;
	int rows = make_tri_constraints( &bunny, triplets );
	int dof = bunny.vertices.size()*3;
	SparseMat A( rows, dof );
	A.setFromTriplets(triplets.begin(), triplets.end());
	A = A.transpose()*A;

	// Try it with stride 1 and stride 3
	for( int stride = 1; stride < 4; stride+=2 ){

		std::cout << "Stride: " << stride << std::endl;

		std::vector< std::vector<int> > colors;
		std::vector<graphcolor::GCNode> nodes;
		graphcolor::make_directed_nodes( A, nodes, stride );
		graphcolor::color_nodes( nodes );
		graphcolor::make_map( nodes, colors );

		// Make sure neighbor nodes don't have the same color
		int n_nodes = nodes.size();
		for( int i=0; i<n_nodes; ++i ){
			graphcolor::GCNode *node = &nodes[i];
			int n_neighbors = node->neighbors.size();
			for(int j=0; j<n_neighbors; ++j ){
				int ca = node->color;
				int cb = nodes[ node->neighbors[j] ].color; 
				if( ca == cb ){
					std::cerr << "**Error: Neighbors with same color" << std::endl;
					return false;
				}
			}
		}

		// Make sure no tris have a node that share a color
		int n_tris = bunny.faces.size();
		for( int i=0; i<n_tris; ++i ){
			Vec3i tri = bunny.faces[i];
			for( int j=0; j<3; ++j ){
			for( int k=0; k<3; ++k ){
			if( stride == 1 ){
			for( int axis=0; axis<3; ++axis ){
				if( j==k ){ continue; }
				int a = nodes[ tri[j]*3+axis ].color;
				int b = nodes[ tri[k]*3+axis ].color;
				if( a == b ){
					std::cerr << "**Error: Coloring failed when checking triangles" << std::endl;
					return false;
				}
			} // end loop axis
			} // end stride 1
			else if( stride == 3 ){
				if( j==k ){ continue; }
				int a = nodes[ tri[j] ].color;
				int b = nodes[ tri[k] ].color;
				if( a == b ){
					std::cerr << "**Error: Coloring failed when checking triangles" << std::endl;
					return false;
				}
			} // end stride 3
			} // end loop k
			} // end loop j
		} // end loop i

		int total_verts = 0;
		int n_colors = colors.size();
		for( int i=0; i<n_colors; ++i ){
			total_verts += colors[i].size();
			std::cout << "\tColor " << i << " has " << colors[i].size() << " nodes" << std::endl;
		}

		if( (total_verts != (int)bunny.vertices.size() && stride == 3 ) ||
			(total_verts != int(bunny.vertices.size()*3) && stride == 1 ) ){
			std::cerr << "**Error: Color groups have more vertices than the mesh" << std::endl;
			return false;
		}

	} // end loop stride

	return true;
}

bool test_dillo(){

	std::cout << "Testing dillo" << std::endl;

	// Load the dillo mesh
	mcl::TetMesh dillo;
	std::stringstream dillofile;
	dillofile << MCLSCENE_ROOT_DIR << "/src/data/armadillo_10k";
	mcl::meshio::load_elenode( &dillo, dillofile.str() );

	// Create the adjacency matrix
	int n_verts = dillo.vertices.size();
	std::cout << "Tet Dillo has " << n_verts << " verts" << std::endl;
	std::vector< Eigen::Triplet<float> > triplets;
	std::vector<graphcolor::GCNode> tempnodes;
	int rows = make_tet_constraints( &dillo, triplets );
	SparseMat A( rows, n_verts*3 );
	A.setFromTriplets(triplets.begin(), triplets.end());
	A = A.transpose()*A;

	// Try with stride 1 and stride 3
	for( int stride = 1; stride < 4; stride+=2 ){

		std::cout << "Stride: " << stride << std::endl;

		// Color the matrix
		std::vector< std::vector<int> > colors;
		std::vector<graphcolor::GCNode> nodes;
		graphcolor::make_directed_nodes( A, nodes, stride );
		graphcolor::color_nodes( nodes );
		graphcolor::make_map( nodes, colors );

		// Make sure no tets have a node that share a color
		int n_tets = dillo.tets.size();
		for( int i=0; i<n_tets; ++i ){
			Vec4i tet = dillo.tets[i];
			for( int j=0; j<4; ++j ){
			for( int k=0; k<4; ++k ){
			if( stride == 1 ){
			for( int axis=0; axis<3; ++axis ){
				if( j==k ){ continue; }
				int a = nodes[ tet[j]*3+axis ].color;
				int b = nodes[ tet[k]*3+axis ].color;
				if( a == b ){
					std::cerr << "**Error: Coloring failed when checking tets" << std::endl;
					return false;
				}
			} // end loop axis
			} // end stride 1
			else if( stride == 3 ){
				if( j==k ){ continue; }
				int a = nodes[ tet[j] ].color;
				int b = nodes[ tet[k] ].color;
				if( a == b ){
					std::cerr << "**Error: Coloring failed when checking tets" << std::endl;
					return false;
				}
			} // end stride 3
			} // end loop k
			} // end loop j
		} // end loop i


		// Make sure neighbor nodes don't have the same color
		int n_nodes = nodes.size();
		for( int i=0; i<n_nodes; ++i ){
			graphcolor::GCNode *node = &nodes[i];
			int n_neighbors = node->neighbors.size();
			for(int j=0; j<n_neighbors; ++j ){
				int ca = node->color;
				int cb = nodes[ node->neighbors[j] ].color; 
				if( ca == cb ){
					std::cerr << "**Error: Neighbors with same color" << std::endl;
					return false;
				}
			}
		}

		// Compute the total number of verts found in the colors
		int n_colors = colors.size();
		int total_verts = 0;
		for( int i=0; i<n_colors; ++i ){
			std::cout << "\tColor " << i << " has " << colors[i].size() << " nodes" << std::endl;
			total_verts += colors[i].size();
		}

		// Test to make sure the number of verts is right
		if( (total_verts != (int)dillo.vertices.size() && stride == 3 ) ||
			(total_verts != int(dillo.vertices.size()*3) && stride == 1 ) ){
			std::cerr << "**Error: Color groups have more vertices than the mesh" << std::endl;
			return false;
		}

	}

	return true;
}

bool test_springs(){

	std::cout << "Testing springs" << std::endl;

	// Do many rounds of coloring. Since it's random it won't return the same
	// result every time, and we want to make sure...
	for( int round=0; round<20; ++round ){

		// Example scene with springs between nodes:
		// (0) --- (1) --- (2) --- (3)
		//      1       2       3 

		// Create an adjacency matrix, e.g. from springs
		int stride = round%3+1; // cycle strides as well
		int n_nodes = 4;
		int dof = n_nodes*stride;
		std::vector< std::pair<int,int> > springs;
		springs.emplace_back( std::make_pair( 0, 1 ) );
		springs.emplace_back( std::make_pair( 1, 2 ) );
		springs.emplace_back( std::make_pair( 2, 3 ) );

		std::vector< Eigen::Triplet<float> > triplets;
		for( int i=0; i<(int)springs.size(); ++i ){
			for( int j=0; j<stride; ++j ){ // 3D springs
				int row = i*stride+j;
				triplets.emplace_back( row, springs[i].first*stride+j, 1 );
				triplets.emplace_back( row, springs[i].second*stride+j, -1 );
			}
		}
		SparseMat A( springs.size()*stride, dof );
		A.setFromTriplets( triplets.begin(), triplets.end() );

		// Make symmetric and color
		A = A.transpose()*A;
		std::vector< std::vector<int> > colors;
		graphcolor::color_matrix( A, colors, stride );

		if( (int)colors.size() == n_nodes ){ continue; } // starting palette is like 6 so...
		for( int i=0; i<(int)colors.size(); ++i ){

			int n_inds = colors[i].size();

			// Do we have indices?
			if( n_inds == 0 ){
				std::cerr << "A color has no indices!" << std::endl;
				return false;
			}

			if( n_inds == 1 ){ continue; }

			if( n_inds == 2 ){
				// Sorting indices to make the check easier
				int idx0 = std::min( colors[i][0], colors[i][1] );
				int idx1 = std::max( colors[i][0], colors[i][1] );
				if( idx0 == idx1 ){
					std::cerr << "An index was added twice!" << std::endl;
					return false;
				}

				if(	(idx0 == 0 && idx1 == 1) || // spring 1
					(idx0 == 1 && idx1 == 2) || // spring 2
					(idx0 == 2 && idx1 == 3) // spring 2
				){
					std::cerr << "Coloring failed: " << idx0 << " and " << idx1 << std::endl;
					return false;
				}
			}

			// Should never have 3 or more with the above setup
			if( n_inds > 2 ){
				std::cerr << "Coloring failed: num inds " << n_inds << std::endl;
				return false;
			}

		} // end loop colors
	}

	return true;
}

//
// This function creates a selection/reduction matrix as used in projective-dynamics
//
int make_tri_constraints( mcl::TriangleMesh *mesh, std::vector< Eigen::Triplet<float> > &trips ){

	using namespace Eigen;
	int n_tris = mesh->faces.size();
	int rows = 0;
	for( int i=0; i<n_tris; ++i ){

		int id0 = mesh->faces[i][0];
		int id1 = mesh->faces[i][1];
		int id2 = mesh->faces[i][2];
		Vec3f x1( mesh->vertices[id0] );
		Vec3f x2( mesh->vertices[id1] );
		Vec3f x3( mesh->vertices[id2] );

		Matrix<float,3,2> D;
		D(0,0) = -1; D(0,1) = -1;
		D(1,0) =  1; D(1,1) =  0;
		D(2,0) =  0; D(2,1) =  1;
		Vec3f e12 = x2 - x1;
		Vec3f e13 = x3 - x1;
		Vec3f n1 = e12.normalized();
		Vec3f n2 = (e13 - e13.dot(n1)*n1).normalized();
		Matrix<float,3,2> basis;
		Matrix<float,3,2> edges;
		basis.col(0) = n1; basis.col(1) = n2;
		edges.col(0) = e12; edges.col(1) = e13;
	
		Matrix<float,2,2> Xg = (basis.transpose() * edges);
		Matrix<float,3,3> X123;
		X123.col(0) = x1; X123.col(1) = x2; X123.col(2) = x3;
		Matrix<float,3,2> B = D * Xg.inverse();

		int cols[3] = { 3*id0, 3*id1, 3*id2 };
		for( int i=0; i<3; ++i ){
			for( int j=0; j<3; ++j ){
				trips.emplace_back( i+rows, cols[j]+i, B(j,0) );
				trips.emplace_back( 3+i+rows, cols[j]+i, B(j,1) );
			}
		}	

		rows += 6;

	} // end loop tris

	return rows;
}


//
// This function creates a selection/reduction matrix as used in projective-dynamics
//
int make_tet_constraints( mcl::TetMesh *mesh, std::vector< Eigen::Triplet<float> > &trips ){

	using namespace Eigen;
	int n_tets = mesh->tets.size();
	int rows = 0;
	for( int i=0; i<n_tets; ++i ){

		int idx[4] = {
			mesh->tets[i][0],
			mesh->tets[i][1],
			mesh->tets[i][2],
			mesh->tets[i][3]
		};

		Matrix<float,3,3> edges;
		Vec3f v0( mesh->vertices[idx[0]] );
		Vec3f v1( mesh->vertices[idx[1]] );
		Vec3f v2( mesh->vertices[idx[2]] );
		Vec3f v3( mesh->vertices[idx[3]] );
		edges.col(0) = v1 - v0;
		edges.col(1) = v2 - v0;
		edges.col(2) = v3 - v0;

		// Xg is beginning state		
		Matrix<float,3,3> Xg = edges;
		Matrix<float,4,3> D;
		D(0,0) = -1; D(0,1) = -1; D(0,2) = -1;
		D(1,0) =  1; D(1,1) =  0; D(1,2) =  0;
		D(2,0) =  0; D(2,1) =  1; D(2,2) =  0;
		D(3,0) =  0; D(3,1) =  0; D(3,2) =  1;

		// B is used to create Ai
		Matrix<float,4,3> B = D * Xg.inverse();
		Matrix<float,3,4> Bt = B.transpose();
		const int col0 = 3 * idx[0];
		const int col1 = 3 * idx[1];
		const int col2 = 3 * idx[2];
		const int col3 = 3 * idx[3];
		const int currrows[3] = { 0, 3, 6 };
		const int cols[4] = { col0, col1, col2, col3 };
		for( int r=0; r<3; ++r ){
			for( int c=0; c<4; ++c ){
				double value = Bt(r,c);
				for( int j=0; j<3; ++j ){
					trips.emplace_back( currrows[r]+j+rows,cols[c]+j, value );
				}
			}
		}

		rows += 9;
	}

	return rows;

}

