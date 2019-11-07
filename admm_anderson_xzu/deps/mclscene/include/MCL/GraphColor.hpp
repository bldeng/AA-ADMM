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

#ifndef MCL_GRAPHCOLOR_H
#define MCL_GRAPHCOLOR_H

#include <Eigen/Sparse>
#include <unordered_set>
#include <cstdlib>
#include <random>
#include "omp.h"

namespace mcl {

// Fast random graph coloring from:
// "Fast Distributed Algorithms for Brooks-Vizing Colourings" by Grable and Panconesi (2000)
// With some minor modifications (mostly that it works on a directed graph).
namespace graphcolor {

	template <typename T> using SparseMat = Eigen::SparseMatrix<T,Eigen::RowMajor>;

	// Interal helper class for setting up the graph and managing colors
	// and node queue (in color_nodes)
	struct GCNode {
		GCNode() : idx(-1), color(-1), conflict(false) {}
		int idx, color; // idx only used for error testing
		bool conflict;
		std::vector<int> neighbors;
		std::unordered_set<int> palette;
	};

	// Creates a directed graph from a symmetric adjacency matrix A (only the upper-triangular is used).
	// Any values in A are considered an edge, no matter the sign/value.
	// So, nodes only see neighbors with a higher index.
	template <typename T>
	static void make_directed_nodes( const SparseMat<T> &A, std::vector<GCNode> &nodes, int stride );

	// Colors the graph created by make_nodes (nodes are modified)
	static void color_nodes( std::vector<GCNode> &nodes );

	// Mapping is color -> list of node indices, e.g.
	//	std::vector<int> color0 = colors[0];
	//	int color0_idx0 = color0[0];
	// where color0_idx0 corresponds to some row in the adjacency matrix
	static void make_map( const std::vector<GCNode> &nodes, std::vector< std::vector<int> > &colors );

	// Colors an adjacency matrix with the above functions
	template <typename T>
	static void color_matrix( const SparseMat<T> &A, std::vector< std::vector<int> > &colors, int stride=1 ){
		std::vector<GCNode> nodes;
		make_directed_nodes( A, nodes, stride );
		color_nodes( nodes );
		make_map( nodes, colors );
	}


} // ns graphcolor

//
// Implementation
//

template <typename T>
static void graphcolor::make_directed_nodes( const SparseMat<T> &A, std::vector<GCNode> &nodes, int stride ){

	// Should also check to make sure it's symmetric but I'll save that for later...
	int dof = A.rows();
	if( dof != A.cols() ){ throw std::runtime_error("**graphcolor::make_nodes Error: Adjacency matrix is not square"); }
	stride = std::max( stride, 1 );
	if( dof % stride != 0 ){ throw std::runtime_error("**graphcolor::make_nodes Error: Bad stride"); }
	const int n_nodes = dof/stride;
	nodes.resize( n_nodes );

	// Create neighbor list
	#pragma omp parallel for schedule(static)
	for( int i=0; i<n_nodes; ++i ){

		GCNode *node = &nodes[i];
		node->idx = i;
		typename Eigen::SparseMatrix<T,Eigen::RowMajor>::InnerIterator it(A,i*stride);

		// Loop until we are past the diagonal
		for( ; it && it.col()/stride <= i; ++it ){}

		// Now add neighbors, avoid adding same neighbor twice
		std::vector<int> already_added = {-1};
		for( ; it; ++it ){
			if( std::abs(it.value()) > T(0) ){
				int idx = it.col()/stride;
				if( already_added.back()==idx ){ continue; }
				already_added.emplace_back(idx);
				node->neighbors.emplace_back(idx);
			}
		}

	} // end loop elements

} // end make nodes


static void graphcolor::color_nodes( std::vector<GCNode> &all_nodes ){

	int init_palette = 6; // based on observations
	int n_nodes = all_nodes.size();
	std::vector<int> nodeq( n_nodes ); // node queue
	std::iota(nodeq.begin(), nodeq.end(), 0);

	// Make queue of nodes that need to be processed
	// and initialize their palette.
	#pragma omp parallel for schedule(static)
	for( int i=0; i<n_nodes; ++i ){
		GCNode *node = &all_nodes[ nodeq[i] ];
		for( int j=0; j<init_palette; ++j ){ node->palette.insert(j); }

		// Make sure graph is directed
		int n_neighbors = node->neighbors.size();
		for( int j=0; j<n_neighbors; ++j ){
			const GCNode *n1 = &all_nodes[ node->neighbors[j] ];
			if( n1->idx == node->idx ){
				std::cerr << "**graphcolor::color_nodes Error: Self is a neighbor" << std::endl;
				throw std::runtime_error("exiting...");
			}
			if( n1->idx < node->idx ){
				std::cerr << "**graphcolor::color_nodes Error: Not a directed graph" << std::endl;
				std::cerr << "\tnode: " << node->idx << ", neighbor: " << n1->idx << std::endl;
				throw std::runtime_error("exiting...");
			}
		}
	}

	// There is a guarantee on max number of loops based on the max degree
	// of the graph. In my tests I rarely hit that, so I'll just pick some large number.
	const int max_iter = all_nodes.size();
	for( int rand_iter=0; nodeq.size()>0 && rand_iter < max_iter; ++rand_iter ){

		#pragma omp parallel
		{
			unsigned int tseed = omp_get_thread_num() * time(NULL);

			// Generate a random color
			#pragma omp for
			for( int i=0; i<n_nodes; ++i ){
				GCNode *node = &all_nodes[ nodeq[i] ];
				// Note about rand_r: posix function, so it may not compile in WIN.
				// Unfortunately there is no std replacement that I know of.
				int c_idx = rand_r(&tseed) % node->palette.size();
				node->color = *std::next( node->palette.begin(), c_idx );
			}

			// Conflict detection
	 		#pragma omp for
			for( int i=0; i<n_nodes; ++i ){
				GCNode *node = &all_nodes[ nodeq[i] ];
				int curr_c = node->color;
				node->conflict = false;

				// Check neighbors
				int n_neighbors = node->neighbors.size();
				for( int j=0; j<n_neighbors && !node->conflict; ++j ){

					// Hungarian heuristic: node with largest index keeps color
					int cn = all_nodes[ node->neighbors[j] ].color;
					if( curr_c == cn ){ node->conflict = true; }
				}

			} // end for-loop conflict resolution

		} // end parallel region

		// Remove color from neighbors
		for( int i=0; i<n_nodes; ++i ){
			GCNode *node = &all_nodes[ nodeq[i] ];
			const int nc = node->color;

			// Look at neighbors. If they are not in conflict,
			// remove their color from your own palette.
			int n_neighbors = node->neighbors.size();
			for( int j=0; j<n_neighbors; ++j ){
				GCNode *n1 = &all_nodes[ node->neighbors[j] ];
				if( !n1->conflict ){ node->palette.erase(n1->color); }

				// This is the part that is not thread safe:
				if( !node->conflict ){ n1->palette.erase(nc); }
			}

		} // end update neighbor palette

		// Remove colored nodes from the queue
		std::vector<int>::iterator node_iter = nodeq.begin();
		for( ; node_iter != nodeq.end(); ){
			GCNode *node = &all_nodes[ *node_iter ];
			if( node->conflict ){ node_iter++; } // keep
			else{ node_iter = nodeq.erase(node_iter); } // remove
		}

		// Feed the hungry
		n_nodes = nodeq.size();
		#pragma omp parallel for schedule(static)
		for( int i=0; i<n_nodes; ++i ){
			GCNode *node = &all_nodes[ nodeq[i] ];
			if( node->palette.size() < 2 ){
				node->palette.insert( init_palette+rand_iter );
			}
		}

	} // end color loop

	// Make sure all nodes are colored
	if( nodeq.size()>0 ){
		throw std::runtime_error("graphcolor::color Error: Nodes remain uncolored");
	}

} // end color


static inline void graphcolor::make_map( const std::vector<GCNode> &nodes, std::vector< std::vector<int> > &colors ){

	// Since the colors are random they are not necessarily 0 to n.
	// So, we'll just remove colors that don't have any nodes but start with a bunch.
	colors.clear();
	for( int i=0; i<14; ++i ){ colors.emplace_back( std::vector<int>() ); }

	int n_nodes = nodes.size();
	for( int i=0; i<n_nodes; ++i ){
		const GCNode *n = &nodes[i];
		// Add colors as needed
		while( n->color >= (int)colors.size() ){ colors.emplace_back( std::vector<int>() ); }
		colors[ n->color ].emplace_back( n->idx );
	}

	// Remove empty colors
	std::vector< std::vector<int> >::iterator it = colors.begin();
	for( ; it != colors.end(); ){ it->size() == 0 ? it = colors.erase(it) : it++; }

} // end make map

} // ns mcl


#endif

