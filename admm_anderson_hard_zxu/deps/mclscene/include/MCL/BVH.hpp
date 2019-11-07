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

#ifndef MCL_BVH_H
#define MCL_BVH_H 1

#include "Visitor.hpp"
#include <memory>
#include <numeric>

namespace mcl {
namespace bvh {

// Binary AABB Tree
// PDIM is the dimension of the primitive,
// i.e. 1=verts, 2=edges, 3=tris, 4=tets
template <typename T, short PDIM>
class AABBTree {
typedef Eigen::AlignedBox<T,3> AABB;
public:
	AABBTree();

	// Create a tree from a list of primitives
	void init( const int *inds, const T *verts, int num_prims );

	// Traverse the tree with a visitor.
	// Returns the result of Visitor::hit_<whatever>()
	// See MCL/Visitor.hpp
	bool traverse( Visitor<T,PDIM> &visitor ) const {
		return traverse_children( root_node.get(), visitor );
	}

private:
	struct Node {
		AABB aabb;
		Node *left, *right;
		int prim;
		Node() : left(nullptr), right(nullptr), prim(-1) {}
		~Node(){
			if( left != nullptr ){ delete left; }
			if( right != nullptr ){ delete right; }
		}
	};

	static void create_children(
		Node *node,
		std::vector<int> &queue,
		const std::vector< AABB > &leaves,
		const std::vector< Vec3<T> > &centroids );

	static bool traverse_children( const Node *node,
		Visitor<T,PDIM> &visitor );

	std::unique_ptr<Node> root_node;

}; // class aabbtree



//
//	Implementation
//


template <typename T, short PDIM>
AABBTree<T,PDIM>::AABBTree(){
	root_node.reset( new Node() );
	if( PDIM < 1 ){
		throw std::runtime_error("AABBTree Error: PDIM must be larger than 0");
	}
}


template <typename T, short PDIM>
void AABBTree<T,PDIM>::init( const int *inds, const T *verts, int num_prims ){

	// Deletes the old tree
	root_node.reset( new Node() );

	// Leaf nodes are copied into the tree, but we'll create
	// them here to make processing faster.
	std::vector< AABB > leaf_aabb( num_prims );
	std::vector< Vec3<T> > leaf_centroids( num_prims, Vec3<T>(0,0,0) );

	// Create leaf AABBS
	#pragma omp parallel for
	for( int i=0; i<num_prims; ++i ){
		for( int j=0; j<PDIM; ++j ){
			int prim_id = inds[i*PDIM+j];
			Vec3<T> p( verts[prim_id*3], verts[prim_id*3+1], verts[prim_id*3+2] );
			leaf_aabb[i].extend( p );
			leaf_centroids[i] += p;
		}
		leaf_centroids[i] /= T(PDIM);
	}

	for( int i=0; i<num_prims; ++i ){ root_node->aabb.extend( leaf_aabb[i] ); }

	// Now do a recursive top down creation
	std::vector<int> queue(num_prims);
	std::iota(queue.begin(), queue.end(), 0);
	create_children(root_node.get(), queue, leaf_aabb, leaf_centroids);

} // end init


template <typename T, short PDIM>
bool AABBTree<T,PDIM>::traverse_children( const Node *node, Visitor<T,PDIM> &visitor ){

	if( !visitor.hit_aabb( node->aabb ) ){ return false; }

	// If it has two children, see which one we should traverse first
	if( node->left != nullptr && node->right != nullptr ){
		bool check_left_first = visitor.check_left_first( node->left->aabb, node->right->aabb );
		if( check_left_first ){
			if( traverse_children( node->left, visitor ) ){ return true; }
			else { return traverse_children( node->right, visitor ); }
		} else {
			if( traverse_children( node->right, visitor ) ){ return true; }
			else { return traverse_children( node->left, visitor ); }
		}	
	}

	// If we just have a left child
	if( node->left != nullptr ){ return traverse_children( node->left, visitor ); }

	// If we just have a right child
	if( node->right != nullptr ){ return traverse_children( node->right, visitor ); }

	// Otherwise we're a leaf!
	if( node->prim == -1 ){
		throw std::runtime_error("AABBTree::traverse Error: Leaf has no primitive");
	}
	return visitor.hit_prim( node->prim );
}


template <typename T, short PDIM>
void AABBTree<T,PDIM>::create_children(
	Node *node,
	std::vector<int> &queue,
	const std::vector< AABB > &leaves,
	const std::vector< Vec3<T> > &centroids ){

	int n_queue = queue.size();
	if( n_queue == 0 ){
		throw std::runtime_error("AABBTree::init Error: Empty queue");
	}

	// One element means we are a leaf
	if( n_queue == 1 ){
		int qidx = queue[0];
		node->prim = qidx;
		node->aabb = leaves[qidx];
		return;
	}

	// Compute the splitting plane
	AABB tempAABB;
	for( int i=0; i<n_queue; ++i ){ tempAABB.extend( centroids[ queue[i] ] ); }	
	Vec3<T> sides = tempAABB.sizes();
	int split = 0;
	if( sides[1] >= sides[0] && sides[1] >= sides[2] ){ split = 1; }
	else if( sides[2] >= sides[0] && sides[2] >= sides[1] ){ split = 2; }

	// If two elements, make left and right
	if( n_queue == 2 ){

		node->left = new Node();
		node->right = new Node();
		int idx0 = queue[0];
		int idx1 = queue[1];
		const Vec3<T> &cent0 = centroids[idx0];
		const Vec3<T> &cent1 = centroids[idx1];

		if( cent0[split] < cent1[split] ){
			node->left->prim = idx0;
			node->left->aabb = leaves[idx0];
			node->right->prim = idx1;
			node->right->aabb = leaves[idx1];
		} else {
			node->left->prim = idx1;
			node->left->aabb = leaves[idx1];
			node->right->prim = idx0;
			node->right->aabb = leaves[idx0];
		}
		return;
	}

	// Split the queue into left and right
	T center = tempAABB.center()[split];
	AABB left_aabb, right_aabb;
	std::vector<int> left_queue, right_queue;
	for( int i=0; i<n_queue; ++i ){
		int idx = queue[i];
		const Vec3<T> &cent = centroids[idx];
		if( cent[split] < center ){
			left_queue.push_back(idx);
			left_aabb.extend(leaves[idx]);
		} else {
			right_queue.push_back(idx);
			right_aabb.extend(leaves[idx]);
		}
	}

	if( left_queue.size()==0 || right_queue.size()==0 ){
		throw std::runtime_error("AABBTree::init Error: problem splitting geometry");
	}

	// Create the left child
	node->left = new Node();
	node->left->aabb = left_aabb;
	create_children( node->left, left_queue, leaves, centroids );

	// Create the right child
	node->right = new Node();
	node->right->aabb = right_aabb;
	create_children( node->right, right_queue, leaves, centroids );

} // end create childrens


} // end ns bvh
} // end ns mcl

#endif
