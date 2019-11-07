// Copyright (c) 2017, University of Minnesota
//
// ADMM-Elastic Uses the BSD 2-Clause License (http://www.opensource.org/licenses/BSD-2-Clause)
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

#ifndef ADMM_COLLIDER_HPP
#define ADMM_COLLIDER_HPP

#include <memory>
#include <iostream>
#include <Eigen/Sparse>

namespace admm {

//
//	Dynamic collision object
//
class DynamicCollision {
public:
	typedef Eigen::Matrix<double,3,1> Vec3;
	typedef Eigen::Matrix<int,3,1> Vec3i;
	typedef Eigen::Matrix<int,4,1> Vec4i;
	typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;

	struct Payload {
		// Input:
		int vert_idx; // colliding vert index
		Vec4i self_tet; // if volumetric, the colliding tet

		// Set by payload:
		double dx; // current lowest signed dist
		Vec3 normal; // collision normal
		Vec3i face; // indices of the collided face
		Vec3 barys; // bary coords of the collision point
		Payload( int idx ) : vert_idx(idx), self_tet(-1,-1,-1,-1),
			dx(std::numeric_limits<double>::max()),
			normal(0,0,0), face(-1,-1,-1), barys(0,0,0) {}
	};

	// Update internal accel structures
	virtual void update( const VecX &x ) = 0;

	// +: no collision
	// 0: on surface
	// -: inside object
	virtual void signed_distance( const Vec3 &x, Payload &p ) const = 0;
};

//
//	Passive collision object
//
class PassiveCollision {
public:
	typedef Eigen::Matrix<double,3,1> Vec3;

	struct Payload {
		int vert_idx; // colliding vert index
		double dx; // current lowest signed dist
		Vec3 point; // point of collision (if needed)
		Vec3 normal; // collision normal (if needed)
		Payload( int idx ) : vert_idx(idx), dx(std::numeric_limits<double>::max()),
			point(0,0,0), normal(0,0,0) {}
	};

	// +: no collision
	// 0: on surface
	// -: inside object
	virtual void signed_distance( const Vec3 &x, Payload &p ) const = 0;
};

//
//	Collision detection interface
//
class Collider {
private:
	typedef Eigen::Matrix<double,3,1> Vec3;
	typedef Eigen::Matrix<int,3,1> Vec3i;
	typedef Eigen::Matrix<double,4,1> Vec4;
	typedef Eigen::Matrix<int,4,1> Vec4i;
	typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMat;
	SparseMat m_C;
	VecX m_c;

public:
	inline void clear_hits(){
		passive_hits.clear();
		dynamic_hits.clear();
		m_C.resize(1,1);
		m_c = VecX::Zero(1);
	}

	// Do a single vertex collision detection test against passive obstacles.
	// Returns true if intersecting something.
	inline bool detect_passive( int idx, const Vec3 &x, Vec3 &n, Vec3 &p );

	// Returns true if collisions have been detected after a detect_whatever call.
	inline bool has_collisions() const { return bool(passive_hits.size()) || bool(dynamic_hits.size()); }

	// Do a round of collision detection and populate the hits vectors.
	// If dynamic objects are added to the collider, the trees are updated with the current vertices.
	// If size of inds=0, all verts are detected.
	inline void detect( const std::vector<int> &inds, const VecX &x, bool with_passive=true );

	// Moves all vertex locations to a state of non-penetration.
	inline void resolve( VecX &x );

	// Passive collisions:
	inline void add_passive_obj( std::shared_ptr<PassiveCollision> obj ){ passive_objs.emplace_back( obj ); }
	const std::vector<PassiveCollision::Payload> &get_passive_hits() const { return passive_hits; }

	// Self collisions:
	inline void add_dynamic_obj( std::shared_ptr<DynamicCollision> obj ){ dynamic_objs.emplace_back( obj ); }
	const std::vector<DynamicCollision::Payload> &get_dynamic_hits() const { return dynamic_hits; }

	// Data:
	std::vector< std::shared_ptr<PassiveCollision> > passive_objs;
	std::vector< std::shared_ptr<DynamicCollision> > dynamic_objs;
	std::vector<PassiveCollision::Payload> passive_hits;
	std::vector<DynamicCollision::Payload> dynamic_hits;
};

inline bool Collider::detect_passive( int idx, const Vec3 &x, Vec3 &n, Vec3 &p ){
	const int num_passive = passive_objs.size();
	if( !num_passive ){ return false; }
	PassiveCollision::Payload p_payload(idx);
	for( int j=0; j<num_passive; ++j ){
		passive_objs[j]->signed_distance(x, p_payload);
		if( p_payload.dx < 0 ){
			n = p_payload.normal;
			p = p_payload.point;
			return true;
		}
	}
	return false;
}

inline void Collider::detect( const std::vector<int> &inds, const VecX &x, bool with_passive ){

	const int num_passive = passive_objs.size();
	const int num_dynamic = dynamic_objs.size();
	if( !num_passive && !num_dynamic ){ return; }
	bool detect_all = inds.size()==0;

	// Thread local results:
	const int nt = omp_get_max_threads();
	std::vector< std::vector<PassiveCollision::Payload> > tl_phits( nt, std::vector<PassiveCollision::Payload>() );
	std::vector< std::vector<DynamicCollision::Payload> > tl_dhits( nt, std::vector<DynamicCollision::Payload>() );
	int n_verts = detect_all ? x.rows()/3 : inds.size();

	// Update trees in dynamic object
	#pragma omp parallel for
	for( int i=0; i<num_dynamic; ++i ){ dynamic_objs[i]->update(x); }

	#pragma omp parallel for
	for( int i=0; i<n_verts; ++i ){

		int idx = i;
		if( !detect_all ){ idx = inds[i]; }
		Vec3 curr_x = x.segment<3>(idx*3);

		//
		// check passive objects
		//
		if( with_passive ){
			PassiveCollision::Payload p_payload(idx);
			for( int j=0; j<num_passive; ++j ){
				passive_objs[j]->signed_distance(curr_x, p_payload);
			}
			if( p_payload.dx < 0 ){
				int thread_idx = omp_get_thread_num();
				tl_phits[ thread_idx ].emplace_back(p_payload);
			}
		}

		//
		// check dynamic objects
		//
		DynamicCollision::Payload d_payload(idx);
		for( int j=0; j<num_dynamic; ++j ){
			dynamic_objs[j]->signed_distance(curr_x, d_payload);
		}
		if( d_payload.dx < 0 ){
			int thread_idx = omp_get_thread_num();
			tl_dhits[ thread_idx ].emplace_back(d_payload);
		}

	} // end per vertex collision

	// Combine thread-local results
	for( int i=0; i<nt; ++i ){
		std::vector<PassiveCollision::Payload> *tlrp = &tl_phits[i];
		passive_hits.insert( std::end(passive_hits), std::begin(*tlrp), std::end(*tlrp) );
		std::vector<DynamicCollision::Payload> *tlrd = &tl_dhits[i];
		dynamic_hits.insert( std::end(dynamic_hits), std::begin(*tlrd), std::end(*tlrd) );
	}

} // end detect

} // namespace admm

#endif
