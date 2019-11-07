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

//
// This class provides sint3 and sint2 types along with hashes so they can be used as
// keys in std::unordered_map. The sint2/3 values are sorted, so only unique pairs/triplets
// are inserted. This allows the counting/finding of outer edges/surfaces of meshes.
//

#ifndef MCL_HASHKEYS_H
#define MCL_HASHKEYS_H 1

#include <unordered_map>

namespace mcl {
namespace hashkey {

	template<typename T>
	void swap_if_greater(T& a, T& b){
		if (a > b){
			T tmp(a);
			a = b;
			b = tmp;
		}
	}

	template<typename T>
	void sort(T& a, T& b, T& c){
		swap_if_greater(a, b);
		swap_if_greater(a, c);
		swap_if_greater(b, c);
	}


	struct sint3 {
		sint3(){}
		sint3( int a, int b, int c ){
			sorted_v[0]=a; sorted_v[1]=b; sorted_v[2]=c;
			swap_if_greater<int>(sorted_v[0], sorted_v[1]);
			swap_if_greater<int>(sorted_v[0], sorted_v[2]);
			swap_if_greater<int>(sorted_v[1], sorted_v[2]);
			orig_v[0]=a; orig_v[1]=b; orig_v[2]=c;
		}
		bool operator==(const sint3 &a) const {
			return (sorted_v[0] == a.sorted_v[0]
				&& sorted_v[1] == a.sorted_v[1]
				&& sorted_v[2] == a.sorted_v[2]
			);
		}
		int operator[](const int i) const {
			return orig_v[i];
		}
		int sorted_v[3]; // vertices SORTED
		int orig_v[3]; // original vertices
	};

	struct sint2 {
		sint2(){}
		sint2( int a, int b ){
			sorted_v[0]=a; sorted_v[1]=b;
			if( b < a ){ sorted_v[0]=b; sorted_v[1]=a; }
			orig_v[0]=a; orig_v[1]=b;
		}
		bool operator==(const sint2 &a) const {
			return (sorted_v[0] == a.sorted_v[0]
				&& sorted_v[1] == a.sorted_v[1]
			);
		}
		int operator[](const int i) const {
			return orig_v[i];
		}
		int sorted_v[2]; // vertices SORTED
		int orig_v[2]; // original vertices
	};

} // end namespace hashkey
} // end namespace mcl

namespace std {
	template <> struct hash<mcl::hashkey::sint3> {
		size_t operator()(const mcl::hashkey::sint3& v) const	{
			int a[3] = { v.sorted_v[0], v.sorted_v[1], v.sorted_v[2] };
			unsigned char *in = reinterpret_cast<unsigned char*>(a);
			unsigned int ret = 2654435761u;
			for(unsigned int i = 0; i < (3 * sizeof(int)); ++i)
				ret = (ret * 2654435761u) ^ *in++;
			return ret;
		}
	};
	template <> struct hash<mcl::hashkey::sint2> {
		size_t operator()(const mcl::hashkey::sint2& v) const	{
			int a[2] = { v.sorted_v[0], v.sorted_v[1] };
			unsigned char *in = reinterpret_cast<unsigned char*>(a);
			unsigned int ret = 2654435761u;
			for(unsigned int i = 0; i < (2 * sizeof(int)); ++i)
				ret = (ret * 2654435761u) ^ *in++;
			return ret;
		}
	};
}

#endif
