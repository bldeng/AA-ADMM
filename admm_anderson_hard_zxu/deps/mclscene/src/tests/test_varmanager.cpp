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
#include "MCL/VarManager.hpp"

using namespace mcl;

class TestVec {
public:
	TestVec(){} // required for VarManager
	TestVec( float x_, float y_, float z_ ) : x(x_), y(y_), z(z_) {}
	float x, y, z;

	// required for VarManager
	friend std::ostream &operator <<(std::ostream &out, TestVec &t){
		out << t.x << ' '; out << t.y << ' '; out << t.z;
		return out;
	}

	// required for VarManager
	friend std::istream &operator >>(std::istream &in, TestVec &t){
		in >> t.x; in >> t.y; in >> t.z;
		return in;
	}
};

int main(void){

	mcl::VarManager vm;

	// Set some vars
	vm.set( "aDouble", 32.0 );
	vm.set( "aInt", 2 );
	vm.set( "aVec", TestVec(1,2,3) );

	// Get some vars
	{
		double aDouble = vm.get<double>( "aDouble" );
		int aInt = vm.get<int>( "aInt" );
		TestVec aVec = vm.get<TestVec>( "aVec" );
		float aFloat = 123.f;
		bool found = vm.get<float>("aFloat_", aFloat );

		// Check
		if( aDouble != 32.0 ){
			std::cerr << "Failed to get double" << std::endl;
			return EXIT_FAILURE;
		}

		if( aInt != 2 ){
			std::cerr << "Failed to get int" << std::endl;
			return EXIT_FAILURE;
		}

		if( found || aFloat != 123.f ){
			std::cerr << "Failed to get float" << std::endl;
			return EXIT_FAILURE;
		}

		if( aVec.x != 1 || aVec.y != 2 || aVec.z != 3 ){
			std::cerr << "Failed to get custom type" << std::endl;
			return EXIT_FAILURE;
		}
	}

	// Test write
	vm.write( "VarManagerOut.txt" );

	// Test read
	vm.clear();
	if( vm.exists( "aDouble" ) ){
		std::cerr << "Should be empty after clear" << std::endl;
		return EXIT_FAILURE;
	}
	vm.read( "VarManagerOut.txt" );

	// Check
	{
		// Get some vars
		double aDouble = vm.get<double>( "aDouble" );
		int aInt = vm.get<int>( "aInt" );
		TestVec aVec = vm.get<TestVec>( "aVec" );
		float aFloat = 123.f;
		bool found = vm.get<float>("aFloat_", aFloat );

		// Check
		if( aDouble != 32.0 ){
			std::cerr << "Failed to get double" << std::endl;
			return EXIT_FAILURE;
		}

		if( aInt != 2 ){
			std::cerr << "Failed to get int" << std::endl;
			return EXIT_FAILURE;
		}

		if( found || aFloat != 123.f ){
			std::cerr << "Failed to get float" << std::endl;
			return EXIT_FAILURE;
		}

		if( aVec.x != 1 || aVec.y != 2 || aVec.z != 3 ){
			std::cerr << "Failed to get custom type" << std::endl;
			return EXIT_FAILURE;
		}
	}


	std::cout << "Success" << std::endl;
	return EXIT_SUCCESS;
}
