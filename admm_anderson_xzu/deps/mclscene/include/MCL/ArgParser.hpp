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

#ifndef MCL_ARGPARSER_H
#define MCL_ARGPARSER_H 1

#include <sstream>
#include <string>
#include <unordered_map>

//
// MCLSCENE ArgParser does simple space-delimited argument parsing
//
// To run a sample program:
//
//	./test_parser --testDouble 0.3 --someflag --testString helloworld --testInt 3 --testFloat 0.33
//
// To load your arguments in the source:
//
//	mcl::ArgParser aParser( argc, argv );
//	double testDouble = aParser.get<double>("--testDouble");
//	float testFloat = aParser.get<float>("--testFloat");
//	int testInt = aParser.get<int>("--testInt");
//	std::string testString = aParser.get<std::string>("--testString");
//	bool someflag_set = aParser.exists("--someflag");
//
// Or if you have defaults that you want to overwrite only if the argument is given:
//
//	double overwriteMe = 3.0;
//	aParser.get<double>("--overwriteMe", &overwriteMe);
//
// And each value will be:
//
//	testDouble: 0.3
//	testFloat: 0.33
//	testInt: 3
//	testString: helloworld
//	someflag_set: true
//	overwriteMe: 3.0
//


namespace mcl {

class ArgParser {
private:
	// Are there faster string-optimized hash tables?
	std::unordered_map< std::string, std::string > args;

public:
	ArgParser( const int &argc, char** argv ) {
		for( int i=1; i<argc-1; ++i ){ args[ argv[i] ] = argv[i+1]; }
		if( argc>1 ){ args[ argv[argc-1] ] = "0"; }
	}

	// Return whether or not an argument exists
	inline bool exists( const std::string label ) const { return ( args.count( label ) > 0 ); }

	// Return the value of the argument (or zero, if not exists)
	template< typename T > const T get( const std::string label ) const {
		if( exists(label) ){
			std::stringstream ss; ss << args.at( label );
			T value; ss >> value;
			return value;
		}
		return T();
	} // end getter

	// Return true on exists and overwrites value, false otherwise.
	template< typename T > bool get( const std::string label, T *result ) const {
		if( exists(label) ){
			std::stringstream ss; ss << args.at( label );
			ss >> *result; return true;
		}
		return false;
	} // end getter to reference
};


}; // end namespace mcl

#endif
