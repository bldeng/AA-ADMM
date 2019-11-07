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
//	OpenGL Material Presets
//	See: http://devernay.free.fr/cours/opengl/materials.html
//	Create a blinn phong material from preset values.
//

#ifndef MCL_MATERIAL_H
#define MCL_MATERIAL_H 1

#include "Vec.hpp"

namespace mcl {
namespace material {

enum class Preset {
	Emerald, Jade, Obsidian, Pearl, Ruby, Turquoise, // gems
	Brass, Bronze, Chrome, Copper, Gold, Silver, Gunmetal, // metals
	BlackPlastic, CyanPlastic, GreenPlastic, RedPlastic, WhitePlastic, YellowPlastic, // plastic
	BlackRubber, CyanRubber, GreenRubber, RedRubber, WhiteRubber, YellowRubber, // rubber
	Cloth, Unknown
};

class Phong {
public:
	Vec3f amb, diff, spec;
	float shini; // In constructor, use a value from 0 to 1.

	static Phong create( Preset m );

	Phong( Vec3f amb_, Vec3f diff_, Vec3f spec_, float shini_ )
		: amb(amb_), diff(diff_), spec(spec_), shini(shini_*128.f) {}
	Phong() : amb(1,0,0), diff(1,0,0), spec(0.3f,0.3f,0.3f), shini(32.f) {}
};

// Create a Phong color from a list of my favorite presets
static Phong autoPhong( unsigned int start=0 ){
	const std::vector<Preset> presets = {	
        Preset::Bronze,
        Preset::Emerald,
		Preset::WhitePlastic,
		Preset::GreenRubber,
        Preset::Ruby,
	};
	static int counter = start;
	int p_idx = counter % presets.size();
	counter++;
	return Phong::create( presets[p_idx] );
}

Phong Phong::create( Preset m ){

	Phong r;
	switch( m ){

	// Gemstones
	case Preset::Emerald:
		r = Phong( Vec3f(0.0215, 0.1745, 0.0215), Vec3f(0.07568, 0.61424, 0.07568), Vec3f(0.633, 0.727811, 0.633), 0.6); break;
	case Preset::Jade:
		r = Phong( Vec3f(0.135, 0.2225, 0.1575), Vec3f(0.54, 0.89, 0.63), Vec3f(0.316228, 0.316228, 0.316228), 0.1); break;
	case Preset::Obsidian:
		r = Phong( Vec3f(0.05375, 0.05, 0.06625), Vec3f(0.18275, 0.17, 0.22525), Vec3f(0.332741, 0.328634, 0.346435), 0.3); break;
	case Preset::Pearl:
		r = Phong( Vec3f(0.25, 0.20725, 0.20725), Vec3f(1.0, 0.829, 0.829), Vec3f(0.296648, 0.296648, 0.296648), 0.088); break;
	case Preset::Ruby:
		r = Phong( Vec3f(0.1745, 0.01175, 0.01175), Vec3f(0.61424, 0.04136, 0.04136), Vec3f(0.727811, 0.626959, 0.626959), 0.6); break;
	case Preset::Turquoise:
		r = Phong( Vec3f(0.1, 0.18725, 0.1745), Vec3f(0.396, 0.74151, 0.69102), Vec3f(0.297254, 0.30829, 0.306678), 0.1); break;

	// Metals
	case Preset::Brass:
		r = Phong( Vec3f(0.329412, 0.223529, 0.027451), Vec3f(0.780392, 0.568627, 0.113725), Vec3f(0.992157, 0.941176, 0.807843), 0.21794872); break;
	case Preset::Bronze:
		r = Phong( Vec3f(0.2125, 0.1275, 0.054), Vec3f(0.714, 0.4284, 0.18144), Vec3f(0.393548, 0.271906, 0.166721), 0.2); break;
	case Preset::Chrome:
		r = Phong( Vec3f(0.25, 0.25, 0.25), Vec3f(0.4, 0.4, 0.4), Vec3f(0.774597, 0.774597, 0.774597), 0.6); break;
	case Preset::Copper:
		r = Phong( Vec3f(0.19125, 0.0735, 0.0225), Vec3f(0.7038, 0.27048, 0.0828), Vec3f(0.256777, 0.137622, 0.086014), 0.6); break;
	case Preset::Gold:
		r = Phong( Vec3f(0.24725, 0.1995, 0.0745), Vec3f(0.75164, 0.60648, 0.22648), Vec3f(0.628281, 0.555802, 0.366065), 0.4); break;
	case Preset::Silver:
		r = Phong( Vec3f(0.19225, 0.19225, 0.19225), Vec3f(0.50754, 0.50754, 0.50754), Vec3f(0.508273, 0.508273, 0.508273), 0.4); break;
	case Preset::Gunmetal: // I made this one up
		r = Phong( Vec3f(0.1, 0.1, 0.1), Vec3f(0.4, 0.4, 0.4), Vec3f(0.1, 0.1, 0.1), 0.1); break;

	// Plastics
	case Preset::BlackPlastic:
		r = Phong( Vec3f(0.0, 0.0, 0.0), Vec3f(0.01, 0.01, 0.01), Vec3f(0.50, 0.50, 0.50), 0.25); break;
	case Preset::CyanPlastic:
		r = Phong( Vec3f(0.0, 0.1, 0.06), Vec3f(0.0, 0.50980392, 0.50980392), Vec3f(0.50196078, 0.50196078, 0.50196078), 0.25); break;
	case Preset::GreenPlastic:
		r = Phong( Vec3f(0.0, 0.0, 0.0), Vec3f(0.1, 0.35, 0.1), Vec3f(0.45, 0.55, 0.45), 0.25); break;
	case Preset::RedPlastic:
		r = Phong( Vec3f(0.0, 0.0, 0.0), Vec3f(0.5, 0.0, 0.0), Vec3f(0.7, 0.6, 0.6), 0.25); break;
	case Preset::WhitePlastic:
		r = Phong( Vec3f(0.0, 0.0, 0.0), Vec3f(0.55, 0.55, 0.55), Vec3f(0.70, 0.70, 0.70), 0.25); break;
	case Preset::YellowPlastic:
		r = Phong( Vec3f(0.0, 0.0, 0.0), Vec3f(0.5, 0.5, 0.0), Vec3f(0.60, 0.60, 0.50), 0.25); break;

	// Rubbers
	case Preset::BlackRubber:
		r = Phong( Vec3f(0.02, 0.02, 0.02), Vec3f(0.01, 0.01, 0.01), Vec3f(0.4, 0.4, 0.4), 0.078125); break;
	case Preset::CyanRubber:
		r = Phong( Vec3f(0.0, 0.05, 0.05), Vec3f(0.4, 0.5, 0.5), Vec3f(0.04, 0.7, 0.7), 0.078125); break;
	case Preset::GreenRubber:
		r = Phong( Vec3f(0.0, 0.05, 0.0), Vec3f(0.4, 0.5, 0.4), Vec3f(0.04, 0.7, 0.04), 0.078125); break;
	case Preset::RedRubber:
		r = Phong( Vec3f(0.05, 0.0, 0.0), Vec3f(0.5, 0.4, 0.4), Vec3f(0.7, 0.04, 0.04), 0.078125); break;
	case Preset::WhiteRubber:
		r = Phong( Vec3f(0.05, 0.05, 0.05), Vec3f(0.5, 0.5, 0.5), Vec3f(0.7, 0.7, 0.7), 0.078125); break;
	case Preset::YellowRubber:
		r = Phong( Vec3f(0.05, 0.05, 0.0), Vec3f(0.5, 0.5, 0.4), Vec3f(0.7, 0.7, 0.04), 0.078125); break;
	case Preset::Cloth:
		r = Phong( Vec3f(0.25, 0.20725, 0.20725), Vec3f(1.0, 0.829, 0.829), Vec3f(0., 0., 0.), 0.088); break;

	default: break;

	} // end switch preset

	// Apply some ambient dampening
	r.amb *= 0.2f;

	return r;

} // end material preset

} // end namespace material
} // end namespace mcl

#endif
