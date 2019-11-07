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

#include "MCL/MeshIO.hpp"
#include "MCL/RenderWindow.hpp"
#include "MCL/ShapeFactory.hpp"

using namespace mcl;

class CustomController : public Controller {
public:
	void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
		Controller::key_callback(window,key,scancode,action,mods); // call base class cb if you want
		std::cout << "KEY PRESSED: " << key << std::endl; // Print keypress for fun
	}
};

int main(void){

	// Load the mesh
	std::shared_ptr<mcl::TriangleMesh> bunny1 = mcl::TriangleMesh::create();
	std::shared_ptr<mcl::TriangleMesh> bunny2 = mcl::TriangleMesh::create();
	std::shared_ptr<mcl::TetMesh> box = mcl::factory::make_tet_blocks( 1, 2, 1 );


	std::stringstream bunnyfile;
	bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny.obj";
	mcl::meshio::load_obj( bunny1.get(), bunnyfile.str() );
	mcl::meshio::load_obj( bunny2.get(), bunnyfile.str() );

	// Move the second bunny lower and both bunnies larger
	mcl::XForm<float> trans = mcl::xform::make_trans(0.f,-0.15f,0.f);
	mcl::XForm<float> scale = mcl::xform::make_scale(10.f,10.f,10.f);
	bunny1->apply_xform( scale );
	bunny2->apply_xform( scale*trans );

	// Create a render mesh
	std::shared_ptr<mcl::RenderMesh> rmBunny1 = mcl::RenderMesh::create( bunny1 );
	std::shared_ptr<mcl::RenderMesh> rmBunny2 = mcl::RenderMesh::create( bunny2 );
	std::shared_ptr<mcl::RenderMesh> rmBox = mcl::RenderMesh::create( box, mcl::RenderMesh::WIREFRAME );

	// Create opengl context
	RenderWindow renderWindow;
	renderWindow.set_controller( std::make_shared<CustomController>() );
	GLFWwindow* window = renderWindow.init();
	if( !window ){ return EXIT_FAILURE; }

	// Add the render mesh to the window
	renderWindow.add_mesh( rmBunny1 );
	renderWindow.add_mesh( rmBunny2 );
	renderWindow.add_mesh( rmBox );

	// Set the camera to a nice location based
	// on the meshes added to the renderWindow.
	renderWindow.nice_camera_location();

	// Game loop
	while( renderWindow.is_open() ){

		renderWindow.draw();
		glfwPollEvents();

	} // end game loop

	return 0;
}

