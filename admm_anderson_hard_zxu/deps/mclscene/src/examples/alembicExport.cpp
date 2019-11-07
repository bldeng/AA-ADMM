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
#include "MCL/AlembicIO.hpp"

using namespace mcl;

void centerize( mcl::TriangleMesh::Ptr mesh );
void generate_animation( std::vector< std::vector<mcl::Vec3f> > &verts, mcl::TriangleMesh::Ptr mesh );

//
// May as well let the controller hold on to the data
//
class TimestepController : public Controller {
public:
	int curr_t;
	bool key_down;
	bool needs_reload;
	std::vector< std::vector<mcl::Vec3f> > per_dt_verts; // verts at each time step
	mcl::TriangleMesh::Ptr tmesh;
	mcl::RenderMesh::Ptr rmesh;

	TimestepController() : curr_t(0), key_down(false), needs_reload(true) {
		tmesh = mcl::TriangleMesh::create();

		std::cout << "Loading the mesh..." << std::endl;
		std::stringstream bunnyfile;
		bunnyfile << MCLSCENE_ROOT_DIR << "/src/data/bunny.obj";
		mcl::meshio::load_obj( tmesh.get(), bunnyfile.str() );
		tmesh->apply_xform( mcl::xform::make_scale(10.f,10.f,10.f) ); // Increase size
		centerize( tmesh );

		std::cout << "Generating an animation..." << std::endl;
		rmesh = mcl::RenderMesh::create( tmesh );
		generate_animation( per_dt_verts, tmesh );
	}

	void load_verts_at_t(){
		if( !needs_reload ){ return; }
		int n_v = tmesh->vertices.size();
		for( int i=0; i<n_v; ++i ){ tmesh->vertices[i] = per_dt_verts[curr_t][i]; }
		tmesh->need_normals(true);
		rmesh->load_buffers();
		needs_reload = false;
	}

	void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
		Controller::key_callback(window,key,scancode,action,mods); // Standard key io
		if( action == GLFW_PRESS ){ key_down = true; }
		else if( action == GLFW_RELEASE ){ key_down = false; }
		if( !key_down ){ return; }

		if( key == GLFW_KEY_LEFT ){
			curr_t = std::max( 0, curr_t-1 );
			needs_reload = true;
		}
		else if( key == GLFW_KEY_RIGHT ){
			curr_t = std::min( int(per_dt_verts.size()-1), curr_t+1 );
			needs_reload = true;
		}

	} // end key callback
};

//
// Basic animation config
//
struct AnimationConfig {
	int frame_rate;
	float end_time;
	mcl::XForm<float> xf;
	AnimationConfig() {
		frame_rate = 20;
		end_time = 3;
		float r_speed = 100.f;
		srand(0);
		mcl::Vec3f rand_vel = Eigen::Vector3f::Random();
		xf = mcl::xform::make_rot( 1.f/float(frame_rate) * r_speed, rand_vel );
	}
} config;

int main(void){

	// Create the controller and generate the animation
	std::shared_ptr<TimestepController> tController = std::make_shared<TimestepController>();
	std::cout << "Press left/right arrow keys to step through the animation" << std::endl;

	// Create opengl context
	RenderWindow renderWindow;
	renderWindow.set_controller( tController );
	GLFWwindow* window = renderWindow.init();
	if( !window ){ return EXIT_FAILURE; }
	renderWindow.add_mesh( tController->rmesh );

	// Set the camera to a nice location based
	// on the meshes added to the renderWindow.
	renderWindow.nice_camera_location();

	// Create the exporter, and export the animation to an alembic file
	AlembicExporter exporter( config.frame_rate, "test.abc" );
	int handle2 = exporter.add_object( "bunny" );
	int handle = exporter.add_object( "bunny2" );

	// Loop through the frames computed in init, and save each one
	int n_frames = tController->per_dt_verts.size();
	for( int i=0; i<n_frames; ++i ){
		mcl::TriangleMesh::Ptr tmesh = tController->tmesh;
		float *verts = &tController->per_dt_verts[i][0][0];

		// The frame is advanced every "set_frame" call on the same handle.
		exporter.add_frame( handle, verts, tmesh->vertices.size(), &tmesh->faces[0][0], tmesh->faces.size() );

		// Every other frame we'll set the second bunny, so it looks like it's spinning faster
		if( i%2 == 0 ){
			std::vector<mcl::Vec3f> bunny2 = tController->per_dt_verts[i];
			int n_v = bunny2.size();
			for( int j=0; j<n_v; ++j ){ bunny2[j][1] -= 2.f; } // move below the other bunny
			exporter.add_frame( handle2, &bunny2[0][0], bunny2.size(), &tmesh->faces[0][0], tmesh->faces.size() );
		}
	}

	// Game loop
	while( renderWindow.is_open() ){
		tController->load_verts_at_t();
		renderWindow.draw();
		glfwPollEvents();
	} // end game loop

	return 0;
}


void centerize( mcl::TriangleMesh::Ptr mesh ){
	Eigen::AlignedBox<float,3> aabb = mesh->bounds();
	mcl::Vec3f offset = -aabb.center();
	int n_v = mesh->vertices.size();
	for( int i=0; i<n_v; ++i ){ mesh->vertices[i] += offset; }
}

void generate_animation( std::vector< std::vector<mcl::Vec3f> > &verts, mcl::TriangleMesh::Ptr mesh ){

	float dt = 1.f / float(config.frame_rate);
	float curr_t = 0.f;
	float end_t = std::abs(config.end_time);

	while( curr_t < end_t ){
		mesh->apply_xform( config.xf );
		int frame_idx = verts.size();
		verts.emplace_back( std::vector<mcl::Vec3f>() );
		int n_v = mesh->vertices.size();

		for( int i=0; i<n_v; ++i ){
			verts[frame_idx].emplace_back( mesh->vertices[i] );
		}
		curr_t += dt;
	}

}














