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

#ifndef MCL_RENDERWINDOW_H
#define MCL_RENDERWINDOW_H 1

#include "Shader.hpp"
#include "RenderMesh.hpp"
#include "Controller.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

namespace mcl {

namespace shaders {
const std::string mcl_vert_src =
#include "shaders/shader.vert"
;
const std::string mcl_frag_src =
#include "shaders/shader.frag"
;
}

// 0 = point (TODO add more later), and a nice class
// that manages things like sampling and whatnot
class Light {
public:
	short type;
	Vec3f color;
	Vec3f pos;
	Vec3f dir;
	float rad;

	static inline Light create_point( const Vec3f &p, const Vec3f &c=Vec3f(1,1,1) ){
		return Light( 0, c, p, Vec3f(0,0,0), 0 );
	}

	Light( short type_, const Vec3f &color_, const Vec3f &pos_, const Vec3f &dir_, float rad_ ) :
		type(type_), color(color_), pos(pos_), dir(dir_), rad(rad_) {}
};


class RenderWindow {
public:
	typedef Eigen::AlignedBox<float,3> AABB;

	RenderWindow();
	~RenderWindow();

	std::shared_ptr<mcl::Camera> m_camera; // manages xforms

	// Initialize OpenGL and the render window
	inline GLFWwindow* init();

	// Use a custom controller
	inline void set_controller( std::shared_ptr<mcl::Controller> c );

	inline bool is_open(){ return !glfwWindowShouldClose(m_window); }

	// Draw all meshes:
	// - Sets lighting
	// - Draws meshes
	// - Calls glfwSwapBuffers
	inline void draw();

	// Add a mesh to be drawn
	inline void add_mesh( std::shared_ptr<mcl::RenderMesh> mesh );

	// Returns bounding box for scene objects (not cam or lights)
	inline AABB bounds() const;

	// Sets the camera to a nice location based on the AABB
	inline void nice_camera_location();

	// Sets up three-point lighting for a scene with
	// camera position (eye) and scene center (c).
	inline void make_3pt_lighting( const Vec3f &eye, const AABB &aabb );

	// Save current frame as a png file
	inline void save_screenshot( const std::string &filename );

protected:

	std::vector< std::shared_ptr<mcl::RenderMesh> > m_meshes;
	std::shared_ptr<mcl::Controller> m_controller;

	std::vector< mcl::Light > m_lights;
	GLFWwindow* m_window;
	std::shared_ptr<Shader> m_shader;
	float screen_dt;
	float screen_dt_old;

	// the scene aabb is used for setting camera position,
	// but also a general "reset the scene" flag on draw
	AABB m_aabb;

	// Pixels is a buffer for saving screens to a file.
	// It's stored to avoid allocation/deallocation overhead
	// when saving sequential frames. temp_pixels is used to store
	// the vertically-flipped data as read from gl.
	std::vector<unsigned char> pixels;
	std::vector<unsigned char> temp_pixels;
};


//
//	Implementation
//


namespace internal {
	static void error_cb(int error, const char* d){ fprintf(stderr, "Error (%d): %s\n", error, d); }
} // ns internal


RenderWindow::RenderWindow() : m_window(nullptr), screen_dt(0), screen_dt_old(0) {
	m_camera = std::make_shared<mcl::Camera>();
	m_shader = std::make_shared<mcl::Shader>();
	m_controller = std::make_shared<mcl::Controller>();
	m_controller->set_camera( m_camera );
}

RenderWindow::~RenderWindow(){
	if( m_window ){
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
		glUseProgram(0); 
	}
	glfwTerminate();
}

inline void RenderWindow::set_controller( std::shared_ptr<mcl::Controller> c ){
	if( m_window ){
		std::cerr << "**RenderWindow Error: Must set controller before init()" << std::endl;
		return;
	}
	m_controller = c;
	m_controller->set_camera( m_camera );
}

inline GLFWwindow* RenderWindow::init(){

	// Set up window
	glfwSetErrorCallback(&internal::error_cb);

	// Initialize the window
	if( !glfwInit() ){ return nullptr; }
	glfwWindowHint(GLFW_SAMPLES, 4); // anti aliasing

	// Ask for OpenGL 3.30 core
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Create the glfw window
	Vec2f init_screen = m_camera->screen();
	m_window = glfwCreateWindow(int(init_screen[0]), int(init_screen[1]), "Viewer", NULL, NULL);
	if( !m_window ){ glfwTerminate(); return nullptr; }

	// Make current
	glfwMakeContextCurrent(m_window);
	glfwSwapInterval(1); // 0=none, 1=60FPS, 2=30FPS

	// Initialize glew AFTER the context creation and before loading the shader.
	// Note we need to use experimental because we're using a modern version of opengl.
	#ifdef MCL_USE_GLEW
	glewExperimental = GL_TRUE;
	glewInit();
	#endif

	// Initialize OpenGL
//	glEnable(GL_CULL_FACE);
//	glCullFace(GL_BACK);
//	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.f,1.f,1.f,1.f);

	// Init stuff
	m_controller->init_glfw(m_window);
	m_shader->init_from_strings( shaders::mcl_vert_src, shaders::mcl_frag_src );
	m_camera->make_default(); // set current settings as default

	// Initialize the camera
	// IMPORTANT: Only call after gl context has been created
	m_controller->framebuffer_size_callback(m_window,int(init_screen[0]),int(init_screen[1]));
	m_controller->key_callback(m_window, 0, 0, GLFW_PRESS, 0); // sets view

	// All done
	screen_dt_old = glfwGetTime();
	return m_window;
}


inline Eigen::AlignedBox<float,3> RenderWindow::bounds() const {
	AABB box;
	int n_meshes = m_meshes.size();
	for( int i=0; i<n_meshes; ++i ){
		box.extend( m_meshes[i]->bounds() );
	}
	return box;
}

inline void RenderWindow::nice_camera_location(){
	AABB box = bounds();
	m_camera->lookat() = box.center();

	Vec3f diag = (box.max() - box.min() )*0.5;
	diag[2] = -diag[2]*4.f; // Move out z a bit
	m_camera->eye() = box.min() + diag;
	m_camera->nearfar()[1] = ( box.max()[2] - m_camera->eye()[2] )*10.f;
}

inline void RenderWindow::draw(){

	m_shader->enable();

	// If the scene aabb is empty, reset:
	// - camera position
	// - lighting
	if( m_aabb.isEmpty() ){

		m_aabb = bounds();
		m_camera->update_view(0.f);

		// Create lights if we haven't already
		if( m_lights.size()==0 ){ make_3pt_lighting( m_camera->eye(), m_aabb ); }

		// Set lighting uniforms
		int n_l = std::min( m_lights.size(), size_t(8) );
		glUniform1i( m_shader->uniform( "num_lights" ), n_l );
		for( int i=0; i<n_l; ++i ){
			std::string light_str = "lights[" + std::to_string(i) + "].";
			const Light &light = m_lights[i];
			glUniform1i( m_shader->uniform( light_str+"type" ), light.type );
			glUniform3f( m_shader->uniform( light_str+"position" ), light.pos[0], light.pos[1], light.pos[2] );
			glUniform3f( m_shader->uniform( light_str+"color" ), light.color[0], light.color[1], light.color[2] );
		}

	} // end aabb is empty

	// Update screen dt
	float t = glfwGetTime();
	screen_dt = t - screen_dt_old;
	screen_dt_old = t;
	m_controller->frame( screen_dt );

	// Clear screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Get camera matrices
	const Eigen::Matrix4f view = m_camera->view().matrix();
	const Eigen::Matrix4f projection = m_camera->projection().matrix();
	const Eigen::Matrix4f inv_trans_view = view.inverse().transpose();
	Eigen::Matrix4f iden = Eigen::Matrix4f::Identity();
	glUniformMatrix4fv( m_shader->uniform("model"), 1, GL_FALSE, iden.data() );
	glUniformMatrix4fv( m_shader->uniform("view"), 1, GL_FALSE, view.data() );
	glUniformMatrix4fv( m_shader->uniform("projection"), 1, GL_FALSE, projection.data() );
	glUniformMatrix4fv( m_shader->uniform("invmv"), 1, GL_FALSE, inv_trans_view.data() );

	// Draw meshes
	int n_m = m_meshes.size();
	for( int i=0; i<n_m; ++i ){

		const Eigen::Matrix4f &model = m_meshes[i]->get_model().matrix();
		bool has_model = ( model.matrix() - iden ).squaredNorm() > 1e-10;
		if( has_model ){ throw std::runtime_error("RenderWindow TODO: Model matrix"); } // ... TODO

		// Set material
		const material::Phong &phong = m_meshes[i]->phong;
		glUniform3f( m_shader->uniform("material.amb"), phong.amb[0], phong.amb[1], phong.amb[2] );
		glUniform3f( m_shader->uniform("material.diff"), phong.diff[0], phong.diff[1], phong.diff[2] );
		glUniform3f( m_shader->uniform("material.spec"), phong.spec[0], phong.spec[1], phong.spec[2] );
		glUniform1f( m_shader->uniform("material.shini"), phong.shini );

		// Draw the mesh
		m_meshes[i]->draw();
	}

	glfwSwapBuffers(m_window);
}


inline void RenderWindow::make_3pt_lighting( const Vec3f &eye, const AABB &aabb ){

	const Vec3f center = aabb.center();
//	const float rad = aabb.sizes().maxCoeff()*0.5; // for falloff
	const float top = aabb.max()[1] - center[1];
	const Vec3f up(0,1,0);

	Vec3f w = eye-center;
	float dist = std::max( 0.1f, w.norm() );
	w /= dist;
	Vec3f u = up.cross(w); // right
	u.normalize();

	Vec3f key = center + w*dist + up*eye[1] - u*dist; // left of camera
	Vec3f fill = center + w*(dist*0.5f) + up*top + u*dist; // right of camera
	Vec3f back = center - w*dist + up*top + u*(dist*0.5f); // opposite of key

	m_lights.clear();
	m_lights.emplace_back( Light::create_point( key, Vec3f(1.f,1.f,1.f) ) );
	m_lights.emplace_back( Light::create_point( fill, Vec3f(0.8,0.8,0.8) ) );
	m_lights.emplace_back( Light::create_point( back, Vec3f(0.6,0.6,0.6) ) );
}


inline void RenderWindow::save_screenshot( const std::string &filename ){

	int w, h;
	glfwGetFramebufferSize(m_window, &w, &h);
	int row_stride = w*sizeof(unsigned char)*3;

	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	if( (int)pixels.size() != 3*w*h || (int)temp_pixels.size() != 3*w*h ){
		pixels.resize(3*w*h);
		temp_pixels.resize(3*w*h);
	}
	glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE, &temp_pixels[0]);

	// loop and swap rows (invert vertical)
	for(int i=0; i < h; ++i){
		std::memcpy( &pixels[i*w*3], &temp_pixels[(h-i-1)*w*3], row_stride );
	}

	int success = stbi_write_png(filename.c_str(), w, h, 3, &pixels[0], row_stride);
	if( !success ){ std::cerr << "**stbi_write_png error" << std::endl; }
}


inline void RenderWindow::add_mesh( std::shared_ptr<mcl::RenderMesh> mesh ){
	m_meshes.emplace_back(mesh);
	m_aabb.setEmpty(); // resized on draw()
}

} // ns mcl

#endif
