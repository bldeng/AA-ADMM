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

#ifndef MCL_CONTROLLER_H
#define MCL_CONTROLLER_H 1

#include "Camera.hpp"

namespace mcl {

// Interface for glfwSetWindowUserPointer
class Controller {
public:
	Controller() : screen_dt(0) {}

	virtual void init_glfw(GLFWwindow* window);

	// Called each frame to update time elapsed
	inline void frame( float dt ){ screen_dt = dt; }

	// If camera is left null, its functions aren't called
	inline void set_camera( std::shared_ptr<mcl::Camera> c ){ camera = c; }

	// glfw3 callback functions
	virtual void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
	virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height);
	virtual void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
	virtual void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos);
	virtual void scroll_callback(GLFWwindow* window, double x, double y);

protected:
	std::shared_ptr<mcl::Camera> camera;
	float screen_dt;
};


//
//	Implementation
//


namespace internal {
	static void key_cb(GLFWwindow* window, int key, int scancode, int action, int mods){
		void *data = glfwGetWindowUserPointer(window);  
		Controller *w = static_cast<Controller*>(data);
		w->key_callback(window,key,scancode,action,mods);
	}
	static void framebuffer_size_cb(GLFWwindow* window, int width, int height){
		void *data = glfwGetWindowUserPointer(window);  
		Controller *w = static_cast<Controller*>(data);
		w->framebuffer_size_callback(window,width,height);
	}
	static void mouse_button_cb(GLFWwindow* window, int button, int action, int mods){
		void *data = glfwGetWindowUserPointer(window);  
		Controller *w = static_cast<Controller*>(data);
		w->mouse_button_callback(window,button,action,mods);
	}
	static void cursor_pos_cb(GLFWwindow* window, double xpos, double ypos){
		void *data = glfwGetWindowUserPointer(window);  
		Controller *w = static_cast<Controller*>(data);
		w->cursor_pos_callback(window,xpos,ypos);
	}
	static void scroll_cb(GLFWwindow* window, double x, double y){
		void *data = glfwGetWindowUserPointer(window);  
		Controller *w = static_cast<Controller*>(data);
		w->scroll_callback(window,x,y);
	}
} // ns internal


void Controller::init_glfw(GLFWwindow* window){
	if( !window ){ return; }
	glfwSetWindowUserPointer(window, this);
	glfwSetKeyCallback(window, &internal::key_cb);
	glfwSetFramebufferSizeCallback(window, &internal::framebuffer_size_cb);
	glfwSetMouseButtonCallback(window,internal::mouse_button_cb);
	glfwSetCursorPosCallback(window,internal::cursor_pos_cb);
	glfwSetScrollCallback(window, internal::scroll_cb);
}


inline void Controller::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
	MCL_UNUSED(mods);
	MCL_UNUSED(scancode);
	if( action != GLFW_PRESS ){ return; }
	bool has_cam = bool(camera);
	switch ( key ){
		case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GL_TRUE); break;
		case GLFW_KEY_R:{
			if( has_cam ){
				camera->set_default(); 
				camera->update_view( screen_dt );
			}
		} break;
		default: break;
	}
}


inline void Controller::framebuffer_size_callback(GLFWwindow* window, int width, int height){
	MCL_UNUSED(window);
	glViewport(0,0,width,height);
	camera->update_projection(width,height);
}


inline void Controller::mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
	MCL_UNUSED(mods);
	if( !camera ){ return; }
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS) {
			camera->mouse(Camera::MOUSE_START | Camera::MOUSE_LEFT,float(x),float(y));
		}
		else if(action == GLFW_RELEASE) {
			camera->mouse(Camera::MOUSE_END | Camera::MOUSE_LEFT,0,0);
		}
	} else if( button == GLFW_MOUSE_BUTTON_RIGHT){
		if(action == GLFW_PRESS) {
			camera->mouse(Camera::MOUSE_START | Camera::MOUSE_RIGHT,float(x),float(y));
		}
		else if(action == GLFW_RELEASE) {
			camera->mouse(Camera::MOUSE_END | Camera::MOUSE_RIGHT,0,0);
		}
	}
}


inline void Controller::scroll_callback(GLFWwindow* window, double x, double y){
	MCL_UNUSED(window);
	MCL_UNUSED(x);
	if( !camera ){ return; }
	float new_fov = camera->fov_deg()-5.f*y; // zoomies
	camera->fov_deg() = std::min( std::max( 1e-8f, new_fov ), 179.f );
	int w=256, h=256;
	glfwGetFramebufferSize(window, &w, &h);
	camera->update_projection(w,h);
}


inline void Controller::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos){
	if( !camera ){ return; }
	bool left_down = glfwGetMouseButton(window,GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
	bool right_down = glfwGetMouseButton(window,GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
	if( left_down ){ camera->mouse( Camera::MOUSE_MOVE | Camera::MOUSE_LEFT, float(xpos), float(ypos) ); }
	if( right_down ){ camera->mouse( Camera::MOUSE_MOVE | Camera::MOUSE_RIGHT, float(xpos), float(ypos) ); }
	if( left_down || right_down ){ camera->update_view( screen_dt ); }
}

} // ns mcl

#endif
