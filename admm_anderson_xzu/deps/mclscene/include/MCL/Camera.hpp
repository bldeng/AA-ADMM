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
// Based on "here be dragons" by Simon Rodriguez (http://simonrodriguez.fr/dragon)

#ifndef MCL_CAMERA_H
#define MCL_CAMERA_H

#include "XForm.hpp"

namespace mcl {

class Camera {
public:
	enum {
		MOUSE_START = 1 << 0,
		MOUSE_MOVE = 1 << 1,
		MOUSE_END = 1 << 2,
		MOUSE_LEFT = 1 << 3,
		MOUSE_RIGHT = 1 << 4,
	};

	Camera();

	// Sets the camera to the default settings
	inline void set_default();

	// Makes the current settings default
	inline void make_default();

	// Update the view matrix from the UVW vectors
	inline void update_view(float dt);

	// Update the projection matrix from the screen size
	inline void update_projection(int w, int h);

	// Perform camera movement from mouse input
	inline void mouse(int mouse, float x, float y);

	inline const XForm<float> &view() const { return m_view; } // returns view matrix
	inline const XForm<float> &projection() const { return m_proj; } // returns projection matrix
	inline Vec2f screen() const { return m_screen; } // returns screen size
	inline Vec2f &nearfar() { return m_nearfar; } // returns near/far for persp matrix
	inline Vec3f &eye() { return m_eye; } // returns eye position
	inline Vec3f &lookat() { return m_lookat; } // returns lookat point
	inline float &fov_deg() { return m_fov_deg; }
	
private:

	bool m_first_person;
	XForm<float> m_view;
	XForm<float> m_proj;
	Vec2f m_screen; // screen size
	Vec2f m_nearfar;

	bool m_left_pressed;
	bool m_right_pressed;
	Vec2f m_mouse; // mouse position (-1,1)
	Vec2f m_mouse_delta; // mouse movement (-1,1)
	float m_move_speed, d_move_speed;
	float m_rot_speed, d_rot_speed;
	float m_fov_deg, d_fov_deg;

	Vec3f m_eye, d_eye;
	Vec3f m_lookat, d_lookat;
	Vec3f m_up, d_up;
	Vec3f m_tangent, d_tangent;
};

//
//	Implementation
//
/*
inline void Camera::update_lookat( const Vec3f &lookat, float rad ){
	m_lookat = lookat;
	Vec3f viewdir = (m_lookat-m_eye);
	viewdir.normalize();
	viewdir *= rad;
	m_eye = m_lookat - viewdir;
	m_view = xform::make_lookat(m_eye, m_lookat, m_up);
}
*/

Camera::Camera() : m_screen(800,600) {

	m_nearfar = Vec2f(0.001, 1000.f);	
	m_left_pressed = false;
	m_right_pressed = false;
	m_first_person = false;

	d_eye = Vec3f(0,0,10);
	d_lookat = Vec3f(0,0,0);
	d_up = Vec3f(0,1,0);
	d_tangent = Vec3f(1,0,0);
	m_mouse = Vec2f(0,0);

	d_move_speed = 1.2f;
	d_rot_speed = 95.f;
	d_fov_deg = 45.f;

	set_default();
}

// Sets camera defaults
inline void Camera::set_default(){

	m_left_pressed = false;
	m_right_pressed = false;
	m_first_person = false;

	m_eye = d_eye;
	m_lookat = d_lookat;
	m_up = d_up;
	m_tangent = d_tangent;
	m_mouse = Vec2f(0,0);

	m_move_speed = d_move_speed;
	m_rot_speed = d_rot_speed;
	m_fov_deg = d_fov_deg;

	m_view = xform::make_lookat(m_eye, m_lookat, m_up);

	float aspect = m_screen[0]/m_screen[1];
	m_proj = xform::make_persp(m_fov_deg, aspect, m_nearfar[0], m_nearfar[1]);
}

inline void Camera::make_default(){
	d_eye = m_eye;
	d_lookat = m_lookat;
	d_up = m_up;
	d_tangent = m_tangent;

	d_move_speed = 1.2f;
	d_rot_speed = 95.f;
	d_fov_deg = 45.f;
}

inline void Camera::update_view(float dt){

	dt = std::max( dt, 1e-6f );
	Vec3f viewdir = m_lookat - m_eye;
	float scene_rad = viewdir.norm();
	viewdir /= scene_rad;

	if(m_right_pressed){
		Vec3f move_delta = (m_mouse_delta[0]*m_tangent+m_mouse_delta[1]*m_up)*m_rot_speed * dt;
		m_lookat -= move_delta;
		m_eye -= move_delta;
	}

	if(m_left_pressed){
		if( m_first_person ){
	  		m_lookat += (m_mouse_delta[0]*m_tangent+m_mouse_delta[1]*m_up)*m_rot_speed * dt;
	  		viewdir = m_lookat - m_eye;
			viewdir.normalize();
		} else {
	  		viewdir += (m_mouse_delta[0]*m_tangent + m_mouse_delta[1]*m_up)*m_rot_speed * dt;
			viewdir.normalize();
			m_eye = m_lookat - viewdir*scene_rad;
		}
	}

	m_tangent = viewdir.cross(m_up);
	m_tangent.normalize();
	m_up = m_tangent.cross(viewdir);
	m_up.normalize();
	m_view = xform::make_lookat(m_eye, m_lookat, m_up);
}


inline void Camera::mouse(int m, float x, float y){

	if( m & MOUSE_END ){
		if( m & MOUSE_LEFT ){ m_left_pressed = false; }
		if( m & MOUSE_RIGHT ){ m_right_pressed = false; }
		return;
	}

	Vec2f screen_xy(
		std::max(std::min(1.f, 2.f * x/m_screen[0] - 1.f), -1.f),
		-std::max(std::min(1.f, 2.f * y/m_screen[1] - 1.f), -1.f)
	);

	if( m & MOUSE_START ){ m_mouse = screen_xy; }
	else if( m & MOUSE_MOVE ){
		if( m & MOUSE_LEFT ){ m_left_pressed = true; }
		if( m & MOUSE_RIGHT ){ m_right_pressed = true; }
		m_mouse_delta = screen_xy - m_mouse;
		m_mouse = screen_xy;
	}

}

inline void Camera::update_projection(int w, int h){
	m_screen[0] = float(w);
	m_screen[1] = float(h);
	float aspect = m_screen[0]/m_screen[1];
	m_proj = xform::make_persp(m_fov_deg, aspect, m_nearfar[0], m_nearfar[1]);
}

} // nm mcl

#endif
