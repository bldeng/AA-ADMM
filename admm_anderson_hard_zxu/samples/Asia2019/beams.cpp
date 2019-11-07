// Original work Copyright (c) 2017 University of Minnesota
// Modified work Copyright 2019 Yue Peng
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
//
// By Matt Overby (http://www.mattoverby.net)

#include "Application.hpp"
#include "MCL/MeshIO.hpp"
#include "MCL/ShapeFactory.hpp"
#include "MCL/XForm.hpp"

using namespace mcl;

std::unique_ptr<Application> app;
//std::vector<int> left_pins; // -x
//std::vector<Eigen::Vector3d> left_points;
//std::vector<int> right_pins; // +x
//std::vector<Eigen::Vector3d> right_points;

std::vector<int> all_pins;
std::vector<Eigen::Vector3d> all_points;
std::vector<int> left_right_label;

void find_pins( const std::vector< std::shared_ptr<mcl::TetMesh> > &meshes );
void stretch_beams();

int main(int argc, char **argv){

	admm::Solver::Settings settings;
    settings.admm_iters =100;

    //Accelerate setting
    //settings.acceleration_type = admm::Solver::Settings::AccelationType::ANDERSON;
    //settings.Anderson_m = 10;

	if( settings.parse_args( argc, argv ) ){ return EXIT_SUCCESS; }

	int dim = 3;
	std::vector< std::shared_ptr<mcl::TetMesh> > meshes = {
        mcl::factory::make_tet_blocks( dim*4, dim, dim ),
        mcl::factory::make_tet_blocks( dim*4, dim, dim ),
        mcl::factory::make_tet_blocks( dim*4, dim, dim ),
	};

	std::vector< int > flags = {
        binding::NOSELFCOLLISION | binding::LINEAR,
        binding::NOSELFCOLLISION | binding::NEOHOOKEAN,
        binding::NOSELFCOLLISION | binding::STVK,
	};

	if( flags.size() != meshes.size() ){
		printf("flags.size() != meshes.size()");
		return EXIT_FAILURE;
	}

	// Center and scale the meshes
	for( int i=0; i<(int)meshes.size(); ++i ){
		Eigen::AlignedBox<float,3> aabb = meshes[i]->bounds();
		mcl::XForm<float> center = mcl::xform::make_trans<float>( -aabb.center() );
        float y = aabb.sizes()[1]; // Make each beam 1m tall.
		mcl::XForm<float> scale = mcl::xform::make_scale<float>(1.f/y, 1.f/y, 1.f/y );
		meshes[i]->apply_xform(scale*center);
	}

	// Now spread them out along y axis
	for( int i=0; i<(int)meshes.size(); ++i ){
		const mcl::XForm<float> xf_up = mcl::xform::make_trans(0.f, 1.75f, 0.f);
		const mcl::XForm<float> xf_down = mcl::xform::make_trans(0.f, -1.75f, 0.f);
		meshes[i]->flags = flags[i];
		if( i==0 ){
			meshes[i]->apply_xform( xf_up );
		} else if( i==1 ){

		} else if( i==2 ){
			meshes[i]->apply_xform( xf_down );
		}	
	}

	app = std::unique_ptr<Application>( new Application(settings) );

	// Add the dynamic meshes
    admm::Lame softRubber(10000000,0.399);
    for( int i=0; i<(int)meshes.size(); ++i ){
		app->add_dynamic_mesh( meshes[i], softRubber );
	}

	// Add a callback to stretch the beams each frame
    app->sim_cb = std::function<void ()>(stretch_beams);
    find_pins(meshes); // initial pin locations

	// Zoom out a bit
	app->renderWindow->m_camera->fov_deg() = 60.f;
    stretch_beams(); // set initial pins

	// display() calls initialize
	bool success = app->display();
	if( !success ){ return EXIT_FAILURE; }

	return EXIT_SUCCESS;
}

void stretch_beams(){

    //int n_left = left_pins.size();
    int n_pins = all_pins.size();//n_left + right_pins.size();
	std::vector<int> pins( n_pins );
	std::vector<Eigen::Vector3d> points( n_pins );

	// Compute movement:
	double dt = app->solver->settings().timestep_s;
	Eigen::Vector3d move = Eigen::Vector3d(1.f,0,0)*dt;

	// Move left and right pins a little
	for( int i=0; i<n_pins; ++i ){
        if (left_right_label[i] == 0) {
            pins[i] = all_pins[i];
            all_points[i] -= move;
            points[i] = all_points[i];
        } else {
            pins[i] = all_pins[i];
            all_points[i] += move;
            points[i] =all_points[i];
        }
	}

	app->solver->set_pins( pins, points );

}

void find_pins( const std::vector< std::shared_ptr<mcl::TetMesh> > &meshes ){

    all_pins.clear();
    all_points.clear();
    left_right_label.clear();

	// Get the left and rightmost vertices of each mesh
    int n_meshes = (int)meshes.size();
	int nv_offset = 0;
	for( int i=0; i<n_meshes; ++i ){
		Eigen::AlignedBox<float,3> aabb = meshes[i]->bounds();
		float min_x = aabb.min()[0] + 1e-2f;
		float max_x = aabb.max()[0] - 1e-2f;
		int n_verts = meshes[i]->vertices.size();
		for( int j=0; j<n_verts; ++j ){
			const Vec3f &v = meshes[i]->vertices[j];
            if( v[0] < min_x ){
                all_pins.emplace_back(j + nv_offset);
                all_points.emplace_back(v.cast<double>());
                left_right_label.emplace_back(0);
            }
			if( v[0] > max_x ){
                all_pins.emplace_back(j + nv_offset);
                all_points.emplace_back(v.cast<double>());
                left_right_label.emplace_back(1);
			}
		} // end loop verts
		nv_offset += n_verts;
	} // end loop meshes
}

