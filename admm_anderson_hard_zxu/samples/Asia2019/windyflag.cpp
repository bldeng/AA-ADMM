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

#include "Application.hpp"
#include <random>
#include "MCL/XForm.hpp"
#include "MCL/ShapeFactory.hpp"

using namespace admm;
using namespace mcl;

void get_pins( mcl::TriangleMesh::Ptr mesh, std::vector<int> &pin_ids, int idx_offset ){

    Eigen::AlignedBox<float,3> aabb = mesh->bounds();

    int up_idx = -1; // min x
    int down_idx = -1; // max x
    float min_y = aabb.min()[0] + 1e-3f;
    float curr_max_y = -99999.f;
    float curr_min_y = -curr_max_y;

    for( size_t i=0; i<mesh->vertices.size(); ++i ){
        float x = mesh->vertices[i][0];
        if( x > min_y ){ continue; }

        float y = mesh->vertices[i][1];
        if( y < curr_min_y ){
            up_idx = static_cast<int>(i);
            curr_min_y = y;
        }
        else if( y > curr_max_y ){
            down_idx = static_cast<int>(i);
            curr_max_y = y;
        }
    }

    if( up_idx < 0 || down_idx < 0 ){
        throw std::runtime_error("Failed to find pin locations");
    }

    pin_ids.emplace_back( up_idx + idx_offset );
    pin_ids.emplace_back( down_idx + idx_offset );

}

int main(int argc, char *argv[]){

    Application app;
    std::shared_ptr<WindForce> wind;
    Eigen::Vector3d orig_wind(10,0,2);

    admm::Solver::Settings settings;
    settings.admm_iters = 100;
    settings.penalty = 1.0;

    if( settings.parse_args( argc, argv ) ){ return EXIT_SUCCESS; }

    //load the mesh
	std::stringstream conf_ss;
    conf_ss << ADMMELASTIC_ROOT_DIR << "/samples/data/cloth.obj";

    mcl::TriangleMesh::Ptr mesh = mcl::TriangleMesh::create();
    mcl::meshio::load_obj(mesh.get(), conf_ss.str());
    mesh->flags = binding::NOSELFCOLLISION | binding::LINEAR;

    // Add meshes to the system
    admm::Lame very_soft_rubber(50,0.1);
    very_soft_rubber.limit_min = 0.95;
    very_soft_rubber.limit_max = 1.05;
    app.add_dynamic_mesh(mesh, very_soft_rubber);

    // Pin corners
    std::vector<int> pins;
    int pin_idx_offset = 0;
    get_pins( mesh, pins, pin_idx_offset );
    pin_idx_offset += mesh->vertices.size();

    std::vector<Eigen::Vector3d> points;
    for (size_t i = 0; i<pins.size(); i++)
    {
        const Vec3f &v = mesh->vertices[static_cast<std::vector<int>::size_type>(pins[i])];
        points.emplace_back(v.cast<double>());
    }
    app.solver->set_pins( pins );

    //
    //	Add wind manually so we can adjust the intensity with a button press
    //	Most of this is just copy-paste from SimContext
    //
    std::vector<int> faces;
    int total_sys_verts = 0;

    // Loop over all dynamic meshes in the scene, and create a vector of all faces.
    // This vector is used to create the wind force.
    std::vector<Application::DynamicMesh>::iterator o_it = app.dynamic_meshes.begin();
    for( ; o_it != app.dynamic_meshes.end(); ++o_it ){
        std::shared_ptr<mcl::TriangleMesh> mesh = o_it->trimesh;
        for( size_t f=0; f<mesh->faces.size(); ++f ){
            faces.push_back( mesh->faces[f][0]+total_sys_verts );
            faces.push_back( mesh->faces[f][1]+total_sys_verts );
            faces.push_back( mesh->faces[f][2]+total_sys_verts );
        }
        total_sys_verts += mesh->vertices.size();

    } // end loop dyanmic meshes

    //set wind force
    wind = std::shared_ptr<admm::WindForce>( new admm::WindForce(faces) );
    wind->direction = orig_wind * 2.5;
    app.solver->ext_forces.push_back(wind);

    // Try to init the solver
    if( !app.solver->initialize(settings) ){ return EXIT_FAILURE; }

    // Create opengl context
    GLFWwindow* window = app.renderWindow->init();
    if( !window ){ return EXIT_FAILURE; }

    // Add render meshes
    //load the mesh
    std::stringstream conf_ss2;
    conf_ss2 << ADMMELASTIC_ROOT_DIR << "/samples/data/pole.obj";

    mcl::TriangleMesh::Ptr mesh2 = mcl::TriangleMesh::create();
    mcl::meshio::load_obj(mesh2.get(), conf_ss2.str());
    app.add_static_mesh(mesh2);

    for( size_t i=0; i<app.dynamic_meshes.size(); ++i ){
        app.renderWindow->add_mesh( app.dynamic_meshes[i].surface );
    }
    for( size_t i=0; i<app.static_meshes.size(); ++i ){
        app.renderWindow->add_mesh( app.static_meshes[i].surface );
    }

    // Game loop
    while( app.renderWindow->is_open() ){

        //
        //	Update
        //
        if( app.controller->sim_running ){
            app.solver->step();

            for( size_t i=0; i<app.dynamic_meshes.size(); ++i ){
                app.dynamic_meshes[i].update( app.solver.get() );
            }

        } // end run continuously
        else if( app.controller->sim_dostep ){
            app.controller->sim_dostep = false;
            app.solver->step();
            for( size_t i=0; i<app.dynamic_meshes.size(); ++i ){
                app.dynamic_meshes[i].update( app.solver.get() );
            }
        } // end do one step

        //
        //	Render
        //
        app.renderWindow->draw();
        glfwPollEvents();

    } // end game loop

    return EXIT_SUCCESS;
}


