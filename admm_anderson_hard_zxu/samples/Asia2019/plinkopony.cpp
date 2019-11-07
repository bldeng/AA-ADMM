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
#include "MCL/MeshIO.hpp"

using namespace mcl;
std::vector<int> all_idx;
std::vector<Eigen::Vector3d> all_points;
void set_collision( const std::vector< std::shared_ptr<mcl::TetMesh> > &meshes );

int main(int argc, char **argv){

    std::vector<mcl::TetMesh::Ptr> meshes = {
        mcl::TetMesh::create(),
    };

    std::stringstream file;
    file << ADMMELASTIC_ROOT_DIR << "/samples/data/horse759";//box768
    for( int i=0; i<(int)meshes.size(); ++i ){
        mcl::meshio::load_elenode( meshes[i].get(), file.str() );
        meshes[i]->flags |= binding::LINEAR;
        mcl::XForm<float> xf = mcl::xform::make_trans(0.25f, 5.0f, .0f) *
                mcl::xform::make_scale(13.0f,13.0f,13.0f);
        meshes[i]->apply_xform(xf);
    }

    admm::Solver::Settings settings;
    settings.admm_iters = 13;
    if( settings.parse_args( argc, argv ) ){ return EXIT_SUCCESS; }
    Application app(settings);

    for( int i=0; i<(int)meshes.size(); ++i ){
        app.add_dynamic_mesh( meshes[i], admm::Lame::rubber() );
    }

    // Add a collision cylinder
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 5; i++)
        {
            float x_move = i*1.5-3.0;
            float y_move = j*3.0-3.0;
            Vector3 center(x_move, y_move, 0.0);
            float radius = 0.4;
            std::shared_ptr<admm::PassiveCollision> cylinder_collider =
                std::make_shared<admm::Cylinder>( admm::Cylinder(center, radius) );

            // Add a floor renderable
            std::shared_ptr<mcl::TriangleMesh> cyl = mcl::factory::make_cyl(30,3,radius);
            mcl::XForm<float> xf = mcl::xform::make_trans(x_move, y_move, 0.f) *
                //mcl::xform::make_rot(-90.f,mcl::Vec3f(1,0,0)) *
                mcl::xform::make_scale(1.f,1.f,1.f);
            cyl->apply_xform(xf);

            app.add_obstacle( cylinder_collider, cyl );
        }
    }

    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < 4; i++)
        {
            float x_move = i*1.5-2.25;
            float y_move = j*3.0-1.5;
            Vector3 center(x_move, y_move, 0.0);
            float radius = 0.4;
            std::shared_ptr<admm::PassiveCollision> cylinder_collider =
                std::make_shared<admm::Cylinder>( admm::Cylinder(center, radius) );

            // Add a floor renderable
            std::shared_ptr<mcl::TriangleMesh> cyl = mcl::factory::make_cyl(30,3,radius);
            mcl::XForm<float> xf = mcl::xform::make_trans(x_move, y_move, 0.f) *
                //mcl::xform::make_rot(-90.f,mcl::Vec3f(1,0,0)) *
                mcl::xform::make_scale(1.f,1.f,1.f);
            cyl->apply_xform(xf);

            app.add_obstacle( cylinder_collider, cyl );
        }
    }

    // Add a collision floor
    float floor_y = -6.5f;
    Vector3 center(0.0, floor_y, 0.0);
    Vector3 normal(0.5, std::sqrt(3.0)/2.0, 0.0);
    std::shared_ptr<admm::PassiveCollision> floor_collider =
        std::make_shared<admm::SlideFloor>( admm::SlideFloor(center, normal) );

    // Add a floor renderable
    std::shared_ptr<mcl::TriangleMesh> floor = mcl::factory::make_plane(2,2);
    mcl::XForm<float> xf = mcl::xform::make_trans(0.f,floor_y,0.f) *
        mcl::xform::make_rot(-90.f,mcl::Vec3f(1,0,0)) *
            mcl::xform::make_rot(30.f,mcl::Vec3f(0,1,0)) *
        mcl::xform::make_scale(8.f,4.f,8.f);
    floor->apply_xform(xf);
    app.add_obstacle( floor_collider, floor );

    set_collision(meshes);

    app.solver->set_collisions( all_idx, all_points );

    app.renderWindow->m_camera->eye() = mcl::Vec3f(0,0,20);
    bool success = app.display();
    if( !success ){ return EXIT_FAILURE; }

    return EXIT_SUCCESS;
}

void set_collision( const std::vector< std::shared_ptr<mcl::TetMesh> > &meshes ){

    all_idx.clear();
    all_points.clear();

    // Get the left and rightmost vertices of each mesh
    int n_meshes = meshes.size();
    int nv_offset = 0;
    for( int i=0; i<n_meshes; ++i ){
        int n_verts = meshes[i]->vertices.size();
        for( int j=0; j<n_verts; ++j ){
            const Vec3f &v = meshes[i]->vertices[j];
            all_idx.emplace_back(j + nv_offset);
            all_points.emplace_back(v.cast<double>());
        } // end loop verts
        nv_offset += n_verts;
    } // end loop meshes

}




