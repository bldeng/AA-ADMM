//  BSD 3-Clause License
//
//  Copyright (c) 2019, Yue Peng
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
//  * Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
        mcl::XForm<float> xf = mcl::xform::make_trans(0.25f, 2.5f, .0f) *
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

    /////////
    //std::unordered_map< hashkey::sint2, int > edge_ids;
    /*double edge_length = 0.0;
    int edge_num = 0;
    for( int id=0; id<(int)meshes.size(); ++id ){
        int n_tets = meshes[id]->tets.size();
        for( int f=0; f<n_tets; ++f )
            for( int i=0; i<4; ++i )
                for( int j=0; j<4; ++j ){
                    if( i==j ){ continue; }
                    Eigen::Vector3f edge = meshes[id]->vertices[meshes[id]->tets[f][i]]-meshes[id]->vertices[meshes[id]->tets[f][j]];
                    edge_length += (double)edge.norm();
                    edge_num++;
                }
    }
    std::cout << "average edge length = " << edge_length/edge_num << std::endl;*/

    // Add a collision floor
    float floor_y = -3.f, radius = 1.0f;
    Vector3 center(0.0, floor_y, 0.0);
    std::shared_ptr<admm::PassiveCollision> PlaneAndHalfSphere_collider =
            std::make_shared<admm::PlaneAndHalfSphere>( admm::PlaneAndHalfSphere(center, radius) );

    Vec3<float> center2(0.0, floor_y+2.5f, 0.0);
    // Add a floor renderable
    std::shared_ptr<mcl::TriangleMesh> PlaneAndHalfSphere = mcl::factory::make_planeHsphere(center2,radius,300,2,2);
    mcl::XForm<float> xf = mcl::xform::make_trans(0.f,-2.5f,0.f);
    PlaneAndHalfSphere->apply_xform(xf);

    app.add_obstacle( PlaneAndHalfSphere_collider, PlaneAndHalfSphere );
    set_collision(meshes);

    app.solver->set_collisions( all_idx, all_points );

    app.renderWindow->m_camera->eye() = mcl::Vec3f(0,10,20);
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




