//  BSD 3-Clause License
//
//  Copyright (c) 2019, Bailin Deng
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

#include "ALMGeometrySolver.h"
#include "Constraint.h"
#include "MeshTypes.h"
#include <iostream>
#include <random>
#include "Parameters.h"

void save_face_index(Eigen::Matrix4Xi & face_vertex_index)
{
    std::string file1;

    file1 = "./result/face_vertex_index.txt";

    std::ofstream ofs1;
    ofs1.open(file1, std::ios::out | std::ios::ate);
    if (!ofs1.is_open())
    {
        std::cout << "Cannot open: " << file1 << std::endl;
        return;
    }

    ofs1 << std::setprecision(16);
    for (size_t i = 0; i < face_vertex_index.cols(); i++)
    {
        ofs1 << face_vertex_index(0, i) << std::endl;
        ofs1 << face_vertex_index(1, i) << std::endl;
        ofs1 << face_vertex_index(2, i) << std::endl;
        ofs1 << face_vertex_index(3, i) << std::endl;
    }

    ofs1.close();
}

void save_error(const VectorX &wiremeshErr, const VectorX &wiremeshErr_after, int save_id)
{
    std::string file1, file2;

    if (save_id == 1)
    {
        file1 = "./result/edge_wiremeshErrBefore.txt";
        file2 = "./result/edge_wiremeshErrAfter.txt";
    } else if (save_id == 2){
        file1 = "./result/angle_wiremeshErrBefore.txt";
        file2 = "./result/angle_wiremeshErrAfter.txt";
    } else if (save_id == 3){
        file1 = "./result/ref_wiremeshErrBefore.txt";
        file2 = "./result/ref_wiremeshErrAfter.txt";
    }


    std::ofstream ofs1, ofs2;
    ofs1.open(file1, std::ios::out | std::ios::ate);
    ofs2.open(file2, std::ios::out | std::ios::ate);
    if (!ofs1.is_open() || (!ofs2.is_open()))
    {
        std::cout << "Cannot open: " << file1 << " or " << file2 << std::endl;
        return;
    }

    ofs1 << std::setprecision(16);
    ofs2 << std::setprecision(16);
    for (size_t i = 0; i < wiremeshErr.size(); i++)
    {
        ofs1 << wiremeshErr[i] << std::endl;
        ofs2 << wiremeshErr_after[i] << std::endl;
    }

    ofs1.close();
    ofs2.close();
}

void check_wiremesh_error(const PolyMesh &mesh, double target_edge_length,
                          double min_angle_radian, double max_angle_radian,
                          VectorX &angle_error, VectorX &edge_error) {
    VectorX edge_err(mesh.n_edges());
    VectorX angle_err(mesh.n_faces() * 4), edge_err_out(mesh.n_faces() * 4);
    angle_error = angle_err;

    for (PolyMesh::ConstEdgeIter ce_it = mesh.edges_begin();
         ce_it != mesh.edges_end(); ++ce_it) {
        edge_err(ce_it->idx()) = std::fabs(
                    mesh.calc_edge_length(*ce_it) - target_edge_length);
    }
    edge_err /= target_edge_length;

    int edge_err_id = 0;
    for (PolyMesh::ConstFaceIter cf_it = mesh.faces_begin();
         cf_it != mesh.faces_end(); ++cf_it) {
        std::vector<Vector3> vtx;
        for (PolyMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*cf_it);
             cfv_it.is_valid(); ++cfv_it) {
            vtx.push_back(to_eigen_vec3(mesh.point(*cfv_it)));
        }

        assert(static_cast<int>(vtx.size()) == 4);
        for (int i = 0; i < 4; ++i) {
            Vector3 e1 = (vtx[(i + 1) % 4] - vtx[i]).normalized(), e2 = (vtx[(i + 3)
                    % 4] - vtx[i]).normalized();
            Scalar angle = std::acos(e1.dot(e2));

            angle_error[4 * cf_it->idx() + i] = abs(angle - 0.5*M_PI);
            if (angle < min_angle_radian) {
                angle_err[4 * cf_it->idx() + i] = min_angle_radian - angle;
            } else if (angle >= max_angle_radian) {
                angle_err[4 * cf_it->idx() + i] = angle - max_angle_radian;
            } else {
                angle_err[4 * cf_it->idx() + i] = 0;
            }
        }

        for (PolyMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*cf_it);
            cfh_it.is_valid(); ++cfh_it) {
            edge_err_out(edge_err_id) = edge_err(mesh.edge_handle(*cfh_it).idx());
            edge_err_id++;
        }
    }

    angle_err *= (180.0 / M_PI);
    angle_error *= (180.0 / M_PI);

    edge_error = edge_err_out;

    std::cout << "Normalized edge length error: max " << edge_err.maxCoeff()
              << ",  average " << edge_err.mean() << std::endl;
    std::cout << "Angle error: max " << angle_err.maxCoeff() << ",  average "
              << angle_err.mean() << std::endl;
}

void check_ref_surface_distance(const PolyMesh &mesh, const TriMesh &ref_mesh, VectorX &ref_error) {
    Matrix3X ref_mesh_points, mesh_pts;
    Eigen::Matrix3Xi ref_mesh_faces;
    get_vertex_points(ref_mesh, ref_mesh_points);
    get_face_vertex_index(ref_mesh, ref_mesh_faces);
    get_vertex_points(mesh, mesh_pts);

    MatrixX3 V_ref = ref_mesh_points.transpose(), P = mesh_pts.transpose();
    Eigen::MatrixX3i F_ref = ref_mesh_faces.transpose();

    igl::AABB<MatrixX3, 3> tree;
    tree.init(V_ref, F_ref);
    VectorX sqrD;
    Eigen::VectorXi I;
    MatrixX3 C;
    tree.squared_distance(V_ref, F_ref, P, sqrD, I, C);
    Scalar edge_length = average_edge_length(mesh);
    VectorX dist_err = sqrD.array().sqrt().matrix() / edge_length;
    std::cout << "Reference surface distance (normalized by edge length): ";
    std::cout << "Max " << dist_err.maxCoeff() << ", ";
    std::cout << "Average " << dist_err.mean() << std::endl;

    ref_error = dist_err;
}


bool setup_quad_laplacian_matrix(const PolyMesh &mesh, Scalar laplacian_weight, ALMGeometrySolver<3> &solver)
{
    std::vector<Scalar> coefs;
    coefs.push_back(Scalar(2));
    coefs.push_back(Scalar(-1));
    coefs.push_back(Scalar(-1));

    PolyMesh::ConstVertexVertexIter vv_it;

    for(PolyMesh::ConstVertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++ v_it){
        int center_idx = v_it->idx();
        std::vector<int> neighbor_idx;

        for(vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++ vv_it){
            neighbor_idx.push_back(vv_it->idx());
        }

        int m = neighbor_idx.size();
        if(m > 4){
            std::cout << "Invalid valence" << std::endl;
            return false;
        }
        else if (m == 4){
            solver.add_laplacian(std::vector<int>({center_idx, neighbor_idx[0], neighbor_idx[2]}), coefs, laplacian_weight);
            solver.add_laplacian(std::vector<int>({center_idx, neighbor_idx[1], neighbor_idx[3]}), coefs, laplacian_weight);
        }
        else if (m == 3){
            if(!mesh.is_boundary(*v_it)){
                std::cout << "Not a regular quad mesh" << std::endl;
                return false;
            }

            std::vector<int> indices;
            indices.push_back(center_idx);
            for(PolyMesh::ConstVertexOHalfedgeIter cvoh_it = mesh.cvoh_iter(*v_it); cvoh_it.is_valid(); ++ cvoh_it){
                if(mesh.is_boundary(mesh.edge_handle(*cvoh_it))){
                    indices.push_back(mesh.to_vertex_handle(*cvoh_it).idx());
                }
            }

            solver.add_laplacian(indices, coefs, laplacian_weight);
        }
    }

    return true;
}

// Anderson_M = 0 means no acceleration
void optimize_mesh(const PolyMesh &mesh, const TriMesh &ref_mesh,
                   int max_iter, int Anderson_m,
                   Scalar penalty_parameter,
                   Scalar min_angle_radian,
                   Scalar max_angle_radian,
                   Scalar edge_length,
                   Scalar closeness_weight,
                   Scalar laplacian_weight,
                   const char *output_filename = NULL)
{


    Matrix3X p;
    get_vertex_points(mesh, p);

    Matrix3X ref_mesh_points;
    Eigen::Matrix3Xi ref_mesh_faces;
    get_vertex_points(ref_mesh, ref_mesh_points);
    get_face_vertex_index(ref_mesh, ref_mesh_faces);

    ALMGeometrySolver<3> solver;

    if (closeness_weight > 0) {
        solver.add_soft_constraint(
                    new ReferenceSurfceConstraint(p.cols(), closeness_weight,
                                                  ref_mesh_points, ref_mesh_faces));
    }

    for (PolyMesh::ConstFaceIter f_it = mesh.faces_begin();
         f_it != mesh.faces_end(); ++f_it) {
        std::vector<int> id_vector;
        for (PolyMesh::ConstFaceVertexIter fv_it = mesh.cfv_iter(*f_it);
             fv_it.is_valid(); ++fv_it) {
            id_vector.push_back(fv_it->idx());
        }

        assert(static_cast<int>(id_vector.size()) == 4);
        for (int i = 0; i < 4; ++i) {
            solver.add_hard_constraint(
                        new AngleConstraint<3>(id_vector[i], id_vector[(i + 1) % 4],
                        id_vector[(i + 3) % 4], 1.0,
                    min_angle_radian, max_angle_radian));
        }
    }

    for (PolyMesh::ConstEdgeIter ce_it = mesh.edges_begin();
         ce_it != mesh.edges_end(); ++ce_it) {
        PolyMesh::HalfedgeHandle heh = mesh.halfedge_handle(*ce_it, 0);
        int v1 = mesh.from_vertex_handle(heh).idx(), v2 = mesh.to_vertex_handle(heh)
                .idx();
        solver.add_hard_constraint(
                    new EdgeLengthConstraint<3>(v1, v2, 1.0, edge_length));
    }

    if (laplacian_weight > 0) {
        setup_quad_laplacian_matrix(mesh, laplacian_weight, solver);
    }

    Scalar eps_ratio = 1e-8;
    Scalar rel_residual_eps = eps_ratio * average_edge_length(mesh);
    std::cout << "Relative residual eps (normalized by edge length): " << eps_ratio << std::endl;

    if (solver.setup_ADMM(p.cols(), penalty_parameter)) {
        solver.solve_ADMM(p, rel_residual_eps, max_iter, Anderson_m);
        //solver.output_iteration_history(AA_SOLVER);

        //    if (Anderson_m > 0)
        //        solver.output_iteration_history(AA_SOLVER);
        //    else
        //        solver.output_iteration_history(SHAPE_UP_SOLVER);

        solver.save(Anderson_m);

        VectorX wiremesh_err_bef(mesh.n_vertices()), wiremesh_err_aft(mesh.n_vertices());
        VectorX angle_err_bef(mesh.n_faces()*4), angle_err_aft(mesh.n_faces()*4);
        VectorX edge_err_bef(mesh.n_faces()*4), edge_err_aft(mesh.n_faces()*4);


        std::cout << "Before optimization:" << std::endl;
        check_wiremesh_error(mesh, edge_length, min_angle_radian,
                             max_angle_radian, angle_err_bef, edge_err_bef);
        check_ref_surface_distance(mesh, ref_mesh, wiremesh_err_bef);

        PolyMesh new_mesh = mesh;
        set_vertex_points(new_mesh, solver.get_solution());
        std::cout << "After optimization:" << std::endl;
        check_wiremesh_error(new_mesh, edge_length, min_angle_radian,
                             max_angle_radian, angle_err_aft, edge_err_aft);
        check_ref_surface_distance(new_mesh, ref_mesh, wiremesh_err_aft);

        save_error(edge_err_bef, edge_err_aft, 1);
        save_error(angle_err_bef, angle_err_aft, 2);
        save_error(wiremesh_err_bef, wiremesh_err_aft, 3);

        if (output_filename) {
            if (!OpenMesh::IO::write_mesh(new_mesh, output_filename,
                                          OpenMesh::IO::Options::Default, 16)) {
                std::cerr << "Error: unable to save perturbed mesh to file "
                          << output_filename << std::endl;
            }
        }
    } else {
        std::cerr << "Error: unable to initialize solver" << std::endl;
    }
}



int main(int argc, char **argv) {
    if (argc != 5) {
        std::cout
                << "Usage:   <WireMeshOpt>  <INPUT_POLY_MESH>  <REF_TRI_MESH>  <OPTIONS_FILE>  <OUTPUT_MESH>"
                << std::endl;
        return 1;
    }

    PolyMesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
        std::cerr << "Error: unable to read input mesh from the file " << argv[1]
                  << std::endl;
        return 1;
    }

    TriMesh ref_mesh;
    if (!OpenMesh::IO::read_mesh(ref_mesh, argv[2])) {
        std::cerr << "Error: unable to read referece mesh from file " << argv[2]
                  << std::endl;
        return 1;
    }

    Scalar edge_length = average_edge_length(mesh);
    Scalar min_angle_radian = M_PI * 0.25, max_angle_radian = M_PI * 0.75;

    PolyMesh sub_mesh = subdivide_and_smooth_mesh(mesh);
    edge_length *= 0.5;
    std::cout << "target length = " << edge_length << std::endl;
    //std::cout << "diag length = " << bbox_diag_length(sub_mesh) << std::endl;

    //For reference mesh error.
//    PolyMesh ori_color_mesh = quad_subdivision(sub_mesh);
//    if (!OpenMesh::IO::write_mesh(ori_color_mesh, "/home/py/Desktop/Wiremesh_com_AA/Data/ori_ref_color_mesh.obj",
//                                  OpenMesh::IO::Options::Default, 16)) {
//        std::cerr << "Error: unable to save perturbed mesh to file "
//                  << "/home/py/Desktop/Wiremesh_com_AA/Data/ori_ref_color_mesh.obj" << std::endl;
//    }

    //  check_wiremesh_error(sub_mesh, edge_length, min_angle_radian,
    //                       max_angle_radian);
    //  check_ref_surface_distance(sub_mesh, ref_mesh);

    //  if (argc == 5) {
    //    if (!OpenMesh::IO::write_mesh(sub_mesh, argv[4],
    //                                  OpenMesh::IO::Options::Default, 16)) {
    //      std::cerr << "Error: unable to save subdvided mesh to file " << argv[4]
    //                << std::endl;
    //    }
    //  }

    Parameters param;
    if (!param.load(argv[3])) {
        std::cerr << "Error: unable to load option file " << argv[3] << std::endl;
        return 1;
    }
    if (!param.valid_parameters()) {
        std::cerr << "Invalid filter options. Aborting..." << std::endl;
        return 1;
    }
    param.output();

    Scalar closeness_weight = 1, laplacian_weight = -1;
    Scalar penalty_parameter = 1000;


    optimize_mesh(sub_mesh, ref_mesh, param.iter, param.anderson_m, penalty_parameter, min_angle_radian,
                  max_angle_radian, edge_length, closeness_weight, laplacian_weight, argv[4]);

    //For reference error color map
//    Eigen::Matrix4Xi face_vertex_idx;
//    get_face_vertex_index(sub_mesh, face_vertex_idx);
//    save_face_index(face_vertex_idx);

//    PolyMesh result_mesh;
//    if (!OpenMesh::IO::read_mesh(result_mesh, "/home/py/Desktop/Wiremesh_com_AA/Data/result_MaleTorso.obj")) {
//        std::cerr << "Error: unable to read input mesh from the file /home/py/Desktop/Wiremesh_com_AA/Data/result_MaleTorso.obj"
//                  << std::endl;
//        return 1;
//    }
//    PolyMesh angle_color_mesh = quad_subdivision(result_mesh);
//    if (!OpenMesh::IO::write_mesh(angle_color_mesh, "/home/py/Desktop/Wiremesh_com_AA/Data/result_ref_color_mesh.obj",
//                                  OpenMesh::IO::Options::Default, 16)) {
//        std::cerr << "Error: unable to save perturbed mesh to file "
//                  << "/home/py/Desktop/Wiremesh_com_AA/Data/result_ref_color_mesh.obj" << std::endl;
//    }


//    TriMesh edge_color_mesh = quad2tri_subdivision(sub_mesh);
//    if (!OpenMesh::IO::write_mesh(edge_color_mesh, "./edge_color_mesh.obj",
//                                  OpenMesh::IO::Options::Default, 16)) {
//        std::cerr << "Error: unable to save edge color mesh to file "
//                  << "./edge_color_mesh.obj" << std::endl;
//    }

//    PolyMesh angle_color_mesh = quad_subdivision(sub_mesh);
//    if (!OpenMesh::IO::write_mesh(angle_color_mesh, "./result/angle_color_mesh.obj",
//                                  OpenMesh::IO::Options::Default, 16)) {
//        std::cerr << "Error: unable to save perturbed mesh to file "
//                  << "./result/angle_color_mesh.obj" << std::endl;
//    }

    return 0;
}

