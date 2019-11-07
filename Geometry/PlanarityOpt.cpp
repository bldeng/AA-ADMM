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
#include "MeshTypes.h"
#include "TriMeshAABB.h"
#include "Constraint.h"
#include <igl/AABB.h>
#include <iostream>
#include "Parameters.h"

void save_error(const VectorX &planarityErr, const VectorX &planarityErr_after)
{
    std::string file1, file2;

    file1 = "./result/planarityErrBefore.txt";
    file2 = "./result/planatityErrAfter.txt";

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
    for (size_t i = 0; i < planarityErr.size(); i++)
    {
        ofs1 << planarityErr[i] << std::endl;
        ofs2 << planarityErr_after[i] << std::endl;
    }

    ofs1.close();
    ofs2.close();
}

void check_planarity_error(const PolyMesh &mesh, VectorX &planarityErr) {
  VectorX planarity_err(mesh.n_faces()), diag_err(mesh.n_faces());
  planarity_err.setZero();
  diag_err.setZero();

  for (PolyMesh::ConstFaceIter cf_it = mesh.faces_begin();
      cf_it != mesh.faces_end(); ++cf_it) {
    Matrix3X points(3, mesh.valence(*cf_it));
    int i = 0;
    for (PolyMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      points.col(i++) = to_eigen_vec3(mesh.point(*cfv_it));
    }

    if (static_cast<int>(points.cols()) == 4) {
      Vector3 d1 = points.col(2) - points.col(0), d2 = points.col(3)
          - points.col(1);
      Vector3 c1 = (points.col(2) + points.col(0)) * 0.5, c2 = (points.col(3)
          + points.col(1)) * 0.5;
      diag_err(cf_it->idx()) = std::fabs(
          d1.cross(d2).normalized().dot(c1 - c2));
    }

    points.colwise() -= points.rowwise().mean().eval();
    Vector3 N = points.jacobiSvd(Eigen::ComputeFullU).matrixU().col(2);
    planarity_err(cf_it->idx()) = (N.transpose() * points).array().abs()
        .maxCoeff();
  }

  Scalar edge_length = average_edge_length(mesh);
  diag_err /= edge_length;
  planarity_err /= edge_length;

  std::cout << "Diagonal error (normalized by edge length): max "
            << diag_err.maxCoeff() << ", average " << diag_err.mean()
            << std::endl;
  std::cout << "Planarity error (normalized by edge length): max "
            << planarity_err.maxCoeff() << ", average " << planarity_err.mean()
            << std::endl;

  planarityErr = planarity_err;
}

void check_ref_surface_distance(const PolyMesh &mesh, const TriMesh &ref_mesh) {
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
}

// Anderson_M = 0 means no acceleration
void optimize_mesh(const PolyMesh &mesh, const TriMesh &ref_mesh, int max_iter,
                   int Anderson_m, Scalar penalty_parameter,
                   Scalar closeness_weight, Scalar laplacian_weight,
                   Scalar relative_laplacian_weight,
                   const char *output_filename = NULL) {
  Matrix3X ref_mesh_points;
  Eigen::Matrix3Xi ref_mesh_faces;
  get_vertex_points(ref_mesh, ref_mesh_points);
  get_face_vertex_index(ref_mesh, ref_mesh_faces);

  Matrix3X p;
  get_vertex_points(mesh, p);

  ALMGeometrySolver<3> solver;
  std::shared_ptr<TriMeshAABB> aabb = std::make_shared<TriMeshAABB>(ref_mesh);

  if (closeness_weight > 0) {
    int n_vtx = p.cols();
    for (int i = 0; i < n_vtx; ++i) {
      solver.add_soft_constraint(
          new PointToRefSurfaceConstraint(i, closeness_weight, aabb));
    }
    //solver.add_constraint(new ReferenceSurfceConstraint(p.cols(), closeness_weight, ref_mesh_points, ref_mesh_faces));
  }

  for (PolyMesh::ConstVertexIter v_it = mesh.vertices_begin();
      v_it != mesh.vertices_end(); ++v_it) {
    if (laplacian_weight <= 0 && relative_laplacian_weight <= 0) {
      continue;
    }

    if (!mesh.is_boundary(*v_it)) {

      std::vector<int> vhs;
      vhs.push_back(v_it->idx());

      for (PolyMesh::ConstVertexVertexIter cvv_it = mesh.cvv_iter(*v_it);
          cvv_it.is_valid(); ++cvv_it) {
        vhs.push_back(cvv_it->idx());
      }

      if (static_cast<int>(vhs.size()) == 5) {
        if (relative_laplacian_weight > 0) {
          solver.add_relative_uniform_laplacian(std::vector<int>( { vhs[0],
              vhs[1], vhs[3] }),
                                                relative_laplacian_weight, p);
          solver.add_relative_uniform_laplacian(std::vector<int>( { vhs[0],
              vhs[2], vhs[4] }),
                                                relative_laplacian_weight, p);
        }

        if (laplacian_weight > 0) {
          solver.add_uniform_laplacian(std::vector<int>( { vhs[0], vhs[1],
              vhs[3] }),
                                       laplacian_weight);
          solver.add_uniform_laplacian(std::vector<int>( { vhs[0], vhs[2],
              vhs[4] }),
                                       laplacian_weight);
        }
      } else {
        if (relative_laplacian_weight > 0) {
          solver.add_relative_uniform_laplacian(vhs, relative_laplacian_weight,
                                                p);
        }

        if (laplacian_weight > 0) {
          solver.add_uniform_laplacian(vhs, laplacian_weight);
        }
      }
    } else {
      std::vector<int> vhs;
      vhs.push_back(v_it->idx());

      std::vector<int> fhs;
      for (PolyMesh::ConstVertexOHalfedgeIter cvoh_it = mesh.cvoh_iter(*v_it);
          cvoh_it.is_valid(); ++cvoh_it) {
        if (mesh.is_boundary(mesh.edge_handle(*cvoh_it))) {
          PolyMesh::HalfedgeHandle heh = *cvoh_it;
          vhs.push_back(mesh.to_vertex_handle(heh).idx());

          if (mesh.is_boundary(heh)) {
            heh = mesh.opposite_halfedge_handle(heh);
          }

          fhs.push_back(mesh.face_handle(heh).idx());
        }
      }

      if (static_cast<int>(fhs.size()) == 2 && fhs[0] != fhs[1]) {
        if (relative_laplacian_weight > 0) {
          solver.add_relative_uniform_laplacian(vhs, relative_laplacian_weight,
                                                p);
        }

        if (laplacian_weight > 0) {
          solver.add_uniform_laplacian(vhs, laplacian_weight);
        }
      }
    }
  }

  for (PolyMesh::ConstFaceIter f_it = mesh.faces_begin();
      f_it != mesh.faces_end(); ++f_it) {
    std::vector<int> id_vector;
    for (PolyMesh::ConstFaceVertexIter fv_it = mesh.cfv_iter(*f_it);
        fv_it.is_valid(); ++fv_it) {
      id_vector.push_back(fv_it->idx());
    }

    if (static_cast<int>(id_vector.size()) > 3) {
      solver.add_hard_constraint(new PlaneConstraint(id_vector, 1.0));
    }
  }

  Scalar eps_ratio = 1e-8;
  Scalar rel_residual_eps = eps_ratio * average_edge_length(mesh);
  std::cout << "Relative residual eps (normalized by edge length): "
            << eps_ratio << std::endl;

  if (solver.setup_ADMM(p.cols(), penalty_parameter)) {
    solver.solve_ADMM(p, rel_residual_eps, max_iter, Anderson_m);

    //    if (Anderson_m > 0)
    //        solver.output_iteration_history(AA_SOLVER);
    //    else
    //        solver.output_iteration_history(SHAPE_UP_SOLVER);

        solver.save(Anderson_m);

    std::cout << "Before optimization:" << std::endl;
    VectorX planarity_error(mesh.n_faces()), planarity_error_after(mesh.n_faces());
    check_planarity_error(mesh, planarity_error);
    check_ref_surface_distance(mesh, ref_mesh);

    //solver.output_iteration_history(AA_SOLVER);
    PolyMesh new_mesh = mesh;
    set_vertex_points(new_mesh, solver.get_solution());
    std::cout << "After optimization:" << std::endl;
    check_planarity_error(new_mesh, planarity_error_after);
    check_ref_surface_distance(new_mesh, ref_mesh);

    save_error(planarity_error, planarity_error_after);

    if (output_filename) {
      if (!OpenMesh::IO::write_mesh(new_mesh, output_filename,
                                    OpenMesh::IO::Options::Default, 16)) {
        std::cerr << "Error: unable to save result mesh to file "
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
        << "Usage:   <PlanarityOpt> <INPUT_MESH> <REFERENCE_MESH> <OPTION_FILES> <OUTPUT_MESH>"
        << std::endl;
    return 1;
  }

  PolyMesh mesh;
  if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
    std::cerr << "Error: unable to read input mesh from file " << argv[1]
              << std::endl;
    return 1;
  }

  TriMesh ref_mesh;
  if (!OpenMesh::IO::read_mesh(ref_mesh, argv[2])) {
    std::cerr << "Error: unable to read reference mesh from file " << argv[2]
              << std::endl;
    return 1;
  }

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

  Scalar closeness_weight = 1;
  Scalar laplacian_weight = 0.0;
  Scalar relative_laplacian_weight = 0.1;
  Scalar penalty_parameter = 100000;

  optimize_mesh(mesh, ref_mesh, param.iter, param.anderson_m, penalty_parameter,
                closeness_weight, laplacian_weight, relative_laplacian_weight,
                argv[4]);

  return 0;
}

