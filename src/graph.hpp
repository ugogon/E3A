#pragma once

#include <Eigen/Core>

struct Graph {
    Eigen::VectorXi B; // Boundary (outer Triangle)
    int root_face;

    Eigen::Matrix<double, -1, 3> V; // Vertices
    Eigen::Matrix<int, -1, 3> F;    // Faces
    Eigen::VectorXi VV;             // Vertex -> Vertex neighbourhood flat list
    Eigen::VectorXi VVi;            // idx in VV indices of neighours of v are VVi[v] to VVi[v+1]
    Eigen::VectorXi VT;             // Vertex -> Face neighbourhood flat list
    Eigen::VectorXi VTi;            // idx in VT indices of neighours of v are VTi[v] to VTi[v+1]
    Eigen::Matrix<int, -1, 3> TT;   // Face -> Face neighbourhood. TT(i,k) is the triangle opposite to F(i,k)
};

Graph build_adjacency_graph(const Eigen::Matrix<int, -1, 3> &F, const Eigen::Matrix<double, -1, 3> &V);

Graph wrap_graph_with_triangle(const Graph &in_graph, const bool padding = false, int keyVertex = -1);

Graph build_graph(const Eigen::Matrix<int, -1, 3> &F, const Eigen::Matrix<double, -1, 3> &V, int keyVertex = -1,
                  bool padding = false);
