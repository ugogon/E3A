#include "graph.hpp"

#include "adjacency.hpp"
#include <cassert>
#include <igl/triangle/triangulate.h>
#include <igl/write_triangle_mesh.h>
#include <iostream>

void add_triangle_around_bdy(const Eigen::VectorXi &B, Eigen::Matrix<int, -1, 3> &F, Eigen::Matrix<double, -1, 3> &V) {
    int nv = V.rows();
    int nbv = B.rows();
    int fc = F.rows();
    F.conservativeResize(fc + nbv + 3, Eigen::NoChange);
    V.conservativeResize(nv + 3, Eigen::NoChange);
    V.bottomRows<3>() = Eigen::Matrix3d::Zero();
    Eigen::Vector3i c = Eigen::Vector3i::Zero();
    for (int n = 0; n < nbv; n++) {
        int vi = B(n);
        if (n == 0) {
            F.row(fc++) << vi, nv, nv + 2;
            V.row(nv + 2) += V.row(vi);
            c[2]++;
        } else if (n == nbv / 3) {
            F.row(fc++) << vi, nv + 1, nv;
            V.row(nv) += V.row(vi);
            c[0]++;
        } else if (n == 2 * nbv / 3) {
            F.row(fc++) << vi, nv + 2, nv + 1;
            V.row(nv + 1) += V.row(vi);
            c[1]++;
        }

        if (n < nbv / 3) {
            F.row(fc++) << B(n), B(n + 1), nv;
            V.row(nv) += V.row(vi);
            c[0]++;
        } else if (n < 2 * nbv / 3) {
            F.row(fc++) << B(n), B(n + 1), nv + 1;
            V.row(nv + 1) += V.row(vi);
            c[1]++;
        } else {
            F.row(fc++) << B(n), B((n + 1) % nbv), nv + 2;
            V.row(nv + 2) += V.row(vi);
            c[2]++;
        }
    }
    Eigen::Vector3d sum = Eigen::Vector3d::Zero(3);
    for (int i = 0; i < 3; i++) {
        V.row(nv + i) *= 1.0 / (double)c[i];
        sum += V.row(nv + i);
    }
    sum /= 3;
    for (int i = 0; i < 3; i++) {
        V.row(nv + i) = (V.row(nv + i) - sum.transpose()) * 5 + sum.transpose();
    }
}

int findIdxInBoundary(const Eigen::VectorXi &oB, int keyVertex) {
    assert(oB.rows() == 3);
    int idx = -1;
    if (oB[0] == keyVertex)
        idx = 0;
    else if (oB[1] == keyVertex)
        idx = 1;
    else if (oB[2] == keyVertex)
        idx = 2;
    return idx;
}

void build_adjacency(const Eigen::Matrix<int, -1, 3> &F, const Eigen::Matrix<double, -1, 3> &V, Eigen::VectorXi &oVT,
                     Eigen::VectorXi &oVTi, Eigen::Matrix<int, -1, 3> &oTT, Eigen::VectorXi &oVV,
                     Eigen::VectorXi &oVVi) {
    vt_adjacency(F, V, oVT, oVTi);
    tt_adjacency(F, oVT, oVTi, oTT);
    vv_adjacency(F, V, oTT, oVV, oVVi);
}

void set_root_face(Graph &graph) {
    // find root face (triangle incident to graph.B[0] and graph.B[1])
    for (int i = graph.VTi(graph.B[0]); i < graph.VTi(graph.B[0] + 1); i++) {
        int nt = graph.VT[i];
        // if v is parent of n, add to stack
        if (graph.F(nt, 0) == graph.B[1] || graph.F(nt, 1) == graph.B[1] || graph.F(nt, 2) == graph.B[1]) {
            graph.root_face = nt;
            break;
        }
    }
}

Graph wrap_graph_with_triangle(const Graph &in_graph, const bool padding, int keyVertex) {
    Graph out_graph;
    out_graph.B.resize(3);
    size_t bsize = in_graph.B.size();
    if (bsize == 3) {
        if (keyVertex == -1)
            keyVertex = in_graph.B[0];
        int idx = findIdxInBoundary(in_graph.B, keyVertex);
        out_graph.F = in_graph.F;
        out_graph.V = in_graph.V;
        out_graph.B << in_graph.B[(idx + 2) % 3], in_graph.B[(idx + 1) % 3], in_graph.B[idx];
        out_graph.VT = in_graph.VT;
        out_graph.VTi = in_graph.VTi;
        out_graph.TT = in_graph.TT;
        out_graph.VV = in_graph.VV;
        out_graph.VVi = in_graph.VVi;
    } else {
        // update V,F and B
        if (!padding) {
            out_graph.V = in_graph.V;
            out_graph.F = in_graph.F;
            add_triangle_around_bdy(in_graph.B, out_graph.F, out_graph.V);
            const int nv = out_graph.V.rows();
            out_graph.B << nv - 1, nv - 2, nv - 3;
        } else {
            const int nv = in_graph.V.rows();
            // calculate wrapping boundary triangle
            Eigen::MatrixXd rV;
            rV.resize(bsize + 3, 2);
            Eigen::Vector2d min = in_graph.V.leftCols(2).colwise().minCoeff();
            Eigen::Vector2d max = in_graph.V.leftCols(2).colwise().maxCoeff();
            Eigen::Vector2d mean = in_graph.V.leftCols(2).colwise().mean();

            Eigen::Vector2d off = (max - min) * 0.1;
            min -= off;
            max += off;
            auto diff = max - min;

            for (size_t i = 0; i < bsize; i++) {
                rV.row(i) << in_graph.V.row(in_graph.B(i)).leftCols(2);
            }
            // triangle to V boundary
            rV.row(bsize + 0) << min(0), max(1) + diff(1);
            rV.row(bsize + 1) << max(0) + diff(0), min(1);
            rV.row(bsize + 2) = min;

            // define boundary edges
            Eigen::MatrixXi E;
            E.resize(bsize + 3, 2);
            for (size_t i = 0; i < bsize; i++) {
                E.row(i) << i, (i + 1) % bsize;
            }
            E.row(bsize + 0) << bsize, bsize + 1;
            E.row(bsize + 1) << bsize + 1, bsize + 2;
            E.row(bsize + 2) << bsize + 2, bsize;

            // insert hole
            Eigen::MatrixXd H;
            H.resize(1, 2);
            H << mean;

            Eigen::MatrixXd V = Eigen::MatrixXd::Zero(rV.rows(), 3);
            V.leftCols(2) = rV;

            Eigen::MatrixXi NewF;
            Eigen::MatrixXd NewV;
            igl::triangle::triangulate(rV, E, H, "pa0.001Y", NewV, NewF);
            Eigen::MatrixX3d test = Eigen::MatrixXd::Zero(NewV.rows(), 3);
            test.leftCols(2) = NewV.leftCols(2);
            igl::write_triangle_mesh("hole.off", test, NewF);

            // add new vertices to graph
            out_graph.V = Eigen::MatrixXd::Zero(nv + NewV.rows() - bsize, 3);
            out_graph.V.block(0, 0, nv, 2) = in_graph.V;
            out_graph.V.block(nv, 0, NewV.rows() - bsize, 2) = NewV.block(bsize, 0, NewV.rows() - bsize, 2);
            // adjust new Vertex ids
            for (size_t i = 0; i < (size_t)NewF.rows(); i++) {
                for (size_t j = 0; j < 3; j++) {
                    size_t id = NewF(i, j);
                    if (id < bsize) {
                        NewF(i, j) = in_graph.B[id];
                    } else if (id >= bsize) {
                        NewF(i, j) = NewF(i, j) - bsize + nv;
                    }
                }
            }
            out_graph.B << nv + 1, nv + 0, nv + 2;

            // add new faces to graph
            int nf = in_graph.F.rows();
            out_graph.F.resize(nf + NewF.rows(), 3);
            out_graph.F.block(0, 0, nf, 3) = in_graph.F;
            out_graph.F.block(nf, 0, NewF.rows(), 3) = NewF;
        }
        // TODO: Do efficiently
        build_adjacency(out_graph.F, out_graph.V, out_graph.VT, out_graph.VTi, out_graph.TT, out_graph.VV,
                        out_graph.VVi);
    }
    set_root_face(out_graph);
    return out_graph;
}

Graph build_adjacency_graph(const Eigen::Matrix<int, -1, 3> &F, const Eigen::Matrix<double, -1, 3> &V) {
    Graph graph;
    graph.F = F;
    graph.V = V;
    build_adjacency(graph.F, graph.V, graph.VT, graph.VTi, graph.TT, graph.VV, graph.VVi);
    bdy_loop(graph.F, graph.TT, graph.VT, graph.VTi, graph.B);
    return graph;
}

Graph build_graph(const Eigen::Matrix<int, -1, 3> &F, const Eigen::Matrix<double, -1, 3> &V, int keyVertex,
                  bool padding) {
    Graph graph;
    graph = build_adjacency_graph(F, V);
    return wrap_graph_with_triangle(graph, padding, keyVertex);
}
