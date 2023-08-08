#include "conversion.hpp"
#include "adjacency.hpp"

#include <iostream>
#include <queue>

#include <Eigen/Dense>

#define UNSEEN -1
#define PROCESSED 3
#define RED 0
#define GREEN 1
#define BLUE 2

#define KEYCOLOR 2

// builds Dual Tree for key color and the according primal tree of the other two colors
void build_trees(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 Eigen::Matrix<int, -1, 2> &Dual_post_order, Eigen::Matrix<int, -1, 6> &Primal_pre_order) {
    const int fc = graph.F.rows();
    const int vc = graph.V.rows();
    Eigen::MatrixXi Primal_parents = Eigen::MatrixXi::Constant(vc, 4, -1);

    Dual_post_order.resize(fc, 2);
    int ti = fc - 1;
    int stack_pointer = 0;
    int f = graph.root_face;
    std::vector<int> stack(ti);
    stack[stack_pointer] = f;
    while (stack_pointer >= 0) // runs once per face
    {
        // current face
        f = stack[stack_pointer];
        stack_pointer--;
        int coloridx = faceindex(Corners, f, KEYCOLOR);
        int commonv = graph.F(f, coloridx);
        Dual_post_order(ti, 0) = f;
        Dual_post_order(ti, 1) = graph.TT(f, coloridx);
        ti--;

        int child0 = graph.TT(f, NEXT(coloridx));
        int child1 = graph.TT(f, PREV(coloridx));

        if (child0 >= 0 && graph.TT(child0, faceindex(Corners, child0, KEYCOLOR)) == f) {
            stack_pointer++;
            stack[stack_pointer] = child0;
            int parent = graph.F(f, PREV(coloridx));
            Primal_parents(commonv, 2) = parent;
            Primal_parents(commonv, 3) = child0;
        }
        if (child1 >= 0 && graph.TT(child1, faceindex(Corners, child1, KEYCOLOR)) == f) {
            stack_pointer++;
            stack[stack_pointer] = child1;
            int parent = graph.F(f, NEXT(coloridx));
            Primal_parents(commonv, 0) = parent;
            Primal_parents(commonv, 1) = child1;
        }
    }
    Primal_pre_order.resize(vc - 2, 6);
    for (size_t j = 0; j < 2; j++) { // could be done in parallel
        // root vertices are given by root face
        int v = graph.F(graph.root_face, faceindex(Corners, graph.root_face, (KEYCOLOR + 1 + j) % 3));
        stack_pointer = 0;
        ti = 0;
        stack.clear();
        stack[stack_pointer] = v;
        while (stack_pointer >= 0) // runs once per vertex
        {
            v = stack[stack_pointer];
            stack_pointer--;
            Primal_pre_order(ti, 0 + j * 3) = v;
            Primal_pre_order(ti, 1 + j * 3) = Primal_parents(v, 0 + j * 2);
            Primal_pre_order(ti, 2 + j * 3) = Primal_parents(v, 1 + j * 2);
            ti++;
            for (int i = graph.VVi(v); i < graph.VVi(v + 1); i++) {
                int n = graph.VV[i];
                // if v is parent of n, add to stack
                if (v == Primal_parents(n, 0 + j * 2)) {
                    stack_pointer++;
                    stack[stack_pointer] = n;
                }
            }
        }
    }
}

template <typename In>
void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 const Eigen::Matrix<int, -1, 2> &Dual_post_order, const Eigen::Matrix<int, -1, 6> &Primal_pre_order,
                 const Eigen::VectorX<In> &weights, Eigen::Matrix<In, -1, 2> &bary) {
    Eigen::VectorX<In> weight_sum = weights;
    const int vc = graph.V.rows();
    const int fc = graph.F.rows();
    for (int i = 0; i < fc; i++) {
        int f = Dual_post_order(i, 0);
        int parent = Dual_post_order(i, 1);
        if (parent < 0)
            continue;
        weight_sum(parent) += weight_sum(f);
    }
    bary = Eigen::Matrix<In, -1, 2>::Zero(vc, 2);
    const int ts = Primal_pre_order.rows();
    for (int j = 0; j < 2; j++)
        for (int i = 1; i < ts; i++) {
            // vertex index
            int vi = Primal_pre_order(i, j * 3 + 0);
            // parent index
            int pi = Primal_pre_order(i, j * 3 + 1);
            // index of triangle with summed weight
            int wi = Primal_pre_order(i, j * 3 + 2);
            bary(vi, 1 - j) = weight_sum(wi) + bary(pi, 1 - j);
        }
    In max = weight_sum(graph.root_face);
    for (size_t i = 0; i < 2; i++) {
        bary(graph.B[i], i) = max;
    }
}

template void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 2> &Dual_post_order,
                          const Eigen::Matrix<int, -1, 6> &Primal_pre_order, const Eigen::VectorXi &weights,
                          Eigen::MatrixX2i &bary);

template void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 2> &Dual_post_order,
                          const Eigen::Matrix<int, -1, 6> &Primal_pre_order, const Eigen::VectorX<int64_t> &weights,
                          Eigen::MatrixX2<int64_t> &bary);

template void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 2> &Dual_post_order,
                          const Eigen::Matrix<int, -1, 6> &Primal_pre_order, const Eigen::VectorXq &weights,
                          Eigen::MatrixX2q &bary);

template void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 2> &Dual_post_order,
                          const Eigen::Matrix<int, -1, 6> &Primal_pre_order, const Eigen::VectorXd &weights,
                          Eigen::MatrixX2d &bary);

void build_neighborhood(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                        Eigen::Matrix<int, -1, 6> &Neighbours) {
    // All cases:
    //
    // 1: /k-1\      2: /k+1\
    //   /     \       /     \
    //  /k   k+1\     /k-1  k \
    // b---------a   b---------a
    //  \k-1 k+1/     \k-1 k+1/
    //   \  f  /       \  f  /
    //    \ k /         \ k /
    //     \ /           \ /
    //      c             c
    //
    // 3:  / \b--------a  4: / \b--------a
    //    /k-1\k-1 k+1/     /k+1\k-1 k+1/
    //   /     \  f  /     /     \  f  /
    //  /k   k+1\ k /     /k-1   k\ k /
    // ----------\ /     ----------\ /
    //            c                 c
    //
    // 5:                 6:
    // b--------a/\       b--------a/\
    // \k-1 k+1/k+1\      \k-1 k+1/k-1\
    //  \  f  /     \      \  f  /     \
    //   \ k /k-1   k\      \ k /k   k+1\
    //    \ /---------\      \ /---------\
    //     c                  c
    //
    // Neigbour Layout
    // a_NEXT(kcolor) a_PREV(kcolor) b_NEXT(kcolor) b_PREV(kcolor) c_NEXT(kcolor)
    // c_PREV(kcolor)
    const int fc = graph.F.rows();
    const int vc = graph.V.rows();
    Neighbours = Eigen::MatrixXi::Zero(fc, 6);
    for (int f = 0; f < graph.F.rows(); f++) {
        const int k = faceindex(Corners, f, KEYCOLOR);
        const int parent = graph.TT(f, k);
        const int a = 0;
        const int b = 2;
        const int c = 4;
        if (parent < 0 || graph.TT(parent, faceindex(Corners, parent, PREV(KEYCOLOR))) == f) {
            // case 1
            // b_PREV - a_PREV
            Neighbours(f, a + 1) = -1;
            Neighbours(f, b + 1) = 1;
        } else {
            // case 2
            // a_NEXT - b_NEXT
            Neighbours(f, b + 0) = -1;
            Neighbours(f, a + 0) = 1;
        }
        const int leftchild = graph.TT(f, NEXT(k));
        if (leftchild >= 0 && graph.TT(leftchild, faceindex(Corners, leftchild, KEYCOLOR)) == f) {
            // case 3
            // -(c_NEXT - b_NEXT)
            Neighbours(f, b + 0) += 1;
            Neighbours(f, c + 0) = -1;
        } else {
            // case 4
            // nothing, leaf
        }
        const int rightchild = graph.TT(f, PREV(k));
        if (rightchild >= 0 && graph.TT(rightchild, faceindex(Corners, rightchild, KEYCOLOR)) == f) {
            // case 5
            // -(c_PREV - a_PREV)
            Neighbours(f, a + 1) += 1;
            Neighbours(f, c + 1) = -1;
        } else {
            // case 6
            // nothing, leaf
        }
    }
}

template <typename In>
void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::Matrix<In, -1, 2> &bary,
                 Eigen::VectorX<In> &weights) {
    const int fc = graph.F.rows();
    weights.resize(fc);
    for (int i = 0; i < fc; i++) {
        const int a = NEXT(faceindex(Corners, i, KEYCOLOR));
        weights(i) = 0;
        for (int j = 0; j < 3; j++) {
            int fac = Neighbours(i, j * 2);
            if (fac != 0) {
                weights(i) += fac * bary(graph.F(i, (a + j) % 3), 0);
            }
            fac = Neighbours(i, j * 2 + 1);
            if (fac != 0) {
                weights(i) += fac * bary(graph.F(i, (a + j) % 3), 1);
            }
        }
    }
}

template void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::MatrixX2i &bary,
                          Eigen::VectorXi &weights);

template void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::MatrixX2<int64_t> &bary,
                          Eigen::VectorX<int64_t> &weights);

template void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::MatrixX2q &bary,
                          Eigen::VectorXq &weights);

template void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                          const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::MatrixX2d &bary,
                          Eigen::VectorXd &weights);

void euclidean2bary(const Eigen::MatrixX3d &V, const Eigen::Vector3i &B, Eigen::MatrixX2d &bary) {
    Eigen::Vector3d a, b, c, p, ab, ac, pa, pb, pc;
    a = V.row(B[0]);
    b = V.row(B[1]);
    c = V.row(B[2]);
    ab = b - a;
    ac = c - a;
    double totArea, aArea, bArea;
    totArea = ab.cross(ac).norm();
    bary.resize(V.rows(), 2);
    for (int i = 0; i < V.rows(); i++) {
        p = V.row(i);
        pa = a - p;
        pb = b - p;
        pc = c - p;
        aArea = pb.cross(pc).norm();
        bArea = pc.cross(pa).norm();
        bary(i, 0) = aArea / totArea;
        bary(i, 1) = bArea / totArea;
    }
}

void bary2euclidean(const Eigen::MatrixX3d &V, const Eigen::MatrixX2d &bary, const Eigen::Vector3i &B,
                    Eigen::MatrixX3d &C) {
    Eigen::Vector3d a, b, c;
    a = V.row(B[0]);
    b = V.row(B[1]);
    c = V.row(B[2]);
    C.resize(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++) {
        C.row(i) = bary(i, 0) * a + bary(i, 1) * b + (1 - bary(i, 0) - bary(i, 1)) * c;
    }
}
