#include "fix.hpp"

#include <iostream>

#include "area.hpp"
#include "conversion.hpp"

#include <igl/write_triangle_mesh.h>
#include <list>
#include <valarray>

int parent_dual_tree(const Graph &graph, const Eigen::MatrixX3i &Corners, int t, int c) {
    for (int i = 0; i < 3; i++) {
        if (Corners(t, i) == c) {
            return graph.TT(t, i);
        }
    }
    return -1;
}

// Updates the areas of all faces along the dual paths of faces in T to list.
// These are the only faces affected by updating each of T's ts weight.
template <typename T, typename S>
void update_dual_path_areas(GraphAndTreesInfo<T> &info, int t, Eigen::VectorX<S> &A, std::list<int> &I) {
    int pc;
    S a, x0, y0, x1, x2, y1, y2, x10, x20, y10, y20;
    x0 = info.fbary(info.graph.F(t, 0), 0), y0 = info.fbary(info.graph.F(t, 0), 1);
    x1 = info.fbary(info.graph.F(t, 1), 0), y1 = info.fbary(info.graph.F(t, 1), 1);
    x2 = info.fbary(info.graph.F(t, 2), 0), y2 = info.fbary(info.graph.F(t, 2), 1);
    x10 = (x1 - x0);
    x20 = (x2 - x0);
    y10 = (y1 - y0);
    y20 = (y2 - y0);
    A[t] = (x10 * y20) - (x20 * y10);
    if (A[t] <= 0)
        I.push_back(t);
    for (int c = 0; c < 3; c++) {
        pc = parent_dual_tree(info.graph, info.Corners, t, c);
        while (pc != -1) {
            x0 = info.fbary(info.graph.F(pc, 0), 0), y0 = info.fbary(info.graph.F(pc, 0), 1);
            x1 = info.fbary(info.graph.F(pc, 1), 0), y1 = info.fbary(info.graph.F(pc, 1), 1);
            x2 = info.fbary(info.graph.F(pc, 2), 0), y2 = info.fbary(info.graph.F(pc, 2), 1);
            x10 = (x1 - x0);
            x20 = (x2 - x0);
            y10 = (y1 - y0);
            y20 = (y2 - y0);
            a = (x10 * y20) - (x20 * y10);
            if (A[pc] > 0 && a <= 0) { // fixing t broke pc (+ => -)
                I.push_back(pc);
            }
            A[pc] = a; // (- => -) (+ => +) (- => +)
            pc = parent_dual_tree(info.graph, info.Corners, pc, c);
        }
    }
}

void fix_triangulation(const Strategy<int64_t, rational> &strategy, const int its, GraphAndTreesInfo<int64_t> &info,
                       Eigen::VectorX<int64_t> &weights) {
    Eigen::VectorX<int64_t> wold, wnew;
    // calculate initial weights
    bary2weight(info.graph, info.Corners, info.Neighbours, info.dbary, wold);
    wnew = wold;
    Eigen::VectorX<rational> A;
    std::list<int> I;
    info.fbary = info.dbary;
    info.f_visits = Eigen::VectorXi::Zero(info.graph.F.rows());
    areas(info.graph.F, info.fbary, A);
    int it = 0;
    while (non_positive_areas(A, I) && it < its) {
        strategy.f(info, A, I, wold, wnew);
        weight2bary(info.graph, info.Corners, info.Dual_post_order, info.Primal_pre_order, wnew, info.fbary);
        {
            std::cout << "one step done" << std::endl;
            Eigen::MatrixX3d fixedEuclidean = Eigen::MatrixX3d::Zero(info.graph.V.rows(), 3);
            fixedEuclidean.leftCols(2) = info.fbary.cast<double>();
            std::string filename = "fixed" + std::to_string(it);
            igl::write_triangle_mesh(filename + ".obj", fixedEuclidean, info.graph.F);
        }
        areas(info.graph.F, info.fbary, A);
        it++;
    }
    info.f_its = it;
    weights = wnew;
}

template <typename T, typename S>
void fix_triangulation_v2(const Strategy<T, S> &strategy, const int its, GraphAndTreesInfo<T> &info,
                          Eigen::VectorX<T> &weights) {
    Eigen::VectorX<T> wold, wnew;
    // calculate initial weights
    bary2weight(info.graph, info.Corners, info.Neighbours, info.dbary, wold);
    wnew = wold;
    Eigen::VectorX<S> A;
    std::list<int> I, idx;
    info.fbary = info.dbary;
    info.f_visits = Eigen::VectorXi::Zero(info.graph.F.rows());
    areas(info.graph.F, info.fbary, A);
    non_positive_areas(A, I);
    int it = 0;
    while (!I.empty() && it < its) {
        idx.clear();
        idx.push_front(I.front());
        I.pop_front();
        strategy.f(info, A, idx, wold, wnew);
        weight2bary(info.graph, info.Corners, info.Dual_post_order, info.Primal_pre_order, wnew, info.fbary);
        update_dual_path_areas(info, idx.front(), A, I);
        it++;
    }
    info.f_its = it;
    weights = wnew;
}

template void fix_triangulation_v2(const Strategy<int, rational> &strategy, const int its, GraphAndTreesInfo<int> &info,
                                   Eigen::VectorXi &weights);

template void fix_triangulation_v2(const Strategy<int64_t, rational> &strategy, const int its,
                                   GraphAndTreesInfo<int64_t> &info, Eigen::VectorX<int64_t> &weights);

template void fix_triangulation_v2(const Strategy<double, rational> &strategy, const int its,
                                   GraphAndTreesInfo<double> &info, Eigen::VectorXd &weights);

template <typename T, typename S>
void negWeights2one(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I,
                    Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew) {
    Eigen::VectorX<T> wnew_copy = wnew;
    wold = wnew;
    for (int f = 0; f < wold.rows(); f++) {
        if (wold[f] < 0)
            wnew[f] = 1;
        else
            wnew[f] = wold[f];
    }
    wold = wnew_copy;
}

template <>
void increase_by_one<int64_t>(GraphAndTreesInfo<int64_t> &info, const Eigen::VectorXq &A, const std::list<int> &I,
                              Eigen::VectorX<int64_t> &wold, Eigen::VectorX<int64_t> &wnew) {
    Eigen::VectorX<int64_t> wnew_copy = wnew;
    wold = wnew;
    int64_t delta = 1;
    for (int idx : I) {
        if (A[idx] > 0)
            continue;
        wnew[idx] = wold[idx] + delta;
    }
    wold = wnew_copy;
}

template <>
void fully_extend<int64_t, rational>(GraphAndTreesInfo<int64_t> &info, const Eigen::VectorXq &A,
                                     const std::list<int> &I, Eigen::VectorX<int64_t> &wold,
                                     Eigen::VectorX<int64_t> &wnew) {
    Eigen::VectorX<int64_t> wnew_copy = wnew;
    wold = wnew;
    int a, b, c;
    int64_t w_sum = wold.sum(), delta_w_idx;
    Eigen::Vector3<int64_t> pa, pb, pc, xs, ys, zs, ds;
    for (int idx : I) {
        if (A[idx] > 0)
            continue;
        info.f_visits[idx]++;
        for (int i = 0; i < 3; i++) {
            if (info.Corners(idx, i) == 0) {
                a = info.graph.F(idx, i);
                b = info.graph.F(idx, (i + 1) % 3);
                c = info.graph.F(idx, (i + 2) % 3);
                break;
            }
        }
        pa << info.fbary(a, 0), info.fbary(a, 1), w_sum - info.fbary(a, 0) - info.fbary(a, 1);
        pb << info.fbary(b, 0), info.fbary(b, 1), w_sum - info.fbary(b, 0) - info.fbary(b, 1);
        pc << info.fbary(c, 0), info.fbary(c, 1), w_sum - info.fbary(c, 0) - info.fbary(c, 1);
        xs << pa[0], pb[0], pc[0];
        ys << pa[1], pb[1], pc[1];
        zs << pa[2], pb[2], pc[2];
        ds << xs.maxCoeff() - xs[0], ys.maxCoeff() - ys[1], zs.maxCoeff() - zs[2];
        delta_w_idx = ds.maxCoeff() + 1;
        wnew[idx] = wold[idx] + delta_w_idx;
    }
    wold = wnew_copy;
}

template <>
void unflip<int64_t, rational>(GraphAndTreesInfo<int64_t> &info, const Eigen::VectorXq &A, const std::list<int> &I,
                               Eigen::VectorX<int64_t> &wold, Eigen::VectorX<int64_t> &wnew) {
    Eigen::VectorX<int64_t> wnew_copy = wnew;
    wold = wnew;
    int a, b, c;
    int64_t delta_w_idx, xa, xc, yb, yc, tr;
    Eigen::Vector2<int64_t> pa, pb, pc;
    Eigen::Vector3<int64_t> xs, ys;
    for (int idx : I) {
        if (A[idx] > 0)
            continue;
        info.f_visits[idx]++;
        for (int i = 0; i < 3; i++) {
            if (info.Corners(idx, i) == 0) {
                a = info.graph.F(idx, i);
                b = info.graph.F(idx, (i + 1) % 3);
                c = info.graph.F(idx, (i + 2) % 3);
                break;
            }
        }
        // A = [ xa - xc    xb - xc ]
        //     [ ya - yc    yb - yc ]
        xa = info.fbary(a, 0);
        xc = info.fbary(c, 0);
        yb = info.fbary(b, 1);
        yc = info.fbary(c, 1);
        tr = ((xa - xc) + (yb - yc));
        delta_w_idx = (int64_t)std::floor(-tr / 2.0 + std::sqrt(std::pow(tr / 2.0, 2) - A[idx].get_d()) + 1.0);
        wnew[idx] = wold[idx] + delta_w_idx;
    }
    wold = wnew_copy;
}
