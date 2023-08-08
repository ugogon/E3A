#pragma once

#include "area.hpp"
#include "conversion.hpp"
#include "graph.hpp"
#include <Eigen/Core>
#include <iostream>
#include <list>

#include <gmpxx.h>

typedef mpq_class rational;

namespace Eigen {
typedef Matrix<rational, -1, -1> MatrixXq;
typedef Matrix<rational, -1, 1> VectorXq;
typedef Matrix<rational, -1, 3> MatrixX3q;
typedef Matrix<rational, -1, 2> MatrixX2q;
} // namespace Eigen

// fix triangulation based on any given strategy

template <typename T> struct GraphAndTreesInfo {
    Graph &graph;
    Eigen::MatrixX3i &Corners;
    Eigen::MatrixX2i &Dual_post_order;
    Eigen::Matrix<int, -1, 6> &Primal_pre_order;
    Eigen::Matrix<int, -1, 6> &Neighbours;
    Eigen::MatrixX2<T> &dbary;
    Eigen::MatrixX2<T> &fbary;
    Eigen::VectorXi &f_visits;
    int f_its;

    GraphAndTreesInfo(Graph &g, Eigen::MatrixX3i &c, Eigen::MatrixX2i &d, Eigen::Matrix<int, -1, 6> &p,
                      Eigen::Matrix<int, -1, 6> &n, Eigen::MatrixX2<T> &b, Eigen::VectorXi &fv)
        : graph(g), Corners(c), Dual_post_order(d), Primal_pre_order(p), Neighbours(n), dbary(b), fbary(b),
          f_visits(fv) {}
};

template <typename T, typename S> struct Strategy {
    std::function<void(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I,
                       Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew)>
        f;

    explicit Strategy(std::function<void(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A,
                                         const std::list<int> &I, Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew)>
                          f)
        : f(f) {}
};

/**
 * Fixes triangulation contained in info.bary with E3A batch version.
 * It updates all weights corresponding to broken triangles at once, then checks which triangles are still broken
 * for every iteration.
 * @tparam T
 * @tparam S
 * @param strategy fixing strategy for updating weights
 * @param its max number of iterations
 * @param info all triangulation information. Contains info.fbary, where the fixed triangulation will be saved
 * @param weights weights that fix the triangulation
 */
void fix_triangulation(const Strategy<int64_t, rational> &strategy, const int its, GraphAndTreesInfo<int64_t> &info,
                       Eigen::VectorX<int64_t> &weights);
/**
 * Fixes triangulation contained in info.bary with E3A single version.
 * It updates one weight corresponding to a broken triangles at a time, then checks which broken
 * triangles are still broken. The iteration finishes when all triangles that started the iteration broken
 * are fixed.
 * @tparam T
 * @tparam S
 * @param strategy fixing strategy for updating weights
 * @param its max number of iterations
 * @param info all triangulation information. Contains info.fbary, where the fixed triangulation will be saved
 * @param weights weights that fix the triangulation
 */
template <typename T, typename S>
void fix_triangulation_v2(const Strategy<T, S> &strategy, const int its, GraphAndTreesInfo<T> &info,
                          Eigen::VectorX<T> &weights);

// negative --> positive weights

/**
 * Sets all negative weights to 1
 * @param wold old weights
 * @param wnew initialized new weights
 * @param A areas of all triangles
 * @param I indices of triangles to be considered
 */
template <typename T, typename S>
void negWeights2one(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I,
                    Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew);

// increase select triangles' weights

/**
 * Increase weights of triangles by one
 * @param wold old weights
 * @param wnew initialized new weights
 * @param A areas of all triangles
 * @param I indices of triangles to be considered
 */
template <typename T, typename S>
void increase_by_one(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I,
                     Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew);

/**
 * Fully extend triangles
 * @param wold old weights
 * @param wnew initialized new weights vector
 * @param A areas of all triangles
 * @param I indices of triangles to be considered
 * @param info graph and trees information
 */
template <typename T, typename S>
void fully_extend(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I,
                  Eigen::VectorX<T> &wold, Eigen::VectorX<T> &wnew);

/**
 * Unflips triangles by minimal amount
 * @param wold old weights
 * @param wnew initialized new weights vector
 * @param A areas of all triangles
 * @param I indices of triangles to be considered
 * @param info graph and trees information
 */
template <typename T, typename S>
void unflip(GraphAndTreesInfo<T> &info, const Eigen::VectorX<S> &A, const std::list<int> &I, Eigen::VectorX<T> &wold,
            Eigen::VectorX<T> &wnew);
