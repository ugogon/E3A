#pragma once

#include "graph.hpp"

#include <gmpxx.h>
typedef mpq_class rational;

namespace Eigen {
typedef Matrix<rational, -1, -1> MatrixXq;
typedef Matrix<rational, -1, 1> VectorXq;
typedef Matrix<rational, -1, 3> MatrixX3q;
typedef Matrix<rational, -1, 2> MatrixX2q;
} // namespace Eigen

void build_trees(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 Eigen::Matrix<int, -1, 2> &Dual_post_order, Eigen::Matrix<int, -1, 6> &Primal_pre_order);

template <typename In>
void weight2bary(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 const Eigen::Matrix<int, -1, 2> &Dual_post_order, const Eigen::Matrix<int, -1, 6> &Primal_pre_order,
                 const Eigen::Vector<In, -1> &weights, Eigen::Matrix<In, -1, 2> &bary);

void build_neighborhood(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                        Eigen::Matrix<int, -1, 6> &Neighbours);

template <typename In>
void bary2weight(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                 const Eigen::Matrix<int, -1, 6> &Neighbours, const Eigen::Matrix<In, -1, 2> &bary,
                 Eigen::Vector<In, -1> &weights);

void euclidean2bary(const Eigen::MatrixX3d &V, const Eigen::Vector3i &B, Eigen::MatrixX2d &bary);

void bary2euclidean(const Eigen::MatrixX3d &V, const Eigen::MatrixX2d &bary, const Eigen::Vector3i &B,
                    Eigen::MatrixX3d &C);
