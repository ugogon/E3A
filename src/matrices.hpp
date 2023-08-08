#include "graph.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

void construct_Ms(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                  const Eigen::Matrix<int, -1, 2> &Dual_post_order, const Eigen::Matrix<int, -1, 6> &Primal_pre_order,
                  Eigen::SparseMatrix<int, Eigen::RowMajor> &M_1_0, Eigen::SparseMatrix<int, Eigen::RowMajor> &M_1_1,
                  Eigen::SparseMatrix<int, Eigen::RowMajor> &M_2);
