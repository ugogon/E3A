#pragma once

#include <Eigen/Core>
#include <list>

#include <gmpxx.h>
typedef mpq_class rational;

namespace Eigen {
typedef Matrix<rational, -1, -1> MatrixXq;
typedef Matrix<rational, -1, 1> VectorXq;
typedef Matrix<rational, -1, 3> MatrixX3q;
typedef Matrix<rational, -1, 2> MatrixX2q;
} // namespace Eigen

template <typename I, typename O>
void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<I> &V, Eigen::MatrixBase<O> &A);

template <typename T> bool non_positive_areas(const Eigen::MatrixBase<T> &A, std::list<int> &I);

template <typename T>
bool non_positive_areas(const Eigen::MatrixBase<T> &A, const std::list<int> &I, std::list<int> &J);

template <typename T> bool non_positive_areas(const Eigen::MatrixBase<T> &A, Eigen::VectorXi &I);
