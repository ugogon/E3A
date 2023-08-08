#include "area.hpp"
#include "graph.hpp"

template <typename In, typename Out>
void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<In> &V, Eigen::MatrixBase<Out> &A) {
    typedef typename Out::Scalar Scalar;
    int n = F.rows();
    A = Out(n, 1);
    if (V.cols() == 2)
        for (int i = 0; i < n; i++) {
            Scalar x0 = V(F(i, 0), 0), y0 = V(F(i, 0), 1);
            Scalar x1 = V(F(i, 1), 0), y1 = V(F(i, 1), 1);
            Scalar x2 = V(F(i, 2), 0), y2 = V(F(i, 2), 1);
            Scalar x10 = (x1 - x0);
            Scalar x20 = (x2 - x0);
            Scalar y10 = (y1 - y0);
            Scalar y20 = (y2 - y0);
            // A(i,0) = (x2*(y0 - y1) + x0*(y1 - y2) + x1*(y2 - y0))/2;
            A(i, 0) = (x10 * y20) - (x20 * y10);
        }
    else
        // info("Vertex coordinate format not supported");
        return;
}

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2d> &V,
                    Eigen::MatrixBase<Eigen::VectorXd> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2d> &V,
                    Eigen::MatrixBase<Eigen::VectorXq> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2f> &V,
                    Eigen::MatrixBase<Eigen::VectorXf> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2i> &V,
                    Eigen::MatrixBase<Eigen::VectorXi> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2i> &V,
                    Eigen::MatrixBase<Eigen::VectorXq> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2<int64_t>> &V,
                    Eigen::MatrixBase<Eigen::VectorX<int64_t>> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2<int64_t>> &V,
                    Eigen::MatrixBase<Eigen::VectorXq> &A);

template void areas(const Eigen::Matrix<int, -1, 3> &F, const Eigen::MatrixBase<Eigen::MatrixX2q> &V,
                    Eigen::MatrixBase<Eigen::VectorXq> &A);

// TODO: only consider possibly affected triangles (i.e. those in dual tree)
template <typename T> bool non_positive_areas(const Eigen::MatrixBase<T> &A, std::list<int> &I) {
    bool found = false;
    I.clear();
    for (int i = 0; i < A.rows(); i++) {
        if (A(i, 0) <= 0) {
            found = true;
            I.push_back(i);
        }
    }
    return found;
}

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXd> &A, std::list<int> &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXf> &A, std::list<int> &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXi> &A, std::list<int> &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorX<int64_t>> &A, std::list<int> &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXq> &A, std::list<int> &I);

// TODO: only consider possibly affected triangles (i.e. those in dual tree)
template <typename T> bool non_positive_areas(const Eigen::MatrixBase<T> &A, Eigen::VectorXi &I) {
    bool found = false;
    I.setZero(A.rows());
    for (int i = 0; i < A.rows(); i++) {
        if (A(i, 0) <= 0) {
            found = true;
            I[i] = 1;
        }
    }
    return found;
}

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXd> &A, Eigen::VectorXi &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXf> &A, Eigen::VectorXi &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXi> &A, Eigen::VectorXi &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorX<int64_t>> &A, Eigen::VectorXi &I);

template bool non_positive_areas(const Eigen::MatrixBase<Eigen::VectorXq> &A, Eigen::VectorXi &I);
