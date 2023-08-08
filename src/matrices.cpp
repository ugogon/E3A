#include "matrices.hpp"

void construct_Ms(const Graph &graph, const Eigen::Matrix<int, -1, 3> &Corners,
                  const Eigen::Matrix<int, -1, 2> &Dual_post_order, const Eigen::Matrix<int, -1, 6> &Primal_pre_order,
                  Eigen::SparseMatrix<int, Eigen::RowMajor> &M_1_0, Eigen::SparseMatrix<int, Eigen::RowMajor> &M_1_1,
                  Eigen::SparseMatrix<int, Eigen::RowMajor> &M_2) {
    const int vc = graph.V.rows();
    const int fc = graph.F.rows();
    M_1_0.resize(vc, fc);
    M_1_1.resize(vc, fc);
    M_2.resize(fc, fc);
    M_2.reserve(fc * fc / 2);
    M_2.setIdentity();
    for (int i = 0; i < fc; i++) {
        int f = Dual_post_order(i, 0);
        int parent = Dual_post_order(i, 1);
        if (parent < 0)
            continue;
        M_2.innerVector(parent) += M_2.innerVector(f);
    }
    const int ts = Primal_pre_order.rows();
    for (int i = 1; i < ts; i++) {
        // vertex index
        int vi = Primal_pre_order(i, 0);
        // parent index
        int pi = Primal_pre_order(i, 1);
        // index of triangle with summed weight
        int wi = Primal_pre_order(i, 2);
        M_1_1.innerVector(vi) = M_1_1.innerVector(pi);
        M_1_1.insert(vi, wi) = 1;
        // vertex index
        vi = Primal_pre_order(i, 3);
        // parent index
        pi = Primal_pre_order(i, 4);
        // index of triangle with summed weight
        wi = Primal_pre_order(i, 5);
        M_1_0.innerVector(vi) = M_1_0.innerVector(pi);
        M_1_0.insert(vi, wi) = 1;
    }
    // the barycentric coordiantes of the boundary vertices are not set correctly!
}
