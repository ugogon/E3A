#include <Eigen/Core>

bool read_triangle_mesh(const std::string str, Eigen::Matrix<double, -1, 3> &V, Eigen::Matrix<int, -1, 3> &F);
