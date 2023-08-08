#include "readmesh.hpp"
#include <igl/read_triangle_mesh.h>

bool read_triangle_mesh(const std::string str, Eigen::Matrix<double, -1, 3> &V, Eigen::Matrix<int, -1, 3> &F) {
    return igl::read_triangle_mesh(str, V, F);
}
