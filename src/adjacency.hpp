#pragma once

#include <Eigen/Core>

#define NEXT(X) ((X + 1) % 3)
#define PREV(X) ((X + 2) % 3)

int faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v);

int next_faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v);

int prev_faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v);

void vt_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V, Eigen::VectorXi &VT, Eigen::VectorXi &VTi);

void tt_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, const Eigen::Ref<const Eigen::VectorXi> VT,
                  const Eigen::Ref<const Eigen::VectorXi> VTi, Eigen::Matrix<int, -1, 3> &TT);

// create unordered vertex-vertex adjacency
// makes use of TT structure to correctly handle boundary vertices
void vv_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V,
                  const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> TT, Eigen::VectorXi &VV, Eigen::VectorXi &VVi);

bool bdy_loop(const Eigen::Ref<Eigen::Matrix<int, -1, 3>> F, const Eigen::Ref<Eigen::Matrix<int, -1, 3>> TT,
              const Eigen::Ref<Eigen::VectorXi> VT, const Eigen::Ref<Eigen::VectorXi> VTi, Eigen::VectorXi &B);

#define FHE(X, Y) (X * 3 + Y)
#define HEF(X) (X / 3)

void he_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V,
                  const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> TT, Eigen::MatrixXi &HEV, Eigen::VectorXi &VHE,
                  Eigen::VectorXi &VHEi, Eigen::VectorXi &HEHE);
