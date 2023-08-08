#include "adjacency.hpp"

int faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v) {
    if (F(f, 0) == v)
        return 0;
    if (F(f, 1) == v)
        return 1;
    if (F(f, 2) == v)
        return 2;
    return -1;
}

int next_faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v) {
    int ret = faceindex(F, f, v);
    if (ret < 0)
        return -1;
    return (ret + 1) % 3;
}

int prev_faceindex(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, int f, int v) {
    int ret = faceindex(F, f, v);
    if (ret < 0)
        return -1;
    return (ret + 3 - 1) % 3;
}

void vt_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V, Eigen::VectorXi &VT, Eigen::VectorXi &VTi) {
    int nv = V.rows() + 1;
    Eigen::VectorXi vfd = Eigen::VectorXi::Zero(nv);
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < 3; j++)
            vfd[F(i, j)]++;
    VTi = Eigen::VectorXi(nv + 1);
    VTi[0] = 0;
    for (int i = 0; i < nv; i++) {
        VTi[i + 1] = VTi[i] + vfd[i];
        vfd[i] = VTi[i];
    }
    VT = Eigen::VectorXi(3 * F.rows());
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < 3; j++)
            VT[vfd[F(i, j)]++] = i;
}

void tt_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F, const Eigen::Ref<const Eigen::VectorXi> VT,
                  const Eigen::Ref<const Eigen::VectorXi> VTi, Eigen::Matrix<int, -1, 3> &TT) {
    TT = Eigen::MatrixXi::Constant(F.rows(), 3, -1);
    for (int f = 0; f < F.rows(); f++)
        for (int k = 0; k < 3; k++) {
            int vi = F(f, (k + 1) % 3), vin = F(f, (k + 2) % 3);
            for (int j = VTi[vi]; j < VTi[vi + 1]; j++) {
                int fn = VT[j];
                if (fn != f)
                    if (F(fn, 0) == vin || F(fn, 1) == vin || F(fn, 2) == vin) {
                        TT(f, k) = fn;
                        break;
                    }
            }
        }
}

// create unordered vertex-vertex adjacency
// makes use of TT structure to correctly handle boundary vertices
void vv_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V,
                  const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> TT, Eigen::VectorXi &VV, Eigen::VectorXi &VVi) {
    int nv = V.rows() + 1;
    Eigen::VectorXi vvd = Eigen::VectorXi::Zero(nv);
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < 3; j++) {
            vvd[F(i, j)]++;
            if (TT(i, j) < 0)
                vvd[F(i, (j + 2) % 3)]++;
        }
    VVi = Eigen::VectorXi(nv + 1);
    VVi[0] = 0;
    for (int i = 0; i < nv; i++) {
        VVi[i + 1] = VVi[i] + vvd[i];
        vvd[i] = VVi[i];
    }
    VV = Eigen::VectorXi(VVi[nv]);
    for (int i = 0; i < F.rows(); i++) {
        VV[vvd[F(i, 0)]++] = F(i, 1);
        VV[vvd[F(i, 1)]++] = F(i, 2);
        VV[vvd[F(i, 2)]++] = F(i, 0);
        if (TT(i, 0) < 0)
            VV[vvd[F(i, 2)]++] = F(i, 1);
        if (TT(i, 1) < 0)
            VV[vvd[F(i, 0)]++] = F(i, 2);
        if (TT(i, 2) < 0)
            VV[vvd[F(i, 1)]++] = F(i, 0);
    }
}

#define FHE(X, Y) (X * 3 + Y)
#define HEF(X) (X / 3)

void he_adjacency(const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> F,
                  const Eigen::Ref<const Eigen::Matrix<double, -1, 3>> V,
                  const Eigen::Ref<const Eigen::Matrix<int, -1, 3>> TT, Eigen::MatrixXi &HEV, Eigen::VectorXi &VHE,
                  Eigen::VectorXi &VHEi, Eigen::VectorXi &HEHE) {
    int nv = V.rows() + 1;
    Eigen::VectorXi vhed = Eigen::VectorXi::Zero(nv);
    int hecnt = F.rows() * 3;
    HEV = Eigen::MatrixXi(hecnt, 2);
    HEHE = Eigen::VectorXi(hecnt);
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            // fill halfedge to vertex map
            HEV.row(i * 3 + j) << F(i, (j + 1) % 3), F(i, (j + 2) % 3);
            // prepare for VEH
            vhed[F(i, j)]++;
            // find twin
            // neighbour triangle
            int f = TT(i, j);
            if (f >= 0) {
                int t = (faceindex(F, f, F(i, (j + 1) % 3)) + 1) % 3;
                HEHE(i * 3 + j) = FHE(f, t);
            } else {
                HEHE(i * 3 + j) = -1;
            }
        }
    }
    // create VEH(VHEi)
    VHEi = Eigen::VectorXi(nv + 1);
    VHEi[0] = 0;
    for (int i = 0; i < nv; i++) {
        VHEi[i + 1] = VHEi[i] + vhed[i];
        vhed[i] = VHEi[i];
    }
    VHE = Eigen::VectorXi(VHEi[nv]);
    for (int i = 0; i < F.rows(); i++) {
        VHE[vhed[F(i, 0)]++] = FHE(i, 2);
        VHE[vhed[F(i, 1)]++] = FHE(i, 0);
        VHE[vhed[F(i, 2)]++] = FHE(i, 1);
        if (TT(i, 0) < 0)
            VHE[vhed[F(i, 2)]++] = HEHE(FHE(i, 0));
        if (TT(i, 1) < 0)
            VHE[vhed[F(i, 0)]++] = HEHE(FHE(i, 1));
        if (TT(i, 2) < 0)
            VHE[vhed[F(i, 1)]++] = HEHE(FHE(i, 2));
    }
}

// generate an (ordered) list of boundary vertices
bool bdy_loop(const Eigen::Ref<Eigen::Matrix<int, -1, 3>> F, const Eigen::Ref<Eigen::Matrix<int, -1, 3>> TT,
              const Eigen::Ref<Eigen::VectorXi> VT, const Eigen::Ref<Eigen::VectorXi> VTi, Eigen::VectorXi &B) {

    // count faces on bdy and remember one
    int bfi = -1, nbf = 0, ki;
    const int fc = F.rows();
    for (int fi = 0; fi < fc; fi++)
        for (int k = 0; k < 3; k++)
            if (TT(fi, k) < 0) {
                nbf++;
                bfi = fi;
                // edge opposite to ki is on boundary so PREV and NEXT are too
                ki = k;
            }
    //   info ("Elements on boundary: " << nbf);
    if (bfi < 0)
        return false;
    B = Eigen::VectorXi(nbf);
    int fi = bfi;
    B(0) = F(fi, PREV(ki));
    for (int i = 1; i < nbf; i++) {
        // emit vertex
        B(i) = F(fi, NEXT(ki));
        // move to next vertex on same face
        int v = F(fi, NEXT(ki));
        // now we need to look for the next edge
        // ie. fi,ki such that TT(fi,ki) = -1 and F(fi,ki+1) = v
        for (int j = VTi[v]; j < VTi[v + 1]; j++)
            for (int k = 0; k < 3; k++)
                if (TT(VT[j], k) < 0)
                    if (F(VT[j], PREV(k)) == v) {
                        fi = VT[j];
                        ki = k;
                        break;
                    }
    }
    return true;
}
