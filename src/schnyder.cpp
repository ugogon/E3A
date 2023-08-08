#include "schnyder.hpp"
#include "adjacency.hpp"

#include <iostream>
#include <queue>

#define UNSEEN -1
#define PROCESSED 3
#define RED 0
#define GREEN 1
#define BLUE 2

void schnyder_label(const Graph &graph, Eigen::Matrix<int, -1, 3> &Corners) {
    Corners = Eigen::MatrixXi::Constant(graph.F.rows(), 3, -1);
    int nv = graph.VVi.size() - 1;
    Eigen::VectorXi vs = Eigen::VectorXi::Constant(nv, -1);
    std::queue<int> av;
    for (int k = 0; k < 3; k++) {
        int v = graph.B[k];
        av.push(v);
        vs[v] = k;
    }
    Eigen::VectorXi face_order = Eigen::VectorXi::Constant(graph.F.rows(), -1);
    Eigen::Vector3i cc = Eigen::Vector3i::Ones(3);

    int vi = 0, fi = 0, nc = 0;
    while (!av.empty()) {
        int cvi = av.front();
        av.pop();

        // if vertex has been processed in the meantime ignore it
        if (vs[cvi] == PROCESSED)
            continue;

        // check if cvi can be expanded
        // this is the case if exactly two of its neighbors are on the border
        int cbo = 0;
        for (int i = graph.VVi[cvi]; i < graph.VVi[cvi + 1]; i++)
            if (vs[graph.VV[i]] != UNSEEN && vs[graph.VV[i]] != PROCESSED)
                cbo++;
        if (cbo != 2)
            continue;
        // expand around cvi
        // 1 adjust vertices and correct color counts
        // 1a start with unseen vertices - they inherit the color of cvi
        cc[vs[cvi]]--;
        for (int i = graph.VVi[cvi]; i < graph.VVi[cvi + 1]; i++) {
            // push every vertex on the new boundary (including the ones already on
            // the boundary)
            if (vs[graph.VV[i]] != PROCESSED)
                av.push(graph.VV[i]);
            // set unseen vertex to color of cvi
            if (vs[graph.VV[i]] == UNSEEN) {
                vs[graph.VV[i]] = vs[cvi];
                cc[vs[cvi]]++;
            }
        }
        // 1b check if color of cvi is gone, if so fix
        if (cc[vs[cvi]] == 0) {
            for (int j = graph.VTi[cvi]; j < graph.VTi[cvi + 1]; j++) {
                int f = graph.VT[j];
                int k = faceindex(graph.F, f, cvi);

                int v1 = graph.F(f, (k + 1) % 3);
                int v2 = graph.F(f, (k + 2) % 3);
                if (vs[v1] == UNSEEN || vs[v2] == UNSEEN || vs[v1] == PROCESSED || vs[v2] == PROCESSED ||
                    vs[v1] == vs[v2])
                    continue;
                if (cc[vs[v2]] > cc[vs[v1]])
                    v1 = v2;
                cc[vs[v1]]--;
                vs[v1] = vs[cvi];
                cc[vs[v1]]++;
            }
            if (cc[vs[cvi]] != 1)
                std::cerr << "Boundary coloring could not be fixed!" << std::endl;
        }
        // 2 color corners that are added by this expansion
        // but first assign order number
        for (int i = graph.VTi[cvi]; i < graph.VTi[cvi + 1]; i++) {
            int f = graph.VT[i];
            if (face_order[f] < 0) {
                face_order[f] = fi;
                fi++;
                int idx = faceindex(graph.F, f, cvi);
                Corners(f, idx) = vs[cvi];
                Corners(f, NEXT(idx)) = NEXT(vs[cvi]);
                Corners(f, PREV(idx)) = PREV(vs[cvi]);
            }
        }
        // 3 adjust cvi to be processed
        vs[cvi] = PROCESSED;
    }
}
