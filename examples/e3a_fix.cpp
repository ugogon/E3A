#include "src/conversion.hpp"
#include "src/fix.hpp"
#include "src/graph.hpp"
#include "src/readmesh.hpp"
#include "src/schnyder.hpp"
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

/**
 * Example script for E3A fixing
 * @param mesh Path to mesh to be fixed
 * @param v E3A version for area checks (batch=1, single=2)
 * @param strat E3A strategy for fixing individual triangles (unflip=0, fully extend=1, increase by one=2)
 * @param res Resolution for E3A discretization (2 ^ res)
 * @param max_its Maximum number of iterations (NOTE: they have different meanings for versions batch and single and
 * thus different effects on the maximum runtime allowed)
 * @param out_mesh Path where to output fixed mesh obj
 */
void e3a_fix(const std::string mesh, int v, int strat, int res, int max_its, const std::string mesh_out) {
    // read mesh
    Eigen::Matrix<double, -1, 3> V;
    Eigen::Matrix<int, -1, 3> F;

    // t_read
    igl::read_triangle_mesh(mesh, V, F);

    // build initial graph. Only necessary for Tutte, so counts towards Tutte embedding time (t_tutte)
    Graph initGraph;
    initGraph = build_graph(F, V);

    // construct tutte if not already embedded
    Eigen::Matrix<double, -1, 2> X_;
    Eigen::MatrixX3d Txyz = Eigen::MatrixXd::Zero(V.rows(), 3);
    X_ = V.leftCols(2);
    Txyz = V;

    // build working graph. Necessary for labeling, so part of building time (t_build)
    Graph graph;
    graph = build_graph(F, Txyz);

    // schnyder label tutte
    Eigen::Matrix<int, -1, 3> Corners;
    schnyder_label(graph, Corners);

    // init trees and other  structures
    Eigen::Matrix<int, -1, 2> Dual_post_order;
    Eigen::Matrix<int, -1, 6> Primal_pre_order;
    build_trees(graph, Corners, Dual_post_order, Primal_pre_order);

    // build neighborhood
    Eigen::Matrix<int, -1, 6> Neighbours;
    build_neighborhood(graph, Corners, Neighbours);

    // convert to barycentric
    Eigen::MatrixX2d bary;
    euclidean2bary(graph.V, graph.B, bary);

    // initialize weights
    Eigen::VectorXd weights;
    bary2weight(graph, Corners, Neighbours, bary, weights);

    // discretize
    Eigen::MatrixX2<int64_t> dbary;
    Eigen::VectorX<int64_t> dweights = (weights * (1LL << res)).cast<int64_t>();
    weight2bary(graph, Corners, Dual_post_order, Primal_pre_order, dweights, dbary);

    int64_t wsum_init = dweights.sum();
    int64_t npow_init = std::ceil(std::log2(abs(wsum_init)));

    // init data structures for fixing
    Eigen::VectorXi visits;
    GraphAndTreesInfo<int64_t> info(graph, Corners, Dual_post_order, Primal_pre_order, Neighbours, dbary, visits);
    Strategy<int64_t, rational> strategy_fe(fully_extend<int64_t, rational>);
    Strategy<int64_t, rational> strategy_u(unflip<int64_t, rational>);
    Strategy<int64_t, rational> strategy_i1(increase_by_one<int64_t, rational>);
    Eigen::MatrixX2<int64_t> fixedBary;
    Eigen::MatrixX3d fixedEuclidean = Eigen::MatrixX3d::Zero(graph.V.rows(), 3);
    int64_t next_power;
    /*====================*/
    if (v == 1) {
        if (strat == 0)
            fix_triangulation(strategy_u, max_its, info, dweights);
        else if (strat == 1)
            fix_triangulation(strategy_fe, max_its, info, dweights);
        else
            fix_triangulation(strategy_i1, max_its, info, dweights);
    } else {
        if (strat == 0)
            fix_triangulation_v2(strategy_u, max_its, info, dweights);
        else if (strat == 1)
            fix_triangulation_v2(strategy_fe, max_its, info, dweights);
        else
            fix_triangulation(strategy_i1, max_its, info, dweights);
    }
    fixedBary = info.fbary;
    next_power = std::ceil(std::log2(abs(dweights.sum())));
    // bary2euclidean(graph.V, fixedBary.cast<double>() / next_power, graph.B, fixedEuclidean);
    fixedEuclidean.leftCols(2) = fixedBary.cast<double>() / (1LL << next_power);

    int64_t wsum_post = dweights.sum();
    int64_t npow_post = std::ceil(std::log2(abs(wsum_post)));
    std::tuple<int64_t, int64_t, int64_t, int64_t> winfo = std::make_tuple(wsum_init, wsum_post, npow_init, npow_post);

    Eigen::MatrixX3d E = fixedEuclidean;
    // write output mes
    igl::write_triangle_mesh(mesh_out, E, graph.F, igl::FileEncoding::Binary);

    std::cout << "[****** fixing stats *******]" << std::endl;
    std::cout << "mesh: " << mesh << std::endl;
    std::cout << "#V: " << V.rows() << std::endl;
    std::cout << "#F: " << F.rows() << std::endl;
    std::cout << "method: E3A" << std::endl;
    std::cout << "version (batch=1,single=2): " << v << std::endl;
    std::cout << "strat: (unflip=0,fully_extend=1): " << strat << std::endl;
    std::cout << "res: " << res << std::endl;
    std::cout << "its: " << info.f_its << " / " << max_its << std::endl;
    std::cout << "wsum_init:npow_init|wsum_post:npow_post " << wsum_init << ":" << npow_init << "|" << wsum_post << ":"
              << npow_post << std::endl;
    std::cout << "========================" << std::endl;
}


int main(int argc, char *argv[]) {
    bool outGiven=false;
    if(argc < 2 || argc > 3) {
        std::cout<<"Usage:\ne3a_fix INPUT [OUTPUT]"<<std::endl;
        return 1;
    }
    std::string path_in(argv[1]);
    std::string path_out = "";
    if(argc == 2) {
        std::string dirname = "";
        std::string basename = "";
        std::string extension = "";
        std::string filename = "";
        igl::pathinfo(path_in, dirname, basename, extension, filename);
        path_out = dirname + "/" + filename + "_OUT" + "." + extension;
    }
    else {
        path_out = argv[3];
    }
    e3a_fix(path_in, 1, 0, 32, 100, path_out);

    return 0;
}

