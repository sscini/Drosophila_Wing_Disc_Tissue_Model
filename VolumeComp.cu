// VolumeComp.cu

#include "System.h"
#include "SystemStructures.h"
#include "VolumeComp.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

// 6 * signed volume of tetrahedron (i,j,k,l) using scalar triple product
static double sixV_tet(
    const std::vector<double>& X,
    const std::vector<double>& Y,
    const std::vector<double>& Z,
    int i, int j, int k, int l) {

    const double ux = X[i] - X[l];
    const double uy = Y[i] - Y[l];
    const double uz = Z[i] - Z[l];

    const double vx = X[j] - X[l];
    const double vy = Y[j] - Y[l];
    const double vz = Z[j] - Z[l];

    const double wx = X[k] - X[l];
    const double wy = Y[k] - Y[l];
    const double wz = Z[k] - Z[l];

    // v × w
    const double cx = vy * wz - vz * wy;
    const double cy = vz * wx - vx * wz;
    const double cz = vx * wy - vy * wx;

    // u · (v × w) = 6 * V_tet_signed
    return ux * cx + uy * cy + uz * cz;
}

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs,
    PrismInfoVecs& prismInfoVecs) {

    // Copy node coordinates from device to host vectors
    const size_t total_nodes = coordInfoVecs.nodeLocX.size();
    std::vector<double> X(total_nodes), Y(total_nodes), Z(total_nodes);

    for (size_t i = 0; i < total_nodes; ++i) {
        X[i] = coordInfoVecs.nodeLocX[i];
        Y[i] = coordInfoVecs.nodeLocY[i];
        Z[i] = coordInfoVecs.nodeLocZ[i];
    }

    // This version assumes exactly 2 layers with the same number of nodes
    if (generalParams.num_layers <= 0 ||
        total_nodes < 8 ||
        total_nodes % generalParams.num_layers != 0) {

        std::cout << "ComputeVolume: unexpected node count " << total_nodes
                  << " or num_layers = " << generalParams.num_layers << "\n";
        generalParams.current_total_volume      = 0.0;
        generalParams.true_current_total_volume = 0.0;
        generalParams.volume_energy             = 0.0;
        return;
    }

//    if (generalParams.num_layers != 2) {
//        std::cout << "ComputeVolume: this implementation assumes 2 layers, got "
//                  << generalParams.num_layers << "\n";
//        generalParams.current_total_volume      = 0.0;
//        generalParams.true_current_total_volume = 0.0;
//        generalParams.volume_energy             = 0.0;
//        return;
//    }

    const int nodes_per_layer = static_cast<int>(total_nodes / generalParams.num_layers);
    const int top_start       = 0;
    const int bot_start       = nodes_per_layer;

    // ---- find "center" node in each layer as the one with max z ----
    auto find_center = [&](int start) {
        int idx   = start;
        double mz = Z[start];
        for (int i = start + 1; i < start + nodes_per_layer; ++i) {
            if (Z[i] > mz) {
                mz = Z[i];
                idx = i;
            }
        }
        return idx;
    };

    const int topC = find_center(top_start);
    const int botC = find_center(bot_start);

    // ---- build ring order by polar angle around each center ----
    struct NodeAngle {
        int idx;
        double ang;
    };

    std::vector<NodeAngle> topRing;
    std::vector<NodeAngle> botRing;

    // top layer ring
    {
        const double cx = X[topC];
        const double cy = Y[topC];

        for (int i = top_start; i < top_start + nodes_per_layer; ++i) {
            if (i == topC) continue;
            const double ang = std::atan2(Y[i] - cy, X[i] - cx);
            topRing.push_back({i, ang});
        }
        std::sort(topRing.begin(), topRing.end(),
                  [](const NodeAngle& a, const NodeAngle& b) {
                      return a.ang < b.ang;
                  });
    }

    // bottom layer ring
    {
        const double cx = X[botC];
        const double cy = Y[botC];

        for (int i = bot_start; i < bot_start + nodes_per_layer; ++i) {
            if (i == botC) continue;
            const double ang = std::atan2(Y[i] - cy, X[i] - cx);
            botRing.push_back({i, ang});
        }
        std::sort(botRing.begin(), botRing.end(),
                  [](const NodeAngle& a, const NodeAngle& b) {
                      return a.ang < b.ang;
                  });
    }

    const int sectors = static_cast<int>(topRing.size());
    if (sectors != static_cast<int>(botRing.size()) || sectors == 0) {
        std::cout << "ComputeVolume: ring mismatch, top sectors = "
                  << sectors << ", bottom sectors = " << botRing.size() << "\n";
        generalParams.current_total_volume      = 0.0;
        generalParams.true_current_total_volume = 0.0;
        generalParams.volume_energy             = 0.0;
        return;
    }

    // ---- 6-sector, 3-tet decomposition (P1..P6) ----
    double sum6V_signed = 0.0;
    double sum6V_abs    = 0.0;

    for (int k = 0; k < sectors; ++k) {
        const int top0 = topRing[k].idx;
        const int top1 = topRing[(k + 1) % sectors].idx;
        const int bot0 = botRing[k].idx;
        const int bot1 = botRing[(k + 1) % sectors].idx;

        // Match your P1..P6 notation:
        // P1 = bottom center, P2/P3 = bottom ring nodes,
        // P4 = top center,    P5/P6 = top ring nodes.
        const int P1 = botC;
        const int P2 = bot0;
        const int P3 = bot1;
        const int P4 = topC;
        const int P5 = top0;
        const int P6 = top1;

        const double s1_signed = sixV_tet(X, Y, Z, P2, P3, P4, P1); // bottom
        const double s2_signed = sixV_tet(X, Y, Z, P2, P4, P6, P1); // side
        const double s3_signed = sixV_tet(X, Y, Z, P4, P5, P6, P1); // top

        sum6V_signed += (s1_signed + s2_signed + s3_signed);

        sum6V_abs += (std::fabs(s1_signed) +
                      std::fabs(s2_signed) +
                      std::fabs(s3_signed));
    }

    const double V_signed = sum6V_signed / 6.0;
    const double V_abs    = sum6V_abs    / 6.0;

    // store both: signed for orientation, abs for "physical" size
    generalParams.current_total_volume      = V_signed;
    generalParams.true_current_total_volume = V_abs;

    const double Kv = generalParams.volume_spring_constant;

    // If you want a volume energy term, you can uncomment and use V_abs, V0:
    // const double V0 = generalParams.eq_total_volume;
    // const double dV = V_abs - V0;
    // generalParams.volume_energy =
    //     (V0 > 0.0) ? (Kv * dV * dV / (2.0 * V0)) : 0.0;

    std::cout << "ComputeVolume: sectors = " << sectors
              << ", V_signed = " << V_signed
              << ", V_abs = "    << V_abs
              << std::endl;
}
