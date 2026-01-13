// VolumeComp.cu
#include "System.h"
#include "SystemStructures.h"
#include "VolumeComp.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

static double sixV_tet(
    const std::vector<double>& X,
    const std::vector<double>& Y,
    const std::vector<double>& Z,
    int i, int j, int k, int l)
{
    const double ux = X[i] - X[l];
    const double uy = Y[i] - Y[l];
    const double uz = Z[i] - Z[l];

    const double vx = X[j] - X[l];
    const double vy = Y[j] - Y[l];
    const double vz = Z[j] - Z[l];

    const double wx = X[k] - X[l];
    const double wy = Y[k] - Y[l];
    const double wz = Z[k] - Z[l];

    const double cx = vy * wz - vz * wy;
    const double cy = vz * wx - vx * wz;
    const double cz = vx * wy - vy * wx;

    return ux * cx + uy * cy + uz * cz;
}

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs,
    PrismInfoVecs& prismInfoVecs)
{
    const size_t total_nodes = coordInfoVecs.nodeLocX.size();

    // Copy node coordinates from device to host vectors
    std::vector<double> X(total_nodes), Y(total_nodes), Z(total_nodes);
    for (size_t i = 0; i < total_nodes; ++i) {
        X[i] = coordInfoVecs.nodeLocX[i];
        Y[i] = coordInfoVecs.nodeLocY[i];
        Z[i] = coordInfoVecs.nodeLocZ[i];
    }

    const int P = prismInfoVecs.num_prisms;
    if (P <= 0) {
        std::cout << "ComputeVolume: no prisms available.\n";
        generalParams.current_total_volume      = 0.0;
        generalParams.true_current_total_volume = 0.0;
        generalParams.volume_energy             = 0.0;
        return;
    }

    // Copy prism connectivity to host
    thrust::host_vector<int> hP1 = prismInfoVecs.P1;
    thrust::host_vector<int> hP2 = prismInfoVecs.P2;
    thrust::host_vector<int> hP3 = prismInfoVecs.P3;
    thrust::host_vector<int> hP4 = prismInfoVecs.P4;
    thrust::host_vector<int> hP5 = prismInfoVecs.P5;
    thrust::host_vector<int> hP6 = prismInfoVecs.P6;

    double sum6V_signed = 0.0;
    double sum6V_abs    = 0.0;
    
    double min6V_tet = std::numeric_limits<double>::infinity();
    int    min6V_tet_prism = -1;
    
    for (int p = 0; p < P; ++p) {
        const int P1 = hP1[p];
        const int P2 = hP2[p];
        const int P3 = hP3[p];
        const int P4 = hP4[p];
        const int P5 = hP5[p]; 
        const int P6 = hP6[p];
    
        const double s1 = sixV_tet(X, Y, Z, P2, P3, P4, P1);
        const double s2 = sixV_tet(X, Y, Z, P2, P4, P6, P1);
        const double s3 = sixV_tet(X, Y, Z, P4, P5, P6, P1);
    
        // signed global volume
        sum6V_signed += (s1 + s2 + s3);
    
        // diagnostic only (mesh inversion / folding)
        sum6V_abs += (std::fabs(s1) + std::fabs(s2) + std::fabs(s3));
    
        // inversion tracking
        const double local_min = std::min({s1, s2, s3});
        if (local_min < min6V_tet) {
            min6V_tet = local_min;
            min6V_tet_prism = p;
        }
    }
    
    const double V_signed = sum6V_signed / 6.0;
    const double V_abs    = sum6V_abs    / 6.0;
    
    generalParams.current_total_volume = V_signed;
    
    // Keep only as diagnostic (if you keep it at all)
    generalParams.true_current_total_volume = V_abs;
    
//    std::cout << "ComputeVolume: prisms=" << P
//              << ", V_signed=" << V_signed
//              << ", V_abs(sum|tets|)=" << V_abs
//              << ", min6V_tet =" << min6V_tet  // if this value is negative then smth is being inverted. Your volume forces are either garbage or NaNs. 
//              << " (prism " << min6V_tet_prism << ")\n";
}
