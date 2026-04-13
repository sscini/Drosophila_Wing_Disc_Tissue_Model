//// VolumeComp.cu
//// Computes total enclosed volume and per-prism volume energy
////
//// IMPORTANT: The energy computed here MUST match the force law in VolumeSprings.cu.
//// VolumeSprings.cu uses per-prism local conservation:
////     E = sum_p  k_v * (V_p - V0_p)^2
//// So the energy reported here must be the same sum, NOT the global
////     E = k_v * (V_total - V0_total)^2
//// which would be inconsistent with the forces.
//
//#include "VolumeComp.h"
//#include "System.h"
//#include "SystemStructures.h"
//
//#include <vector>
//#include <algorithm>
//#include <cmath>
//#include <iostream>
//#include <limits>
//
///**
// * Compute 6 times the signed volume of tetrahedron (i,j,k,l)
// *
// * Uses scalar triple product: 6V = u . (v x w)
// * where u = r_i - r_l, v = r_j - r_l, w = r_k - r_l
// *
// * Positive volume when vertices are in right-hand orientation.
// */
//static double sixV_tet(
//    const std::vector<double>& X,
//    const std::vector<double>& Y,
//    const std::vector<double>& Z,
//    int i, int j, int k, int l)
//{
//    const double ux = X[i] - X[l];
//    const double uy = Y[i] - Y[l];
//    const double uz = Z[i] - Z[l];
//
//    const double vx = X[j] - X[l];
//    const double vy = Y[j] - Y[l];
//    const double vz = Z[j] - Z[l];
//
//    const double wx = X[k] - X[l];
//    const double wy = Y[k] - Y[l];
//    const double wz = Z[k] - Z[l];
//
//    const double cx = vy * wz - vz * wy;
//    const double cy = vz * wx - vx * wz;
//    const double cz = vx * wy - vy * wx;
//
//    return ux * cx + uy * cy + uz * cz;
//}
//
//void ComputeVolume(
//    GeneralParams& generalParams,
//    CoordInfoVecs& coordInfoVecs,
//    LinearSpringInfoVecs& linearSpringInfoVecs,
//    LJInfoVecs& ljInfoVecs,
//    PrismInfoVecs& prismInfoVecs)
//{
//    const size_t total_nodes = coordInfoVecs.nodeLocX.size();
//
//    // Copy node coordinates from device to host
//    std::vector<double> X(total_nodes), Y(total_nodes), Z(total_nodes);
//    for (size_t i = 0; i < total_nodes; ++i) {
//        X[i] = coordInfoVecs.nodeLocX[i];
//        Y[i] = coordInfoVecs.nodeLocY[i];
//        Z[i] = coordInfoVecs.nodeLocZ[i];
//    }
//
//    const int P = prismInfoVecs.num_prisms;
//    if (P <= 0) {
//        std::cout << "ComputeVolume: no prisms available.\n";
//        generalParams.current_total_volume = 0.0;
//        generalParams.true_current_total_volume = 0.0;
//        generalParams.volume_energy = 0.0;
//        return;
//    }
//
//    // Copy prism connectivity to host
//    thrust::host_vector<int> hP1 = prismInfoVecs.P1;
//    thrust::host_vector<int> hP2 = prismInfoVecs.P2;
//    thrust::host_vector<int> hP3 = prismInfoVecs.P3;
//    thrust::host_vector<int> hP4 = prismInfoVecs.P4;
//    thrust::host_vector<int> hP5 = prismInfoVecs.P5;
//    thrust::host_vector<int> hP6 = prismInfoVecs.P6;
//
//    // Check if per-prism equilibrium volumes are available
//    const bool have_eq_prism_volumes =
//        ((int)prismInfoVecs.eq_prism_volume.size() == P);
//
//    // If available, copy equilibrium volumes to host
//    thrust::host_vector<double> hV0;
//    if (have_eq_prism_volumes) {
//        hV0 = prismInfoVecs.eq_prism_volume;
//    }
//
//    const double kv = generalParams.volume_spring_constant;
//
//    double sum6V_signed = 0.0;
//    double sum6V_abs = 0.0;
//    double per_prism_energy_sum = 0.0;
//
//    double min6V_tet_val = std::numeric_limits<double>::infinity();
//    int min6V_tet_prism = -1;
//
//    for (int p = 0; p < P; ++p) {
//        const int p1 = hP1[p];
//        const int p2 = hP2[p];
//        const int p3 = hP3[p];
//        const int p4 = hP4[p];
//        const int p5 = hP5[p];
//        const int p6 = hP6[p];
//
//        // Decompose prism into 3 tetrahedra with P1 as common apex
//        // This decomposition MUST match VolumeSprings.cu exactly!
//        const double s1 = sixV_tet(X, Y, Z, p2, p3, p4, p1);
//        const double s2 = sixV_tet(X, Y, Z, p2, p4, p6, p1);
//        const double s3 = sixV_tet(X, Y, Z, p4, p5, p6, p1);
//
//        const double prism_6V = s1 + s2 + s3;
//        const double prism_vol = prism_6V / 6.0;
//
//        // Accumulate global signed volume
//        sum6V_signed += prism_6V;
//
//        // Accumulate absolute values for diagnostic
//        sum6V_abs += (std::fabs(s1) + std::fabs(s2) + std::fabs(s3));
//
//        // Per-prism energy: E_p = k_v * (V_p - V0_p)^2
//        // This matches the force law in VolumeSprings.cu exactly
//        if (have_eq_prism_volumes) {
//            const double V0p = hV0[p];
//            if (std::fabs(V0p) > 1e-12) {
//                const double dV = prism_vol - V0p;
//                per_prism_energy_sum += kv * dV * dV;
//            }
//        }
//
//        // Track minimum tet volume for inversion detection
//        const double local_min = std::min({s1, s2, s3});
//        if (local_min < min6V_tet_val) {
//            min6V_tet_val = local_min;
//            min6V_tet_prism = p;
//        }
//    }
//
//    // Convert from 6V to V
//    const double V_signed = sum6V_signed / 6.0;
//    const double V_abs = sum6V_abs / 6.0;
//
//    generalParams.current_total_volume = V_signed;
//    generalParams.true_current_total_volume = V_abs;
//
//    // Use per-prism energy (consistent with VolumeSprings.cu force law)
//    if (have_eq_prism_volumes) {
//        generalParams.volume_energy = per_prism_energy_sum;
//    } else {
//        // Fallback to global energy before equilibrium volumes are set
//        const double dV = V_signed - generalParams.eq_total_volume;
//        generalParams.volume_energy = kv * dV * dV;
//    }
//
//    // Warn if mesh is inverting
//    if (V_signed < 0) {
//        std::cout << "  WARNING: Negative total volume. "
//                  << V_signed << std::endl;
//    }
//}



// previous working global volume code. 03/09/2026
// VolumeComp.cu
// Computes total enclosed volume by summing prism volumes

#include "VolumeComp.h"
#include "System.h"
#include "SystemStructures.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

/**
 * Compute 6 times the signed volume of tetrahedron (i,j,k,l)
 * 
 * Uses scalar triple product: 6V = u · (v × w)
 * where u = r_i - r_l, v = r_j - r_l, w = r_k - r_l
 * 
 * Positive volume when vertices are in right-hand orientation.
 */
static double sixV_tet(
    const std::vector<double>& X,
    const std::vector<double>& Y,
    const std::vector<double>& Z,
    int i, int j, int k, int l)
{
    // Edge vectors relative to vertex l (apex)
    const double ux = X[i] - X[l];
    const double uy = Y[i] - Y[l];
    const double uz = Z[i] - Z[l];

    const double vx = X[j] - X[l];
    const double vy = Y[j] - Y[l];
    const double vz = Z[j] - Z[l];

    const double wx = X[k] - X[l];
    const double wy = Y[k] - Y[l];
    const double wz = Z[k] - Z[l];

    // Cross product v × w
    const double cx = vy * wz - vz * wy;
    const double cy = vz * wx - vx * wz;
    const double cz = vx * wy - vy * wx;

    // Scalar triple product u · (v × w)
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
        generalParams.current_total_volume = 0.0;
        generalParams.true_current_total_volume = 0.0;
        generalParams.volume_energy = 0.0;
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
    double sum6V_abs = 0.0;
    
    double min6V_tet = std::numeric_limits<double>::infinity();
    int min6V_tet_prism = -1;
    
    for (int p = 0; p < P; ++p) {
        const int p1 = hP1[p];
        const int p2 = hP2[p];
        const int p3 = hP3[p];
        const int p4 = hP4[p];
        const int p5 = hP5[p]; 
        const int p6 = hP6[p];
    
        // Decompose prism into 3 tetrahedra with P1 as common apex
        // This decomposition MUST match VolumeSprings.cu exactly!
        const double s1 = sixV_tet(X, Y, Z, p2, p3, p4, p1);
        const double s2 = sixV_tet(X, Y, Z, p2, p4, p6, p1);
        const double s3 = sixV_tet(X, Y, Z, p4, p5, p6, p1);
    
        // Accumulate signed volume
        sum6V_signed += (s1 + s2 + s3);
    
        // Accumulate absolute values for diagnostic
        sum6V_abs += (std::fabs(s1) + std::fabs(s2) + std::fabs(s3));
    
        // Track minimum tet volume for inversion detection
        const double local_min = std::min({s1, s2, s3});
        if (local_min < min6V_tet) {
            min6V_tet = local_min;
            min6V_tet_prism = p;
        }
    }
    
    // Convert from 6V to V
    const double V_signed = sum6V_signed / 6.0;
    const double V_abs = sum6V_abs / 6.0;
    
    generalParams.current_total_volume = V_signed;
    generalParams.true_current_total_volume = V_abs;
    
    // Compute volume energy
    const double dV = V_signed - generalParams.eq_total_volume;
    generalParams.volume_energy = generalParams.volume_spring_constant * dV * dV;
    
//    // Diagnostic output (can be commented out for production)
//    std::cout << "ComputeVolume: prisms=" << P
//              << ", V_signed=" << V_signed
//              << ", V_target=" << generalParams.eq_total_volume
//              << ", dV=" << dV
//              << ", E_vol=" << generalParams.volume_energy
//              << std::endl;
//              
    // Warn if mesh is inverting
    if (V_signed < 0) {
        std::cout << "  WARNING: Negative total volume. " 
                  << V_signed << std::endl;
    }
}