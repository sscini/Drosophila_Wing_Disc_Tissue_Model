#include "System.h"
#include "SystemStructures.h"
#include "VolumeSprings.h"

#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <cmath>
#include <iostream>

// ============================================================================
// GPU kernel: compute volume spring force for each node using atomicAdd
// ============================================================================

struct VolumeSpringAtomicFunctor
{
    double prefactor;   // = -2 * k_v * (Omega - Omega_0)
    int num_prisms;
    int num_nodes;

    // Prism connectivity (read-only)
    const int* P1;
    const int* P2;
    const int* P3;
    const int* P4;
    const int* P5;
    const int* P6;

    // Node coordinates (read-only)
    const double* X;
    const double* Y;
    const double* Z;

    // Node forces (write via atomicAdd - same as LinearSprings!)
    double* nodeForceX;
    double* nodeForceY;
    double* nodeForceZ;

    __host__ __device__
    VolumeSpringAtomicFunctor(
        double _prefactor, int _num_prisms, int _num_nodes,
        const int* _P1, const int* _P2, const int* _P3,
        const int* _P4, const int* _P5, const int* _P6,
        const double* _X, const double* _Y, const double* _Z,
        double* _fX, double* _fY, double* _fZ)
        : prefactor(_prefactor), num_prisms(_num_prisms), num_nodes(_num_nodes),
          P1(_P1), P2(_P2), P3(_P3), P4(_P4), P5(_P5), P6(_P6),
          X(_X), Y(_Y), Z(_Z),
          nodeForceX(_fX), nodeForceY(_fY), nodeForceZ(_fZ)
    {}

    // Compute gradient of tet volume w.r.t. one vertex and accumulate
    __device__ __forceinline__
    void accumulate_tet_grad(
        int i, int j, int k, int l,
        int node_id,
        double& gx, double& gy, double& gz) const
    {
        if (i < 0 || j < 0 || k < 0 || l < 0) return;
        if (i >= num_nodes || j >= num_nodes || k >= num_nodes || l >= num_nodes) return;

        const double uix = X[i] - X[l];
        const double uiy = Y[i] - Y[l];
        const double uiz = Z[i] - Z[l];

        const double vjx = X[j] - X[l];
        const double vjy = Y[j] - Y[l];
        const double vjz = Z[j] - Z[l];

        const double wkx = X[k] - X[l];
        const double wky = Y[k] - Y[l];
        const double wkz = Z[k] - Z[l];

        // d(6V)/dr_i = v x w
        const double gix = vjy * wkz - vjz * wky;
        const double giy = vjz * wkx - vjx * wkz;
        const double giz = vjx * wky - vjy * wkx;

        // d(6V)/dr_j = w x u
        const double gjx = wky * uiz - wkz * uiy;
        const double gjy = wkz * uix - wkx * uiz;
        const double gjz = wkx * uiy - wky * uix;

        // d(6V)/dr_k = u x v
        const double gkx = uiy * vjz - uiz * vjy;
        const double gky = uiz * vjx - uix * vjz;
        const double gkz = uix * vjy - uiy * vjx;

        // d(6V)/dr_l = -(g_i + g_j + g_k)
        const double glx = -(gix + gjx + gkx);
        const double gly = -(giy + gjy + gky);
        const double glz = -(giz + gjz + gkz);

        const double inv6 = 1.0 / 6.0;

        if (node_id == i) { gx += inv6 * gix; gy += inv6 * giy; gz += inv6 * giz; }
        if (node_id == j) { gx += inv6 * gjx; gy += inv6 * gjy; gz += inv6 * gjz; }
        if (node_id == k) { gx += inv6 * gkx; gy += inv6 * gky; gz += inv6 * gkz; }
        if (node_id == l) { gx += inv6 * glx; gy += inv6 * gly; gz += inv6 * glz; }
    }

    __device__
    void operator()(int node_id) const
    {
        if (node_id >= num_nodes) return;

        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int p = 0; p < num_prisms; ++p) {
            const int a = P1[p];
            const int b = P2[p];
            const int c = P3[p];
            const int A = P4[p];
            const int B = P5[p];
            const int C = P6[p];

            // Skip if this node isn't in this prism
            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // 3-tet decomposition (matches VolumeComp.cu exactly)
            accumulate_tet_grad(b, c, A, a, node_id, gx, gy, gz);
            accumulate_tet_grad(b, A, C, a, node_id, gx, gy, gz);
            accumulate_tet_grad(A, B, C, a, node_id, gx, gy, gz);
        }

        // Force = prefactor * gradient
        double fx = prefactor * gx;
        double fy = prefactor * gy;
        double fz = prefactor * gz;

        // CRITICAL FIX: Use atomicAdd (same as LinearSprings!)
        // This ACCUMULATES onto existing forces instead of overwriting them
        if (fx != 0.0 || fy != 0.0 || fz != 0.0) {
            atomicAdd(&nodeForceX[node_id], fx);
            atomicAdd(&nodeForceY[node_id], fy);
            atomicAdd(&nodeForceZ[node_id], fz);
        }
    }
};

// ============================================================================
// ComputeVolumeSprings - main entry point
// ============================================================================

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs)
{
    const int P = prismInfoVecs.num_prisms; 
    if (P <= 0) {
        return;
    }

    const double kv = generalParams.volume_spring_constant;
    if (kv == 0.0) {
        return;
    }

    const double Omega_s = generalParams.current_total_volume;
    const double Omega0 = generalParams.eq_total_volume;
    
    // Safeguard: NaN/Inf check
    if (std::isnan(Omega_s) || std::isinf(Omega_s)) {
        std::cout << "WARNING: Volume is NaN/Inf (" << Omega_s << "), skipping volume springs." << std::endl;
        return;
    }
    
    // Safeguard: Negative volume (mesh inversion)
    if (Omega_s < 0.0) {
        std::cout << "WARNING: Negative volume (" << Omega_s 
                  << "). Skipping volume springs." << std::endl;
        return;
    }
    
    double volume_diff = Omega_s - Omega0;
    double volume_ratio = Omega_s / Omega0;
    double prefactor;
    
    // Clamp extreme volume changes
    if (volume_ratio < 0.5 || volume_ratio > 2.0) {
        std::cout << "WARNING: Extreme volume change (ratio=" << volume_ratio 
                  << "). Clamping." << std::endl;
        double max_diff = 0.5 * std::fabs(Omega0);
        double clamped_diff = std::max(-max_diff, std::min(max_diff, volume_diff));
        prefactor = -2.0 * kv * clamped_diff;
    } else {
        prefactor = -2.0 * kv * volume_diff;
    }
    
    // Safeguard: NaN/Inf prefactor
    if (std::isnan(prefactor) || std::isinf(prefactor)) {
        std::cout << "WARNING: Volume spring prefactor is NaN/Inf. Skipping." << std::endl;
        return;
    }
    
    // Clamp magnitude
    const double max_prefactor = 1e6;
    if (std::fabs(prefactor) > max_prefactor) {
        prefactor = (prefactor > 0) ? max_prefactor : -max_prefactor;
    }
    
    generalParams.volume_energy = kv * volume_diff * volume_diff;
    const int Nnodes = (int)coordInfoVecs.nodeLocX.size();

    // Use for_each with atomicAdd (NOT thrust::transform!)
    VolumeSpringAtomicFunctor functor(
        prefactor, P, Nnodes,
        thrust::raw_pointer_cast(prismInfoVecs.P1.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P2.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P3.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P4.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P5.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P6.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data()));

    // for_each calls functor(node_id) for each node
    // functor uses atomicAdd to accumulate forces
    thrust::for_each(
        thrust::device,
        thrust::counting_iterator<int>(0),
        thrust::counting_iterator<int>(Nnodes),
        functor);
}
