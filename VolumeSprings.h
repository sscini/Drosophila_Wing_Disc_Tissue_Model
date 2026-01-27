#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

/**
 * ComputeVolumeSprings - Apply volume-conserving forces to all nodes
 * 
 * Implements the force: F_n = -?E_vol/?r_n = -2*k_v*(O - O0)*?O/?r_n
 * 
 * The gradient ?O/?r_n is computed by summing contributions from all
 * prisms that contain node n, decomposed into tetrahedra.
 */
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs);

/**
 * VolumeSpringPrismFunctor - GPU functor for volume spring force computation
 * 
 * For each node, this functor:
 * 1. Loops over all prisms containing the node
 * 2. Decomposes each prism into 3 tetrahedra
 * 3. Accumulates gradient contributions from each tet
 * 4. Multiplies by prefactor to get force
 * 
 * Tetrahedron volume gradient derivation:
 * For V_tet = (1/6) * u · (v × w) where:
 *   u = r_i - r_l (vector from apex l to vertex i)
 *   v = r_j - r_l
 *   w = r_k - r_l
 * 
 * The gradients are:
 *   ?(6V)/?r_i = v × w
 *   ?(6V)/?r_j = w × u  
 *   ?(6V)/?r_k = u × v
 *   ?(6V)/?r_l = -[?(6V)/?r_i + ?(6V)/?r_j + ?(6V)/?r_k]
 */
struct VolumeSpringPrismFunctor
{
    double prefactor;  // = -2 * k_v * (Omega - Omega_0)
    int num_prisms;
    int num_nodes;

    // Prism connectivity arrays
    const int* P1;
    const int* P2;
    const int* P3;
    const int* P4;
    const int* P5;
    const int* P6;

    // Node coordinate arrays
    const double* X;
    const double* Y;
    const double* Z;

    __host__ __device__
    VolumeSpringPrismFunctor(
        double _prefactor,
        int _num_prisms,
        int _num_nodes,
        const int* _P1, const int* _P2, const int* _P3,
        const int* _P4, const int* _P5, const int* _P6,
        const double* _X, const double* _Y, const double* _Z)
        : prefactor(_prefactor),
          num_prisms(_num_prisms),
          num_nodes(_num_nodes),
          P1(_P1), P2(_P2), P3(_P3), P4(_P4), P5(_P5), P6(_P6),
          X(_X), Y(_Y), Z(_Z) 
    {}

    /**
     * Accumulate gradient contribution from one tetrahedron
     * 
     * Tetrahedron vertices: (i, j, k, l) where l is the apex
     * 
     * @param i,j,k,l  Vertex indices
     * @param node_id  The node we're computing gradient for
     * @param gx,gy,gz Gradient accumulator (modified in place)
     */
    __device__ inline void accumulate_tet_grad(
        int i, int j, int k, int l,
        int node_id,
        double& gx, double& gy, double& gz) const
    {
        // Bounds checking
        if (i == INT_MAX || j == INT_MAX || k == INT_MAX || l == INT_MAX) return;
        if (i < 0 || j < 0 || k < 0 || l < 0) return;
        if (i >= num_nodes || j >= num_nodes || k >= num_nodes || l >= num_nodes) return;

        // Edge vectors relative to apex l
        const double uix = X[i] - X[l];
        const double uiy = Y[i] - Y[l];
        const double uiz = Z[i] - Z[l];

        const double vjx = X[j] - X[l];
        const double vjy = Y[j] - Y[l];
        const double vjz = Z[j] - Z[l];

        const double wkx = X[k] - X[l];
        const double wky = Y[k] - Y[l];
        const double wkz = Z[k] - Z[l];

        // ?(6V)/?r_i = v × w
        const double gix = vjy * wkz - vjz * wky;
        const double giy = vjz * wkx - vjx * wkz;
        const double giz = vjx * wky - vjy * wkx;

        // ?(6V)/?r_j = w × u
        const double gjx = wky * uiz - wkz * uiy;
        const double gjy = wkz * uix - wkx * uiz;
        const double gjz = wkx * uiy - wky * uix;

        // ?(6V)/?r_k = u × v
        const double gkx = uiy * vjz - uiz * vjy;
        const double gky = uiz * vjx - uix * vjz;
        const double gkz = uix * vjy - uiy * vjx;

        // ?(6V)/?r_l = -(g_i + g_j + g_k)
        const double glx = -(gix + gjx + gkx);
        const double gly = -(giy + gjy + gky);
        const double glz = -(giz + gjz + gkz);

        // Convert from ?(6V) to ?V
        const double inv6 = 1.0 / 6.0;

        // Accumulate gradient for the relevant node
        if (node_id == i) { gx += inv6 * gix; gy += inv6 * giy; gz += inv6 * giz; }
        if (node_id == j) { gx += inv6 * gjx; gy += inv6 * gjy; gz += inv6 * gjz; }
        if (node_id == k) { gx += inv6 * gkx; gy += inv6 * gky; gz += inv6 * gkz; }
        if (node_id == l) { gx += inv6 * glx; gy += inv6 * gly; gz += inv6 * glz; }
    }

    /**
     * Main operator - compute volume force for one node
     */
    __device__
    U2CVec3 operator()(const U2CVec3& u2d3) const
    {
        const int node_id = thrust::get<0>(u2d3);
        const int bucketId = thrust::get<1>(u2d3);

        // Start with existing forces
        double fx = thrust::get<2>(u2d3);
        double fy = thrust::get<3>(u2d3);
        double fz = thrust::get<4>(u2d3);

        // Accumulate gradient from all prisms containing this node
        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int p = 0; p < num_prisms; ++p) {
            const int a = P1[p];  // Apex for all 3 tets
            const int b = P2[p];
            const int c = P3[p];
            const int A = P4[p];
            const int B = P5[p];
            const int C = P6[p];

            // Skip if this node isn't in this prism
            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // Decompose prism into 3 tets (MUST match VolumeComp.cu!)
            // Tet 1: sixV_tet(p2, p3, p4, p1) => (i=b, j=c, k=A, l=a)
            accumulate_tet_grad(b, c, A, a, node_id, gx, gy, gz);
            
            // Tet 2: sixV_tet(p2, p4, p6, p1) => (i=b, j=A, k=C, l=a)
            accumulate_tet_grad(b, A, C, a, node_id, gx, gy, gz);
            
            // Tet 3: sixV_tet(p4, p5, p6, p1) => (i=A, j=B, k=C, l=a)
            accumulate_tet_grad(A, B, C, a, node_id, gx, gy, gz);
        }

        // Apply force: F = prefactor * gradient
        // where prefactor = -2 * k_v * (Omega - Omega_0)
        fx += prefactor * gx;
        fy += prefactor * gy;
        fz += prefactor * gz;
        
        

        return thrust::make_tuple(node_id, bucketId, fx, fy, fz);
    }
};

#endif // VOLUMESPRINGS_H_