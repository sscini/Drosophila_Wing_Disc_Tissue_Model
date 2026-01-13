#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

/**
 * VOLUME CONSERVATION SPRINGS - CORRECTED VERSION
 * 
 * Mathematical Background:
 * ========================
 * 
 * For a tetrahedron with vertices (r_i, r_j, r_k, r_l), we define:
 *   u = r_i - r_l
 *   v = r_j - r_l
 *   w = r_k - r_l
 * 
 * The signed volume is:
 *   6V = u · (v × w) = det([u, v, w])
 * 
 * The gradients with respect to each vertex position are:
 *   ?(6V)/?r_i = v × w
 *   ?(6V)/?r_j = w × u  (equivalently: -(u × w))
 *   ?(6V)/?r_k = u × v  (equivalently: -(v × u))
 *   ?(6V)/?r_l = -[?(6V)/?r_i + ?(6V)/?r_j + ?(6V)/?r_k]
 * 
 * CRITICAL: The order in cross products matters!
 *   v × w ? w × v (they differ by a sign)
 * 
 * Energy and Force:
 * =================
 * E_vol = k_v * (O - O_0)²
 * 
 * F_n = -?E_vol/?r_n = -2 * k_v * (O - O_0) * ?O/?r_n
 * 
 * Where O is the total enclosed volume (sum over all prisms).
 */

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs);

/**
 * CORRECTED Volume Spring Functor
 * 
 * This functor computes volume-preserving forces on each node by summing
 * contributions from all prisms that contain that node.
 */
struct VolumeSpringPrismFunctor
{
    double prefactor;  // = -2 * k_v * (Omega - Omega_0)
    int num_prisms;
    int num_nodes;

    // Prism connectivity: each prism has 6 vertices (P1-P6)
    // Convention: P1,P2,P3 form bottom triangle, P4,P5,P6 form top triangle
    const int* P1;
    const int* P2;
    const int* P3;
    const int* P4;
    const int* P5;
    const int* P6;

    // Node coordinates
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
     * CORRECTED gradient accumulation for a single tetrahedron
     * 
     * For tetrahedron (i, j, k, l) where l is the "apex" vertex:
     *   u = r_i - r_l
     *   v = r_j - r_l
     *   w = r_k - r_l
     *   6V = u · (v × w)
     * 
     * Gradients:
     *   ?(6V)/?r_i = v × w
     *   ?(6V)/?r_j = w × u
     *   ?(6V)/?r_k = u × v
     *   ?(6V)/?r_l = -(g_i + g_j + g_k)
     */
    __device__ inline void accumulate_tet_grad_corrected(
        int i, int j, int k, int l,  // Vertex indices for tet (i,j,k,l)
        int node_id,                  // Which node we're computing gradient for
        double& gx, double& gy, double& gz) const
    {
        // Bounds checking
        if (i == INT_MAX || j == INT_MAX || k == INT_MAX || l == INT_MAX) return;
        if (i < 0 || j < 0 || k < 0 || l < 0) return;
        if (i >= num_nodes || j >= num_nodes || k >= num_nodes || l >= num_nodes) return;

        // Compute edge vectors relative to vertex l (the apex)
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

        // ?(6V)/?r_j = w × u  (CORRECTED: was incorrectly u × w in original)
        // w × u = (wky*uiz - wkz*uiy, wkz*uix - wkx*uiz, wkx*uiy - wky*uix)
        const double gjx = wky * uiz - wkz * uiy;
        const double gjy = wkz * uix - wkx * uiz;
        const double gjz = wkx * uiy - wky * uix;

        // ?(6V)/?r_k = u × v  (CORRECTED: was incorrectly v × u in original)
        // u × v = (uiy*vjz - uiz*vjy, uiz*vjx - uix*vjz, uix*vjy - uiy*vjx)
        const double gkx = uiy * vjz - uiz * vjy;
        const double gky = uiz * vjx - uix * vjz;
        const double gkz = uix * vjy - uiy * vjx;

        // ?(6V)/?r_l = -(g_i + g_j + g_k)
        const double glx = -(gix + gjx + gkx);
        const double gly = -(giy + gjy + gky);
        const double glz = -(giz + gjz + gkz);

        // Convert from ?(6V) to ?V by dividing by 6
        const double inv6 = 1.0 / 6.0;

        // Accumulate gradient for the specific node
        if (node_id == i) { gx += inv6 * gix; gy += inv6 * giy; gz += inv6 * giz; }
        if (node_id == j) { gx += inv6 * gjx; gy += inv6 * gjy; gz += inv6 * gjz; }
        if (node_id == k) { gx += inv6 * gkx; gy += inv6 * gky; gz += inv6 * gkz; }
        if (node_id == l) { gx += inv6 * glx; gy += inv6 * gly; gz += inv6 * glz; }
    }

    /**
     * Main operator: compute volume-preserving force contribution for a single node
     * 
     * For each prism containing this node, we decompose it into 3 tetrahedra
     * and accumulate the gradients.
     * 
     * Prism decomposition (using P1 as common apex):
     *   Tet 1: (P2, P3, P4, P1) - connects bottom edge P2-P3 to top vertex P4
     *   Tet 2: (P2, P4, P6, P1) - connects P2 to diagonal P4-P6
     *   Tet 3: (P4, P5, P6, P1) - top triangle
     * 
     * This decomposition must match exactly what's used in VolumeComp.cu!
     */
    __device__
    U2CVec3 operator()(const U2CVec3& u2d3) const
    {
        const int node_id  = thrust::get<0>(u2d3);
        const int bucketId = thrust::get<1>(u2d3);

        // Start with existing forces
        double fx = thrust::get<2>(u2d3);
        double fy = thrust::get<3>(u2d3);
        double fz = thrust::get<4>(u2d3);

        // Accumulate gradient contributions from all prisms containing this node
        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int p = 0; p < num_prisms; ++p) {
            const int a = P1[p], b = P2[p], c = P3[p];
            const int A = P4[p], B = P5[p], C = P6[p];

            // Skip if this node isn't part of this prism
            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // Decompose prism into 3 tetrahedra (same as VolumeComp.cu)
            // Tet 1: sixV_tet(b, c, A, a) => tet with vertices (i=b, j=c, k=A, l=a)
            accumulate_tet_grad_corrected(b, c, A, a, node_id, gx, gy, gz);
            
            // Tet 2: sixV_tet(b, A, C, a) => tet with vertices (i=b, j=A, k=C, l=a)
            accumulate_tet_grad_corrected(b, A, C, a, node_id, gx, gy, gz);
            
            // Tet 3: sixV_tet(A, B, C, a) => tet with vertices (i=A, j=B, k=C, l=a)
            accumulate_tet_grad_corrected(A, B, C, a, node_id, gx, gy, gz);
        }

        // Apply prefactor: F = -2 * k_v * (Omega - Omega_0) * grad(Omega)
        fx += prefactor * gx;
        fy += prefactor * gy;
        fz += prefactor * gz;

        return thrust::make_tuple(node_id, bucketId, fx, fy, fz);
    }
};

#endif // VOLUMESPRINGS_H_
