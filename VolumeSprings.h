#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

// Compute global volume *forces* using your prism-based volume O
// and the force formula: F_n = -2 k_v (O_abs - O0) sign(O_s) ?O/?r_n.
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

// Functor applied per node: takes (node_id, dummy, Fx, Fy, Fz),
// adds prism-volume forces, and returns the same 5-tuple.
struct VolumeSpringPrismFunctor
    : public thrust::unary_function<U2CVec3, U2CVec3> {

    double prefactor;        // -2 k_v (O_abs - O0) sign(O_s)
    int nodes_per_layer;     // N
    int tris_per_layer;      // Tper
    int num_layers;          // L
    int num_prisms;          // (L-1)*Tper

    const int* tri1;
    const int* tri2;
    const int* tri3;
    const double* X;
    const double* Y;
    const double* Z;

    __host__ __device__
    VolumeSpringPrismFunctor(double _prefactor,
                             int _nodes_per_layer,
                             int _tris_per_layer,
                             int _num_layers,
                             int _num_prisms,
                             const int* _tri1,
                             const int* _tri2,
                             const int* _tri3,
                             const double* _X,
                             const double* _Y,
                             const double* _Z)
        : prefactor(_prefactor),
          nodes_per_layer(_nodes_per_layer),
          tris_per_layer(_tris_per_layer),
          num_layers(_num_layers),
          num_prisms(_num_prisms),
          tri1(_tri1), tri2(_tri2), tri3(_tri3),
          X(_X), Y(_Y), Z(_Z) {}

    // Helper: accumulate gradient contribution from one tetrahedron (i,j,k,l)
    // to node 'node_id'. We work with "6V" (det) gradients, then divide by 6
    // to get ?O_tet/?r_n.
    __device__ inline void accumulate_tet_grad(int i, int j, int k, int l,
                                               int node_id,
                                               double& gx,
                                               double& gy,
                                               double& gz) const {

        if (i == INT_MAX || j == INT_MAX || k == INT_MAX || l == INT_MAX)
            return;

        // Positions
        const double rix = X[i], riy = Y[i], riz = Z[i];
        const double rjx = X[j], rjy = Y[j], rjz = Z[j];
        const double rkx = X[k], rky = Y[k], rkz = Z[k];
        const double rlx = X[l], rly = Y[l], rlz = Z[l];

        // Differences with respect to i
        const double bkx = rkx - rix;
        const double bky = rky - riy;
        const double bkz = rkz - riz;

        const double ckx = rlx - rix;
        const double cky = rly - riy;
        const double ckz = rlz - riz;

        const double ajx = rjx - rix;
        const double ajy = rjy - riy;
        const double ajz = rjz - riz;

        // g_j = ?(6V)/?r_j = (rk - ri) × (rl - ri)
        const double gjx = bky * ckz - bkz * cky;
        const double gjy = bkz * ckx - bkx * ckz;
        const double gjz = bkx * cky - bky * ckx;

        // g_k = (rl - ri) × (rj - ri)
        const double gkx = cky * ajz - ckz * ajy;
        const double gky = ckz * ajx - ckx * ajz;
        const double gkz = ckx * ajy - cky * ajx;

        // g_l = (rj - ri) × (rk - ri)
        const double glx = ajy * bkz - ajz * bky;
        const double gly = ajz * bkx - ajx * bkz;
        const double glz = ajx * bky - ajy * bkx;

        // g_i = - (g_j + g_k + g_l)
        const double gix = -(gjx + gkx + glx);
        const double giy = -(gjy + gky + gly);
        const double giz = -(gjz + gkz + glz);

        // Convert 6V gradient to O_tet gradient: divide by 6.
        const double inv6 = 1.0 / 6.0;

        if (node_id == i) {
            gx += inv6 * gix;
            gy += inv6 * giy;
            gz += inv6 * giz;
        }
        if (node_id == j) {
            gx += inv6 * gjx;
            gy += inv6 * gjy;
            gz += inv6 * gjz;
        }
        if (node_id == k) {
            gx += inv6 * gkx;
            gy += inv6 * gky;
            gz += inv6 * gkz;
        }
        if (node_id == l) {
            gx += inv6 * glx;
            gy += inv6 * gly;
            gz += inv6 * glz;
        }
    }

    __device__
    U2CVec3 operator()(const U2CVec3& u2d3) const {

        const int node_id  = thrust::get<0>(u2d3);
        const int bucketId = thrust::get<1>(u2d3);  // not used, but preserved

        double fx = thrust::get<2>(u2d3);
        double fy = thrust::get<3>(u2d3);
        double fz = thrust::get<4>(u2d3);

        // Accumulate ?O/?r_n over all prisms touching node_id
        double gx = 0.0;
        double gy = 0.0;
        double gz = 0.0;

        // Loop over all prisms (same indexing as in ComputeVolume)
        for (int p = 0; p < num_prisms; ++p) {
            const int ell = p / tris_per_layer;  // layer index
            const int t   = p % tris_per_layer;  // triangle index within layer
            const int tri_idx = ell * tris_per_layer + t;

            const int a = tri1[tri_idx];
            const int b = tri2[tri_idx];
            const int c = tri3[tri_idx];

            // top nodes
            const int A = (a == INT_MAX) ? INT_MAX : (a + nodes_per_layer);
            const int B = (b == INT_MAX) ? INT_MAX : (b + nodes_per_layer);
            const int C = (c == INT_MAX) ? INT_MAX : (c + nodes_per_layer);

            if (a == INT_MAX || b == INT_MAX || c == INT_MAX)
                continue;

            // If node_id is nowhere in this prism, skip quickly
            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // Prism decomposition: 3 tets
            // T1: (P1,P2,P3,P4) = (a,b,c,A)
            accumulate_tet_grad(a, b, c, A, node_id, gx, gy, gz);
            // T2: (P2,P3,P4,P6) = (b,c,A,C)
            accumulate_tet_grad(b, c, A, C, node_id, gx, gy, gz);
            // T3: (P4,P5,P6,P1) = (A,B,C,a)
            accumulate_tet_grad(A, B, C, a, node_id, gx, gy, gz);
        }

        // Apply global prefactor from energy derivative
        fx += prefactor * gx;
        fy += prefactor * gy;
        fz += prefactor * gz;

        // Return full tuple so zip_iterator types match
        return thrust::make_tuple(node_id, bucketId, fx, fy, fz);
    }
};

#endif // VOLUMESPRINGS_H_
