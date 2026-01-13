#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    //AuxVecs& auxVecs,
    PrismInfoVecs& prismInfoVecs);

struct VolumeSpringPrismFunctor
    : public thrust::unary_function<U2CVec3, U2CVec3>
{
    double prefactor;
    int num_prisms;
    int num_nodes;

    const int* P1;
    const int* P2;
    const int* P3;
    const int* P4;
    const int* P5;
    const int* P6;

    const double* X;
    const double* Y;
    const double* Z;

    __host__ __device__
    VolumeSpringPrismFunctor(double _prefactor,
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

    __device__ inline void accumulate_tet_grad(int i, int j, int k, int l,
                                               int node_id,
                                               double& gx, double& gy, double& gz) const
    {
        // basic validity checks
        if (i == INT_MAX || j == INT_MAX || k == INT_MAX || l == INT_MAX) return;
        if (i < 0 || j < 0 || k < 0 || l < 0) return;
        if (i >= num_nodes || j >= num_nodes || k >= num_nodes || l >= num_nodes) return;
    
        // vectors relative to l (THIS matches sixV_tet)
        const double uix = X[i] - X[l], uiy = Y[i] - Y[l], uiz = Z[i] - Z[l];
        const double vjx = X[j] - X[l], vjy = Y[j] - Y[l], vjz = Z[j] - Z[l];
        const double wkx = X[k] - X[l], wky = Y[k] - Y[l], wkz = Z[k] - Z[l];
    
        // g_i = cross(v, w)
        const double gix = vjy * wkz - vjz * wky;
        const double giy = vjz * wkx - vjx * wkz;
        const double giz = vjx * wky - vjy * wkx;
    
        // g_j = cross(w, u)
        const double gjx = wky * uiz - wkz * uiy;
        const double gjy = wkz * uix - wkx * uiz;
        const double gjz = wkx * uiy - wky * uix;
    
        // g_k = cross(u, v)
        const double gkx = uiy * vjz - uiz * vjy;
        const double gky = uiz * vjx - uix * vjz;
        const double gkz = uix * vjy - uiy * vjx;
    
        // g_l = -(g_i + g_j + g_k)
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
    U2CVec3 operator()(const U2CVec3& u2d3) const
    {
        const int node_id  = thrust::get<0>(u2d3);
        const int bucketId = thrust::get<1>(u2d3);

        double fx = thrust::get<2>(u2d3);
        double fy = thrust::get<3>(u2d3);
        double fz = thrust::get<4>(u2d3);

        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int p = 0; p < num_prisms; ++p) {
            const int a = P1[p], b = P2[p], c = P3[p];
            const int A = P4[p], B = P5[p], C = P6[p];

            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // SAME tet decomposition pattern you used (just using explicit A,B,C):
            accumulate_tet_grad(b, c, A, a, node_id, gx, gy, gz);
            accumulate_tet_grad(b, A, C, a, node_id, gx, gy, gz);
            accumulate_tet_grad(A, B, C, a, node_id, gx, gy, gz);
        }

        fx += prefactor * gx;
        fy += prefactor * gy;
        fz += prefactor * gz;

        return thrust::make_tuple(node_id, bucketId, fx, fy, fz);
    }
};

#endif
