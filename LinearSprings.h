#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"
#include <stdio.h>
#include <atomic>

// Declare the function ComputeLinearSprings, which is responsible for calculating linear spring forces.
void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs);

// In LinearSprings.h, modify the functor to write forces directly:

struct LinearSpringDirectFunctor {
    double spring_constant;
    double* rest_length;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;
    double* nodeForceX;  // Direct output
    double* nodeForceY;
    double* nodeForceZ;
    int maxNodeCount;

    __host__ __device__
    LinearSpringDirectFunctor(
        double _k,
        double* _rest,
        double* _x, double* _y, double* _z,
        double* _fx, double* _fy, double* _fz,
        int _max)
        : spring_constant(_k), rest_length(_rest),
          locXAddr(_x), locYAddr(_y), locZAddr(_z),
          nodeForceX(_fx), nodeForceY(_fy), nodeForceZ(_fz),
          maxNodeCount(_max) {}

    __device__ double operator()(const Tuuu &u3d) {
        int counter = thrust::get<0>(u3d);
        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        // Validate indices
        if (edgeL < 0 || edgeL >= maxNodeCount ||
            edgeR < 0 || edgeR >= maxNodeCount ||
            edgeL == INT_MAX || edgeR == INT_MAX) {
            return 0.0;
        }

        double length_zero = rest_length[counter];
        if (length_zero <= 0.0 || isnan(length_zero)) return 0.0;

        // Compute displacement
        double dx = locXAddr[edgeL] - locXAddr[edgeR];
        double dy = locYAddr[edgeL] - locYAddr[edgeR];
        double dz = locZAddr[edgeL] - locZAddr[edgeR];

        double length_sq = dx*dx + dy*dy + dz*dz;
        if (length_sq < 1e-20) return 0.0;

        double length_current = sqrt(length_sq);
        double length_diff = length_current - length_zero;

        // Force magnitude
        double magnitude = -spring_constant * length_diff;

        // Unit vector
        double inv_len = 1.0 / length_current;
        double fx = magnitude * dx * inv_len;
        double fy = magnitude * dy * inv_len;
        double fz = magnitude * dz * inv_len;

        // Accumulate forces using atomicAdd (no sort needed!)
        atomicAdd(&nodeForceX[edgeL], fx);
        atomicAdd(&nodeForceY[edgeL], fy);
        atomicAdd(&nodeForceZ[edgeL], fz);

        atomicAdd(&nodeForceX[edgeR], -fx);
        atomicAdd(&nodeForceY[edgeR], -fy);
        atomicAdd(&nodeForceZ[edgeR], -fz);

        // Return energy
        return 0.5 * spring_constant * length_diff * length_diff;
    }
};

#endif