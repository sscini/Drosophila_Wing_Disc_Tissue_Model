#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"
#include <stdio.h>

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs);

// LinearSpringFunctor — branchless per-edge spring constant lookup + atomicAdd
struct LinearSpringFunctor {
    double* edge_spring_k;
    double* rest_length;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;
    double* nodeForceX;
    double* nodeForceY;
    double* nodeForceZ;
    int maxNodeCount;

    __host__ __device__ LinearSpringFunctor(
        double* _edge_spring_k, double* _rest_length,
        double* _locXAddr, double* _locYAddr, double* _locZAddr,
        double* _nodeForceX, double* _nodeForceY, double* _nodeForceZ,
        int _maxNodeCount) :
        edge_spring_k(_edge_spring_k), rest_length(_rest_length),
        locXAddr(_locXAddr), locYAddr(_locYAddr), locZAddr(_locZAddr),
        nodeForceX(_nodeForceX), nodeForceY(_nodeForceY), nodeForceZ(_nodeForceZ),
        maxNodeCount(_maxNodeCount) {}

    __device__ double operator()(const Tuuu &u3d) {
        int counter = thrust::get<0>(u3d);
        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        if (edgeL == INT_MAX || edgeL < 0 || edgeL >= maxNodeCount ||
            edgeR == INT_MAX || edgeR < 0 || edgeR >= maxNodeCount)
            return 0.0;

        double length_zero = rest_length[counter];
        if (length_zero <= 0.0 || isnan(length_zero)) return 0.0;

        double k = edge_spring_k[counter];

        double xLR = locXAddr[edgeL] - locXAddr[edgeR];
        double yLR = locYAddr[edgeL] - locYAddr[edgeR];
        double zLR = locZAddr[edgeL] - locZAddr[edgeR];
        double len_sq = xLR*xLR + yLR*yLR + zLR*zLR;
        if (len_sq < 1e-20) return 0.0;

        double len_cur = sqrt(len_sq);
        double len_diff = len_cur - length_zero;
        if (fabs(len_diff) < 1e-15) return 0.0;

        double mag = -k * len_diff;
        double inv = 1.0 / len_cur;
        double fx = mag * xLR * inv;
        double fy = mag * yLR * inv;
        double fz = mag * zLR * inv;

        atomicAdd(&nodeForceX[edgeL],  fx);
        atomicAdd(&nodeForceY[edgeL],  fy);
        atomicAdd(&nodeForceZ[edgeL],  fz);
        atomicAdd(&nodeForceX[edgeR], -fx);
        atomicAdd(&nodeForceY[edgeR], -fy);
        atomicAdd(&nodeForceZ[edgeR], -fz);

        return 0.5 * k * len_diff * len_diff;
    }
};

#endif /* LINEARSPRINGS_H_ */