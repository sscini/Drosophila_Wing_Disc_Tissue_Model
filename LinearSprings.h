#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"
#include <stdio.h>

// Declare the function ComputeLinearSprings, which is responsible for calculating linear spring forces.
void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs);

// LinearSpringFunctor using atomicAdd for direct force accumulation
// This eliminates the need for sort_by_key and reduces synchronization issues
struct LinearSpringFunctor {
    double spring_constant;           // Normal spring constant for body/basal edges
    double spring_constant_weak;      // Weak spring constant for apical (layer 1) edges
    double spring_constant_vertical;  // Spring constant for vertical edges (layer -1)
    double* rest_length;              // Rest length array for each edge
    int* edges_in_upperhem;           // Edge layer classification
    
    // Node positions (read-only)
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;
    
    // Node forces (write via atomicAdd)
    double* nodeForceX;
    double* nodeForceY;
    double* nodeForceZ;
    
    int maxNodeCount;

    // Constructor
    __host__ __device__ LinearSpringFunctor(
        double _spring_constant,
        double _spring_constant_weak,
        double _spring_constant_vertical,
        double* _rest_length,
        int* _edges_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        double* _nodeForceX,
        double* _nodeForceY,
        double* _nodeForceZ,
        int _maxNodeCount) :
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        spring_constant_vertical(_spring_constant_vertical),
        rest_length(_rest_length),
        edges_in_upperhem(_edges_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        nodeForceX(_nodeForceX),
        nodeForceY(_nodeForceY),
        nodeForceZ(_nodeForceZ),
        maxNodeCount(_maxNodeCount) {}

    __device__ double operator()(const Tuuu &u3d) {
        // counter ranges from 0 to num_edges
        int counter = thrust::get<0>(u3d);
        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        // Validate node indices
        if (edgeL == INT_MAX || edgeL < 0 || edgeL >= maxNodeCount ||
            edgeR == INT_MAX || edgeR < 0 || edgeR >= maxNodeCount) {
            return 0.0;
        }

        double length_zero = rest_length[counter];
        
        // Check for valid rest length
        if (length_zero <= 0.0 || isnan(length_zero)) {
            return 0.0;
        }
        
        // Determine spring constant based on edge layer
        double what_spring_constant;// = spring_constant;
        int edge_layer = edges_in_upperhem[counter];
        
        if (edge_layer >= 0) {
            // Horizontal edges in upper hemisphere - weak springs
            what_spring_constant = spring_constant;
        }
        else if (edge_layer == -1) {
            // Vertical edges - use vertical spring constant
            what_spring_constant = spring_constant_vertical;
        }
        else {
            // All other edges - normal spring constant
            what_spring_constant = spring_constant;
        }

        // Compute displacement vector
        double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
        double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
        double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];

        // Compute current edge length
        double length_sq = xLoc_LR * xLoc_LR + yLoc_LR * yLoc_LR + zLoc_LR * zLoc_LR;
        
        // Avoid division by zero for degenerate edges
        if (length_sq < 1e-20) {
            return 0.0;
        }
        
        double length_current = sqrt(length_sq);
        double length_diff = length_current - length_zero;

        // Only compute forces if there's actual deformation
        if (fabs(length_diff) < 1e-15) {
            return 0.0;
        }

        // Compute force magnitude: F = -k * (L - L0)
        // Negative sign because force opposes extension
        double magnitude = -what_spring_constant * length_diff;

        // Compute unit vector from edgeR to edgeL
        double inv_length = 1.0 / length_current;
        double ux = xLoc_LR * inv_length;
        double uy = yLoc_LR * inv_length;
        double uz = zLoc_LR * inv_length;

        // Force components on edgeL node
        double fx = magnitude * ux;
        double fy = magnitude * uy;
        double fz = magnitude * uz;

        // Accumulate forces using atomicAdd (thread-safe, no sort needed!)
        // Force on edgeL
        atomicAdd(&nodeForceX[edgeL], fx);
        atomicAdd(&nodeForceY[edgeL], fy);
        atomicAdd(&nodeForceZ[edgeL], fz);

        // Equal and opposite force on edgeR (Newton's 3rd law)
        atomicAdd(&nodeForceX[edgeR], -fx);
        atomicAdd(&nodeForceY[edgeR], -fy);
        atomicAdd(&nodeForceZ[edgeR], -fz);

        // Compute and return energy: E = 0.5 * k * (L - L0)^2
        double energy = 0.5 * what_spring_constant * length_diff * length_diff;
        
        return energy;
    }
};

#endif /* LINEARSPRINGS_H_ */