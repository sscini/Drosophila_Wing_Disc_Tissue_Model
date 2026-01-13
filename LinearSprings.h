#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"
#include <stdio.h>

// Declare the function ComputeLinearSprings, which is responsible for calculating linear spring forces.
void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);

// Define a functor named LinearSpringFunctor for calculating linear spring forces
struct LinearSpringFunctor {
    // Member variables representing parameters required by the functor
    int SCALE_TYPE;
    bool nonuniform_wall_weakening_linear;
    double maxSpringScaler_linear;
    double scaling_pow;
    double gausssigma;
    double hilleqnconst;
    double hilleqnpow;
    double* scaling_per_edge;
    double spring_constant;
    double spring_constant_weak;
    double spring_constant_vertical;
    double* rest_length;
    int* edges_in_upperhem;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
    // Store max node count for bounds checking
    int maxNodeCount;
    
    // Constructor for the LinearSpringFunctor
    __host__ __device__ LinearSpringFunctor(
        int& _SCALE_TYPE,
        bool& _nonuniform_wall_weakening_linear,
        double& _maxSpringScaler_linear,
        double& _scaling_pow,
        double& _gausssigma,
        double& _hilleqnconst,
        double& _hilleqnpow,
        double* _scaling_per_edge,
        double& _spring_constant,
        double& _spring_constant_weak,
        double& _spring_constant_vertical,
        double* _rest_length,
        int* _edges_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr,
        int _maxNodeCount = 0) :  // Added parameter with default
        SCALE_TYPE(_SCALE_TYPE),
        nonuniform_wall_weakening_linear(_nonuniform_wall_weakening_linear),
        maxSpringScaler_linear(_maxSpringScaler_linear),
        scaling_pow(_scaling_pow),
        gausssigma(_gausssigma),
        hilleqnconst(_hilleqnconst),
        hilleqnpow(_hilleqnpow),
        scaling_per_edge(_scaling_per_edge),
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        spring_constant_vertical(_spring_constant_vertical),
        rest_length(_rest_length),
        edges_in_upperhem(_edges_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr),
        maxNodeCount(_maxNodeCount) {}

    // Functor operator for the LinearSpringFunctor
    __device__ double operator()(const Tuuu &u3d) {
        
        // Counter ranges from 0 to num_edges
        int counter = thrust::get<0>(u3d);
        int place = 2 * counter;  // Location in write-to vector

        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        // ========================================================================
        // CRITICAL FIX: Always initialize output arrays to prevent NaN propagation
        // Use node index 0 as a safe default (forces will be zero, so no effect)
        // ========================================================================
        idKey[place] = 0;
        idKey[place + 1] = 0;
        forceXAddr[place] = 0.0;
        forceYAddr[place] = 0.0;
        forceZAddr[place] = 0.0;
        forceXAddr[place + 1] = 0.0;
        forceYAddr[place + 1] = 0.0;
        forceZAddr[place + 1] = 0.0;

        // Validate edge node indices
        // Check for INT_MAX, negative values, and out-of-bounds
        bool validEdge = (edgeL != INT_MAX && edgeL >= 0 && 
                          edgeR != INT_MAX && edgeR >= 0);
        
        // Additional bounds check if maxNodeCount is set
        if (validEdge && maxNodeCount > 0) {
            validEdge = (edgeL < maxNodeCount && edgeR < maxNodeCount);
        }

        if (!validEdge) {
            // Edge is invalid - return zero energy (arrays already initialized above)
            // Uncomment for debugging:
            // printf("SKIPPING invalid edge %d: edgeL=%d, edgeR=%d\n", counter, edgeL, edgeR);
            return 0.0;
        }

        // Get rest length for this edge
        double length_zero = rest_length[counter];
        
        // Validate rest length
        if (length_zero <= 0.0 || isnan(length_zero) || isinf(length_zero)) {
            // Invalid rest length - skip this edge
            // printf("SKIPPING edge %d: invalid rest_length=%f\n", counter, length_zero);
            return 0.0;
        }

        // Determine spring constant based on edge type
        double what_spring_constant = spring_constant;
        
        // You can uncomment and modify this section if you need layer-specific constants:
        /*
        int edgeLayer = edges_in_upperhem[counter];
        if (edgeLayer == -1) {
            // Vertical edge
            what_spring_constant = spring_constant_vertical;
        } else if (edgeLayer >= 0) {
            // Horizontal edge (apical or basal layer)
            what_spring_constant = spring_constant;
        }
        */

        // Compute displacement vector
        double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
        double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
        double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];

        // Compute current edge length
        double length_squared = (xLoc_LR * xLoc_LR) + 
                                (yLoc_LR * yLoc_LR) + 
                                (zLoc_LR * zLoc_LR);
        
        // Check for degenerate edge (zero length)
        if (length_squared < 1e-20) {
            // Degenerate edge - nodes are coincident
            // printf("WARNING: Degenerate edge %d, length^2=%e\n", counter, length_squared);
            return 0.0;
        }

        double length_current = sqrt(length_squared);

        // Validate computed length
        if (isnan(length_current) || isinf(length_current)) {
            // printf("ERROR: Invalid length_current for edge %d\n", counter);
            return 0.0;
        }

        double energy = 0.0;

        // Only compute forces if there's a length difference
        double length_diff = length_current - length_zero;
        
        if (fabs(length_diff) > 1e-15) {
            // Compute force magnitude: F = -k * (L - L0)
            // Negative sign because force opposes extension
            double magnitude = -what_spring_constant * length_diff;

            // Validate magnitude
            if (isnan(magnitude) || isinf(magnitude)) {
                // printf("ERROR: Invalid magnitude for edge %d\n", counter);
                return 0.0;
            }

            // Compute unit vector components
            double inv_length = 1.0 / length_current;
            double ux = xLoc_LR * inv_length;
            double uy = yLoc_LR * inv_length;
            double uz = zLoc_LR * inv_length;

            // Assign forces to left node (edgeL)
            idKey[place] = edgeL;
            forceXAddr[place] = magnitude * ux;
            forceYAddr[place] = magnitude * uy;
            forceZAddr[place] = magnitude * uz;

            // Assign equal and opposite forces to right node (edgeR)
            idKey[place + 1] = edgeR;
            forceXAddr[place + 1] = -magnitude * ux;
            forceYAddr[place + 1] = -magnitude * uy;
            forceZAddr[place + 1] = -magnitude * uz;

            // Compute spring potential energy: E = 0.5 * k * (L - L0)^2
            energy = 0.5 * what_spring_constant * length_diff * length_diff;
        }
        else {
            // Edge is at rest length - set proper node IDs but zero forces
            idKey[place] = edgeL;
            idKey[place + 1] = edgeR;
            // Forces already set to 0.0 above
        }

        return energy;
    }
};

#endif