#include "gradientRelax.h"
#include <vector>            // for std::vector
#include <thrust/copy.h>     // for thrust::copy
#include <limits>            // for std::numeric_limits
#include <cmath>             // for std::sqrt
#include <iostream>          // for std::cout

int relaxUntilConverged(System& system)
{
    // Pull references to the pieces we need
    auto& coordInfoVecs = system.coordInfoVecs;
    auto& generalParams = system.generalParams;
    auto& domainParams  = system.domainParams;

    // Number of nodes
    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
    
    if (N == 0) {
        std::cout << "WARNING: No nodes in system, cannot relax." << std::endl;
        return 0;
    }

    // Host-side buffers for before/after positions
    std::vector<double> x_old(N), y_old(N), z_old(N);
    std::vector<double> x_new(N), y_new(N), z_new(N);

    // Force-movement accumulator
    generalParams.dx = std::numeric_limits<double>::infinity();
    int iter = 0;
    
    // Maximum iterations to prevent infinite loops
    const int max_iter = 1000000;
    
    // Progress reporting interval
    const int report_interval = 1000;

    while (generalParams.dx > generalParams.tol && iter < max_iter) {
        // 1) Snapshot old positions (device -> host)
        thrust::copy(
            coordInfoVecs.nodeLocX.begin(),
            coordInfoVecs.nodeLocX.end(),
            x_old.begin());
        thrust::copy(
            coordInfoVecs.nodeLocY.begin(),
            coordInfoVecs.nodeLocY.end(),
            y_old.begin());
        thrust::copy(
            coordInfoVecs.nodeLocZ.begin(),
            coordInfoVecs.nodeLocZ.end(),
            z_old.begin());

        // 2) Build forces, then move nodes
        system.Solve_Forces();  // member in System
        AdvancePositions(
            coordInfoVecs,
            generalParams,
            domainParams);       // free function

        // 3) Snapshot new positions 
        thrust::copy(
            coordInfoVecs.nodeLocX.begin(),
            coordInfoVecs.nodeLocX.end(),
            x_new.begin());
        thrust::copy(
            coordInfoVecs.nodeLocY.begin(),
            coordInfoVecs.nodeLocY.end(),
            y_new.begin());
        thrust::copy(
            coordInfoVecs.nodeLocZ.begin(),
            coordInfoVecs.nodeLocZ.end(),
            z_new.begin());

        // 4) Compute average displacement per node
        // This matches the paper's convergence criterion: avg_movement = sum(|displacement|) / N
        double dx_sum = 0.0;
        for (int i = 0; i < N; ++i) {
            double dx = x_new[i] - x_old[i];
            double dy = y_new[i] - y_old[i];
            double dz = z_new[i] - z_old[i];
            dx_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        
        // FIX: Divide by actual node count N, not a magic number
        // This is the average displacement per node per iteration
        generalParams.dx = dx_sum / static_cast<double>(N);

        ++iter;
        
        // Progress reporting
        if (iter % report_interval == 0) {
            std::cout << "  Relaxation iter " << iter 
                      << ": avg_displacement = " << generalParams.dx 
                      << " (tol = " << generalParams.tol << ")" << std::endl;
        }
    }
    
    if (iter >= max_iter) {
        std::cout << "WARNING: Relaxation did not converge after " << max_iter 
                  << " iterations. Final dx = " << generalParams.dx << std::endl;
    }

    return iter;
}
