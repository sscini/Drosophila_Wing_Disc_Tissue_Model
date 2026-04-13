// ============================================================================
// gradientRelax.cu - Force-Based Convergence
//
// FIXES from previous version:
//   1. avgDisplacement variable shadowing bug — was staying at infinity
//   2. Device vector allocation moved OUTSIDE the loop (was reallocating
//      GPU memory every iteration — massive performance hit)
//   3. Removed redundant std::vector declarations
// ============================================================================

#include "gradientRelax.h"
#include <vector>
#include <thrust/copy.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <limits>
#include <cmath>
#include <iostream>

// Functor to compute displacement between old and new positions
struct DisplacementFunctor {
    __host__ __device__
    double operator()(const thrust::tuple<double,double,double,double,double,double>& t) const {
        double dx = thrust::get<0>(t) - thrust::get<3>(t);
        double dy = thrust::get<1>(t) - thrust::get<4>(t);
        double dz = thrust::get<2>(t) - thrust::get<5>(t);
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
};

// Functor to compute force magnitude from (fx, fy, fz)
struct ForceMagnitudeFunctor {
    __host__ __device__
    double operator()(const thrust::tuple<double, double, double>& f) const {
        double fx = thrust::get<0>(f);
        double fy = thrust::get<1>(f);
        double fz = thrust::get<2>(f);
        return sqrt(fx*fx + fy*fy + fz*fz);
    }
};

int relaxUntilConvergedWithParams(
    System& system,
    double force_tolerance,
    double displacement_tolerance,
    int max_iterations,
    int print_every)
{
    auto& coordInfoVecs = system.coordInfoVecs;
    auto& generalParams = system.generalParams;
    auto& domainParams  = system.domainParams;

    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());

    // Allocate device vectors ONCE outside the loop
    thrust::device_vector<double> x_old(N);
    thrust::device_vector<double> y_old(N);
    thrust::device_vector<double> z_old(N);

    int iter = 0;
    double maxForce = std::numeric_limits<double>::infinity();
    double avgDisplacement = std::numeric_limits<double>::infinity();

    while (iter < max_iterations) {

        // Snapshot current positions (device-to-device copy, no realloc)
        x_old = coordInfoVecs.nodeLocX;
        y_old = coordInfoVecs.nodeLocY;
        z_old = coordInfoVecs.nodeLocZ;

        // Compute forces
        system.Solve_Forces();

        // Compute max force magnitude on GPU
        maxForce = thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(
                coordInfoVecs.nodeForceX.end(),
                coordInfoVecs.nodeForceY.end(),
                coordInfoVecs.nodeForceZ.end())),
            ForceMagnitudeFunctor(),
            0.0,
            thrust::maximum<double>());

        // Advance positions
        AdvancePositions(coordInfoVecs, generalParams, domainParams);

        // Compute average displacement on GPU
        double totalDisp = thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(
                coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocZ.begin(),
                x_old.begin(), y_old.begin(), z_old.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(
                coordInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocZ.end(),
                x_old.end(), y_old.end(), z_old.end())),
            DisplacementFunctor(),
            0.0,
            thrust::plus<double>());
        avgDisplacement = totalDisp / N;  // NO 'double' keyword — assigns to outer variable

        // Check convergence
        if (maxForce < force_tolerance) {
            generalParams.dx = avgDisplacement;
            if (print_every > 0) {
                std::cout << "  Converged: iter=" << iter
                          << ", maxF=" << maxForce << std::endl;
            }
            return iter;
        }

        // Progress output
        if (print_every > 0 && iter % print_every == 0) {
            std::cout << "  iter=" << iter << ", maxF=" << maxForce
                      << ", avgDisp=" << avgDisplacement << std::endl;
        }

        ++iter;
    }

    // Hit max iterations
    std::cout << "  WARNING: Max iterations reached, maxF=" << maxForce << std::endl;
    generalParams.dx = avgDisplacement;
    return iter;
}



// uncomment whole block below for previous working version of gradient relax


//// ============================================================================
//// gradientRelax.cu - Force-Based Convergence Version
////
//// This version uses maximum force magnitude as the convergence criterion
//// instead of displacement. This is more physically meaningful because:
////   - Equilibrium means forces are balanced (F ˜ 0), not that motion stopped
////   - Displacement depends on timestep; force doesn't
////   - Avoids "false convergence" when dt is small
////
//// Replace your existing gradientRelax.cu with this file.
//// ============================================================================
//
//#include "gradientRelax.h"
//#include <vector>
//#include <thrust/copy.h>
//#include <thrust/transform_reduce.h>
//#include <thrust/functional.h>
//#include <limits>
//#include <cmath>
//#include <iostream>
//
//struct DisplacementFunctor {
//    __host__ __device__
//    double operator()(const thrust::tuple<double,double,double,double,double,double>& t) const {
//        double dx = thrust::get<0>(t) - thrust::get<3>(t);
//        double dy = thrust::get<1>(t) - thrust::get<4>(t);
//        double dz = thrust::get<2>(t) - thrust::get<5>(t);
//        return sqrt(dx*dx + dy*dy + dz*dz);
//    }
//};
//
//// Functor to compute force magnitude from (fx, fy, fz)
//struct ForceMagnitudeFunctor {
//    __host__ __device__
//    double operator()(const thrust::tuple<double, double, double>& f) const {
//        double fx = thrust::get<0>(f);
//        double fy = thrust::get<1>(f);
//        double fz = thrust::get<2>(f);
//        return sqrt(fx*fx + fy*fy + fz*fz);
//    }
//};
//
////int relaxUntilConverged(System& system)
////{
////    // Pull references to the pieces we need
////    auto& coordInfoVecs = system.coordInfoVecs;
////    auto& generalParams = system.generalParams;
////    auto& domainParams  = system.domainParams;
////
////    // Number of nodes
////    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
////
////    // ========================================================================
////    // CONVERGENCE PARAMETERS
////    // ========================================================================
////    
////    // Force-based convergence threshold
////    // Converge when max|F| < force_tol
////    double force_tol = 1.0;  // Adjust based on your spring constants
////    
////    // Also keep displacement check as a secondary criterion
////    double disp_tol = generalParams.tol;  // Use the existing tolerance
////    
////    // Maximum iterations to prevent infinite loops
////    int max_iterations = 100000;
////    
////    // How often to print progress (0 = never)
////    int print_interval = 0;  // Set to 1000 or so for debugging
////    
////    // ========================================================================
////    // RELAXATION LOOP
////    // ========================================================================
////    
////    // Host-side buffers for position tracking (optional, for displacement check)
////    std::vector<double> x_old(N), y_old(N), z_old(N);
////    
////    int iter = 0;
////    double maxForce = std::numeric_limits<double>::infinity();
////    double avgDisplacement = std::numeric_limits<double>::infinity();
////    
////    while (iter < max_iterations) {
////        
////        // 1) Snapshot old positions (for displacement calculation)
////        thrust::copy(coordInfoVecs.nodeLocX.begin(),
////                     coordInfoVecs.nodeLocX.end(),
////                     x_old.begin());
////        thrust::copy(coordInfoVecs.nodeLocY.begin(),
////                     coordInfoVecs.nodeLocY.end(),
////                     y_old.begin());
////        thrust::copy(coordInfoVecs.nodeLocZ.begin(),
////                     coordInfoVecs.nodeLocZ.end(),
////                     z_old.begin());
////
////        // 2) Compute forces
////        system.Solve_Forces();
////        
////        // 3) Compute maximum force magnitude (PRIMARY convergence criterion)
////        maxForce = thrust::transform_reduce(
////            thrust::make_zip_iterator(thrust::make_tuple(
////                coordInfoVecs.nodeForceX.begin(),
////                coordInfoVecs.nodeForceY.begin(),
////                coordInfoVecs.nodeForceZ.begin())),
////            thrust::make_zip_iterator(thrust::make_tuple(
////                coordInfoVecs.nodeForceX.end(),
////                coordInfoVecs.nodeForceY.end(),
////                coordInfoVecs.nodeForceZ.end())),
////            ForceMagnitudeFunctor(),
////            0.0,
////            thrust::maximum<double>());
////        
////        // 4) Check force-based convergence BEFORE advancing positions
////        if (maxForce < force_tol) {
////            // Store final values for external access
////            generalParams.dx = avgDisplacement;
////            
////            if (print_interval > 0) {
////                std::cout << "  [relaxUntilConverged] Converged! iter=" << iter 
////                          << ", maxF=" << maxForce << std::endl;
////            }
////            return iter;
////        }
////        
////        // 5) Advance positions
////        AdvancePositions(coordInfoVecs, generalParams, domainParams);
////
////        // 6) Compute displacement (SECONDARY criterion, for monitoring)
////        double dx_sum = 0.0;
////        std::vector<double> x_new(N), y_new(N), z_new(N);
////        thrust::copy(coordInfoVecs.nodeLocX.begin(),
////                     coordInfoVecs.nodeLocX.end(),
////                     x_new.begin());
////        thrust::copy(coordInfoVecs.nodeLocY.begin(),
////                     coordInfoVecs.nodeLocY.end(),
////                     y_new.begin());
////        thrust::copy(coordInfoVecs.nodeLocZ.begin(),
////                     coordInfoVecs.nodeLocZ.end(),
////                     z_new.begin());
////        
////        for (int i = 0; i < N; ++i) {
////            double dx = x_new[i] - x_old[i];
////            double dy = y_new[i] - y_old[i];
////            double dz = z_new[i] - z_old[i];
////            dx_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
////        }
////        avgDisplacement = dx_sum / N;
////        
////        // 7) Optional: Also check displacement-based convergence
////        //    (useful if forces are tiny but non-zero due to numerical noise)
////        if (avgDisplacement < disp_tol && maxForce < force_tol * 10) {
////            generalParams.dx = avgDisplacement;
////            
////            if (print_interval > 0) {
////                std::cout << "  [relaxUntilConverged] Converged (disp)! iter=" << iter 
////                          << ", maxF=" << maxForce << ", avgDisp=" << avgDisplacement << std::endl;
////            }
////            return iter;
////        }
////        
////        // 8) Progress output
////        if (print_interval > 0 && iter % print_interval == 0) {
////            std::cout << "  [relaxUntilConverged] iter=" << iter 
////                      << ", maxF=" << maxForce 
////                      << ", avgDisp=" << avgDisplacement << std::endl;
////        }
////
////        ++iter;
////    }
////    
////    // If we hit max iterations, warn and return
////    std::cout << "  [relaxUntilConverged] WARNING: Hit max iterations (" << max_iterations 
////              << "), maxF=" << maxForce << ", avgDisp=" << avgDisplacement << std::endl;
////    
////    generalParams.dx = avgDisplacement;
////    return iter;
////}
////
//// ============================================================================
//// Alternative version with configurable parameters
//// ============================================================================
//
//int relaxUntilConvergedWithParams(
//    System& system,
//    double force_tolerance,      // Stop when max|F| < this
//    double displacement_tolerance, // Secondary criterion
//    int max_iterations,          // Safety limit
//    int print_every)             // Print progress every N iterations (0 = silent)
//{
//    auto& coordInfoVecs = system.coordInfoVecs;
//    auto& generalParams = system.generalParams;
//    auto& domainParams  = system.domainParams;
//
//    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
//    
//    std::vector<double> x_old(N), y_old(N), z_old(N);
//    
//    int iter = 0;
//    double maxForce = std::numeric_limits<double>::infinity();
//    double avgDisplacement = std::numeric_limits<double>::infinity();
//    
//    while (iter < max_iterations) {
//        
//                thrust::device_vector<double> x_old = coordInfoVecs.nodeLocX;
//        thrust::device_vector<double> y_old = coordInfoVecs.nodeLocY;
//        thrust::device_vector<double> z_old = coordInfoVecs.nodeLocZ;
//
//        system.Solve_Forces();
//        
//        // SINGLE max force computation using GPU reduction
//        maxForce = thrust::transform_reduce(
//            thrust::make_zip_iterator(thrust::make_tuple(
//                coordInfoVecs.nodeForceX.begin(),
//                coordInfoVecs.nodeForceY.begin(),
//                coordInfoVecs.nodeForceZ.begin())),
//            thrust::make_zip_iterator(thrust::make_tuple(
//                coordInfoVecs.nodeForceX.end(),
//                coordInfoVecs.nodeForceY.end(),
//                coordInfoVecs.nodeForceZ.end())),
//            ForceMagnitudeFunctor(),
//            0.0,
//            thrust::maximum<double>());
//        
//        
//        // Check convergence BEFORE position advance
////        if (maxForce < force_tolerance) {
////            if (print_every > 0) {
////                std::cout << "  Converged: iter=" << iter << ", maxF=" << maxForce << std::endl;
////            }
////            return iter;
////        }
//        
//        //std::cout << "DEBUG: dt = " << generalParams.dt << ", maxF = " << maxForce << std::endl;
//   
//        
//        //Advance positions
//        AdvancePositions(coordInfoVecs, generalParams, domainParams);
//        
////        std::cout << "Max|F|=" << maxF  << std::endl;
////        
////        // Compute maximum force magnitude
////        maxForce = thrust::transform_reduce(
////            thrust::make_zip_iterator(thrust::make_tuple(
////                coordInfoVecs.nodeForceX.begin(),
////                coordInfoVecs.nodeForceY.begin(),
////                coordInfoVecs.nodeForceZ.begin())),
////            thrust::make_zip_iterator(thrust::make_tuple(
////                coordInfoVecs.nodeForceX.end(),
////                coordInfoVecs.nodeForceY.end(),
////                coordInfoVecs.nodeForceZ.end())),
////            ForceMagnitudeFunctor(),
////            0.0,
////            thrust::maximum<double>());
////            
////        static double cumulative_max_displacement = 0.0;
////        static double initial_x0 = -999999;
////        
////        if (initial_x0 < -999998) {
////            initial_x0 = coordInfoVecs.nodeLocX[0];
////        }
////        
//        
//        
////         //Advance positions
////        AdvancePositions(coordInfoVecs, generalParams, domainParams);
////        
//        
//        // Track cumulative displacement of node 0
////        double current_x0 = coordInfoVecs.nodeLocX[0];
////        double node0_total_disp = fabs(current_x0 - initial_x0);
////        
////        if (iter % 100 == 0) {
////            std::cout << "Node 0: initial_x=" << initial_x0 
////                      << " current_x=" << current_x0 
////                      << " total_displacement=" << node0_total_disp << std::endl;
////        }
//        // Check convergence
//        if (maxForce < force_tolerance) {
//            generalParams.dx = avgDisplacement;
//            if (print_every > 0) {
//                std::cout << "  Converged: iter=" << iter << ", maxF=" << maxForce << std::endl;
//                double totalDisp = thrust::transform_reduce(
//                    thrust::make_zip_iterator(thrust::make_tuple(
//                        coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocZ.begin(),
//                        x_old.begin(), y_old.begin(), z_old.begin())),
//                    thrust::make_zip_iterator(thrust::make_tuple(
//                        coordInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocZ.end(),
//                        x_old.end(), y_old.end(), z_old.end())),
//                    DisplacementFunctor(),
//                    0.0,
//                    thrust::plus<double>());
//                double avgDisplacement = totalDisp / N;
//                
//                generalParams.dx = avgDisplacement;
//            }
//            return iter;
//        
//        }
//        
//
////        // Compute displacement
////        double dx_sum = 0.0;
////        for (int i = 0; i < N; ++i) {
////            double dx = coordInfoVecs.nodeLocX[i] - x_old[i];
////            double dy = coordInfoVecs.nodeLocY[i] - y_old[i];
////            double dz = coordInfoVecs.nodeLocZ[i] - z_old[i];
////            dx_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
////        }
////        avgDisplacement = dx_sum / N;
//
//        double totalDisp = thrust::transform_reduce(
//            thrust::make_zip_iterator(thrust::make_tuple(
//                coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocZ.begin(),
//                x_old.begin(), y_old.begin(), z_old.begin())),
//            thrust::make_zip_iterator(thrust::make_tuple(
//                coordInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocZ.end(),
//                x_old.end(), y_old.end(), z_old.end())),
//            DisplacementFunctor(),
//            0.0,
//            thrust::plus<double>());
//        double avgDisplacement = totalDisp / N;
//        
////       //  Progress output
////        if (print_every > 0 && iter % print_every == 0) {
////            std::cout << "  iter=" << iter << ", maxF=" << maxForce 
////                      << ", avgDisp=" << avgDisplacement << std::endl;
////        }
//
//        ++iter;
//    }
//    
//    std::cout << "  WARNING: Max iterations reached, maxF=" << maxForce << std::endl;
//    generalParams.dx = avgDisplacement;
//    return iter;
//}
//
//
//
//
//
//
//
////
////#include "gradientRelax.h"
////#include <vector>            // for std::vector
////#include <thrust/copy.h>     // for thrust::copy
////#include <limits>            // for std::numeric_limits
////#include <cmath>             // for std::sqrt
////
////int relaxUntilConverged(System& system)
////{
////    // Pull references to the pieces we need
////    auto& coordInfoVecs = system.coordInfoVecs;
////    auto& generalParams = system.generalParams;
////    auto& domainParams  = system.domainParams;
////
////    // Number of nodes
////    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
////
////    // Host-side buffers for before/after positions
////    std::vector<double> x_old(N), y_old(N), z_old(N);
////    std::vector<double> x_new(N), y_new(N), z_new(N);
////
////    // Force-movement accumulator
////    generalParams.dx = std::numeric_limits<double>::infinity();
////    int iter = 0;
////
////    while (generalParams.dx > generalParams.tol) {
////        // 1) Snapshot old positions (device ? host)
////        thrust::copy(
////            coordInfoVecs.nodeLocX.begin(),
////            coordInfoVecs.nodeLocX.end(),
////            x_old.begin());
////        thrust::copy(
////            coordInfoVecs.nodeLocY.begin(),
////            coordInfoVecs.nodeLocY.end(),
////            y_old.begin());
////        thrust::copy(
////            coordInfoVecs.nodeLocZ.begin(),
////            coordInfoVecs.nodeLocZ.end(),
////            z_old.begin());
////
////        // 2) Build forces, then move nodes
////        system.Solve_Forces();  // member in System
////        // After Solve_Forces() in relaxUntilConverged:
////        double max_F = 0.0;
////        thrust::host_vector<double> hfx = coordInfoVecs.nodeForceX;
////        thrust::host_vector<double> hfy = coordInfoVecs.nodeForceY;
////        thrust::host_vector<double> hfz = coordInfoVecs.nodeForceZ;
////        for (int i = 0; i < N; i++) {
////            double F = sqrt(hfx[i]*hfx[i] + hfy[i]*hfy[i] + hfz[i]*hfz[i]);
////            max_F = std::max(max_F, F);
////        }
////        std::cout << "Max |F| = " << max_F << std::endl;
////        AdvancePositions(
////            coordInfoVecs,
////            generalParams,
////            domainParams);       // free function
////
////        // 3) Snapshot new positions 
////        thrust::copy(
////            coordInfoVecs.nodeLocX.begin(),
////            coordInfoVecs.nodeLocX.end(),
////            x_new.begin());
////        thrust::copy(
////            coordInfoVecs.nodeLocY.begin(),
////            coordInfoVecs.nodeLocY.end(),
////            y_new.begin());
////        thrust::copy(
////            coordInfoVecs.nodeLocZ.begin(),
////            coordInfoVecs.nodeLocZ.end(),
////            z_new.begin());
////
////        // 4) Compute total L2-movement across all nodes
////        double dx_sum = 0.0;
////        for (int i = 0; i < N; ++i) {
////            double dx = x_new[i] - x_old[i];
////            double dy = y_new[i] - y_old[i];
////            double dz = z_new[i] - z_old[i];
////            dx_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
////        }
////        generalParams.dx = dx_sum/N;
////
////        ++iter;
////    }
////
////    return iter;
////}
