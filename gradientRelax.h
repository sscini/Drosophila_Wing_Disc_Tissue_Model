// ============================================================================
// gradientRelax.h - Force-Based Convergence Version
// ============================================================================

#ifndef GRADIENTRELAX_H_
#define GRADIENTRELAX_H_

#include "System.h"
#include "NodeAdvance.h"

/// Run overdamped (gradient-descent) relaxation on the system
/// until max force magnitude falls below threshold.
/// Returns number of iterations taken.
///
/// Default force tolerance is 1.0 - adjust in the .cu file if needed.
int relaxUntilConverged(System& system);

/// Version with configurable parameters
/// @param force_tolerance      Stop when max|F| < this value
/// @param displacement_tolerance Secondary criterion (usually can be loose)
/// @param max_iterations       Safety limit to prevent infinite loops
/// @param print_every          Print progress every N iterations (0 = silent)
int relaxUntilConvergedWithParams(
    System& system,
    double force_tolerance,
    double displacement_tolerance,
    int max_iterations,
    int print_every);

#endif // GRADIENTRELAX_H_



//#ifndef GRADIENTRELAX_H_
//#define GRADIENTRELAX_H_
//
//#include "System.h"        // for class System and Solve_Forces()
//#include "NodeAdvance.h"   // for AdvancePositions(...)
//  
///// Run overdamped (gradient-descent) relaxation on the system
///// until forces fall below generalParams.tol. Returns number of iterations.
//int relaxUntilConverged(System& system);
//
//#endif // GRADIENTRELAX_H_

// remove the force relaxation. make it based on node movement