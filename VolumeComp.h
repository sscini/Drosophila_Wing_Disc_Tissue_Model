#ifndef VOLUMECOMP_H_
#define VOLUMECOMP_H_

#include "SystemStructures.h"
#include <cmath>

// Forward declarations
struct GeneralParams;
struct CoordInfoVecs;
struct LinearSpringInfoVecs;
struct LJInfoVecs;
struct PrismInfoVecs;

/**
 * ComputeVolume - Calculate the total enclosed volume of the mesh
 * 
 * This function computes the global volume by summing signed tetrahedron
 * volumes from all prisms. Each prism is decomposed into 3 tetrahedra
 * sharing a common apex (P1).
 * 
 * Tetrahedron decomposition for prism (P1,P2,P3,P4,P5,P6):
 *   Tet 1: (P2, P3, P4, P1) - bottom triangle edge to top vertex
 *   Tet 2: (P2, P4, P6, P1) - diagonal connection
 *   Tet 3: (P4, P5, P6, P1) - top triangle
 * 
 * The signed volume uses the scalar triple product:
 *   6V = (r_i - r_l) · [(r_j - r_l) × (r_k - r_l)]
 * 
 * Stores results in:
 *   generalParams.current_total_volume - signed volume (can be negative if inverted)
 *   generalParams.true_current_total_volume - diagnostic only
 * 
 * @param generalParams  General simulation parameters
 * @param coordInfoVecs  Node coordinate vectors
 * @param linearSpringInfoVecs  Spring information
 * @param ljInfoVecs  Lennard-Jones parameters
 * @param prismInfoVecs  Prism connectivity (P1-P6 for each prism)
 */
void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs,
    PrismInfoVecs& prismInfoVecs);

#endif // VOLUMECOMP_H_