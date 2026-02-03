#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

/**
 * ComputeVolumeSprings - Apply volume-conserving forces to all nodes
 * 
 * Implements the force: F_n = -dE_vol/dr_n = -2*k_v*(V - V0)*dV/dr_n
 * 
 * The gradient dV/dr_n is computed by summing contributions from all
 * prisms that contain node n, decomposed into tetrahedra.
 *
 * Forces are accumulated via atomicAdd (matching LinearSprings pattern)
 * to avoid overwriting forces from other springs.
 */
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs);

// Note: The VolumeSpringAtomicFunctor is defined internally in VolumeSprings.cu.
// The old VolumeSpringPrismFunctor (thrust::transform-based) has been removed 
// because it caused a force overwrite bug when combined with LinearSprings' 
// atomicAdd pattern.

#endif // VOLUMESPRINGS_H_
