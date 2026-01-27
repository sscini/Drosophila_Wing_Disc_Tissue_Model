// VolumeSprings.cu - FIXED VERSION
// Applies volume-conserving forces to maintain enclosed volume
// With safeguards against numerical instability

#include "System.h"
#include "SystemStructures.h"
#include "VolumeSprings.h"

#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <cmath>
#include <iostream>

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs)
{
    const int P = prismInfoVecs.num_prisms; 
    if (P <= 0) {
        return;
    }

    const double kv = generalParams.volume_spring_constant;
    if (kv == 0.0) {
        return;
    }

    const double Omega_s = generalParams.current_total_volume;
    const double Omega0 = generalParams.eq_total_volume;
    
    // ==================== SAFEGUARD 1: Check for NaN/Inf ====================
    if (std::isnan(Omega_s) || std::isinf(Omega_s)) {
        std::cout << "WARNING: Volume is NaN/Inf (" << Omega_s << "), skipping volume springs." << std::endl;
        return;
    }
    
    // ==================== SAFEGUARD 2: Check for negative volume ====================
    // Negative volume indicates mesh inversion - catastrophic state
    if (Omega_s < 0.0) {
        std::cout << "WARNING: Negative volume detected (" << Omega_s 
                  << "). Mesh may be inverted. Skipping volume springs to prevent further instability." << std::endl;
        return;
    }
    
    // ==================== SAFEGUARD 3: Check for extreme volume change ====================
    double volume_ratio = Omega_s / Omega0;
    double volume_diff = Omega_s - Omega0;
    double prefactor;
    
    // If volume changed by more than 50%, clamp the force to prevent explosion
    if (volume_ratio < 0.5 || volume_ratio > 2.0) {
        std::cout << "WARNING: Extreme volume change (ratio=" << volume_ratio 
                  << "). Clamping volume spring force." << std::endl;
        
        // Clamp the effective volume difference to at most 50% of equilibrium volume
        double max_diff = 0.5 * std::fabs(Omega0);
        double clamped_diff = std::max(-max_diff, std::min(max_diff, volume_diff));
        prefactor = -2.0 * kv * clamped_diff;
    }
    else {
        // Normal calculation
        prefactor = -2.0 * kv * volume_diff;
    }
    
    // ==================== SAFEGUARD 4: Check for NaN/Inf prefactor ====================
    if (std::isnan(prefactor) || std::isinf(prefactor)) {
        std::cout << "WARNING: Volume spring prefactor is NaN/Inf. Skipping." << std::endl;
        return;
    }
    
    // ==================== SAFEGUARD 5: Clamp maximum force magnitude ====================
    // This prevents runaway forces from destroying the simulation
    const double max_prefactor = 1e6;  // Adjust based on your simulation scale
    if (std::fabs(prefactor) > max_prefactor) {
        std::cout << "WARNING: Volume spring prefactor too large (" << prefactor 
                  << "). Clamping to " << max_prefactor << std::endl;
        prefactor = (prefactor > 0) ? max_prefactor : -max_prefactor;
    }
    
    generalParams.volume_energy = kv * volume_diff * volume_diff ;
    const int Nnodes = (int)coordInfoVecs.nodeLocX.size();

    // Create counting iterator for node IDs
    auto ids0 = thrust::make_counting_iterator<int>(0);

    // Zip iterator over (node_id, bucket_id, Fx, Fy, Fz)
    auto begin = thrust::make_zip_iterator(
        thrust::make_tuple(
            ids0,
            ids0,  // dummy bucket ID
            coordInfoVecs.nodeForceX.begin(),
            coordInfoVecs.nodeForceY.begin(),
            coordInfoVecs.nodeForceZ.begin()));

    auto end = begin + Nnodes;

    // Create functor
    VolumeSpringPrismFunctor functor(
        prefactor,
        P,
        Nnodes,
        thrust::raw_pointer_cast(prismInfoVecs.P1.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P2.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P3.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P4.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P5.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P6.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()));

    // Apply transform to compute and add volume forces
    thrust::transform(begin, end, begin, functor);
}
