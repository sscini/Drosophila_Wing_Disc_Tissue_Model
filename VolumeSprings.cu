// VolumeSprings.cu
// Applies volume-conserving forces to maintain enclosed volume

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
        std::cout << "ComputeVolumeSprings: No prisms, skipping." << std::endl;
        return;
    }

    const double kv = generalParams.volume_spring_constant;
    if (kv == 0.0) {
        std::cout << "ComputeVolumeSprings: kv=0, skipping." << std::endl;
        return;
    }

    const double Omega_s = generalParams.current_total_volume;
    const double Omega0 = generalParams.eq_total_volume;
    
    // Avoid division issues if volume is near zero (mesh inversion)
    if (std::fabs(Omega_s) < 1e-12) {
        std::cout << "ComputeVolumeSprings: Volume near zero, skipping." << std::endl;
        return;
    }

    // Force prefactor: F = -?E/?r = -2*k_v*(O - O0)*?O/?r
    const double prefactor = -2.0 * kv * (Omega_s - Omega0);
//
//    // Debug output
//    std::cout << "ComputeVolumeSprings: Omega=" << Omega_s 
//              << ", Omega0=" << Omega0
//              << ", dOmega=" << (Omega_s - Omega0)
//              << ", kv=" << kv
//              << ", prefactor=" << prefactor << std::endl;

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
//    
//    // Debug: Print force on first few nodes
//    std::cout << "  Sample volume forces:" << std::endl;
//    for (int i = 0; i < std::min(3, Nnodes); ++i) {
//        std::cout << "    Node " << i << ": F=(" 
//                  << coordInfoVecs.nodeForceX[i] << ", "
//                  << coordInfoVecs.nodeForceY[i] << ", "
//                  << coordInfoVecs.nodeForceZ[i] << ")" << std::endl;
//    }
}
