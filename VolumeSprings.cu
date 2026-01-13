// VolumeSprings.cu
#include "System.h"
#include "SystemStructures.h"
#include "VolumeSprings.h"

#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <cmath>

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    //AuxVecs& auxVecs,
    PrismInfoVecs& prismInfoVecs)
{
    const int P = prismInfoVecs.num_prisms; 
    if (P <= 0) return;

    const double Omega_s = generalParams.current_total_volume; // signed global volume
    const double Omega0  = generalParams.eq_total_volume;
    const double kv      = generalParams.volume_spring_constant;
    
    if (kv == 0.0) return;
    
    // (optional) avoid insane forces if volume is near zero due to inversion
    if (std::fabs(Omega_s) < 1e-12) return;
    if (kv ==0.0) return; 
    
    const double prefactor = -2.0 * kv * (Omega_s - Omega0);

    //const double prefactor = -2.0 * kv * (Omega_abs - Omega0) * signOmega;

    const int Nnodes = (int)coordInfoVecs.nodeLocX.size(); // true node count

    // node ids: 0..Nnodes-1
    auto ids0 = thrust::make_counting_iterator<int>(0);

    // Zip: (node_id, dummy_bucket, Fx, Fy, Fz)
    auto begin = thrust::make_zip_iterator(
        thrust::make_tuple(
            ids0,
            ids0, // dummy; keep your tuple shape unchanged
            coordInfoVecs.nodeForceX.begin(),
            coordInfoVecs.nodeForceY.begin(),
            coordInfoVecs.nodeForceZ.begin()));

    auto end = begin + Nnodes;  // IMPORTANT: match force vector length

    VolumeSpringPrismFunctor functor(
        prefactor,
        P,
        Nnodes, // pass for bounds checks
        thrust::raw_pointer_cast(prismInfoVecs.P1.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P2.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P3.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P4.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P5.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P6.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()));

    thrust::transform(begin, end, begin, functor);
}
