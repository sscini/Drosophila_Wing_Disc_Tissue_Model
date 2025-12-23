#include "System.h"
#include "SystemStructures.h"
#include "VolumeSprings.h"

#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <cmath>
#include <iostream>

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs) {

    // ---- infer layout as in ComputeVolume ----
    const int N = generalParams.maxNodeCount;               // nodes per layer
    const size_t total_nodes = coordInfoVecs.nodeLocX.size();

    if (N <= 0 || total_nodes % N != 0) {
        std::cout << "VolumeSprings: bad N or total_nodes, skipping.\n";
        return;
    }

    const int L = generalParams.num_layers;        // # layers
    if (L < 2) {
        std::cout << "VolumeSprings: L < 2, skipping.\n";
        return;
    }

    const int Ttot = coordInfoVecs.num_triangles;
////    if (Ttot <= 0 || Ttot != 0) {
////        std::cout << "VolumeSprings: bad num_triangles = "
////                  << Ttot << ", L = " << L << ", skipping.\n";
////        return;
////    }
//
    const int Tper      = Ttot / L;                         // triangles per layer
    const int numPrisms = (L - 1) * Tper;                   // prisms between layers

    // ---- global volume values (must be set by ComputeVolume) ----
    // Use signed volume for orientation info:
    const double Omega_s   = generalParams.current_total_volume;
    // Use absolute for physical magnitude:
    const double Omega_abs = std::fabs(Omega_s);
    const double Omega0    = generalParams.eq_total_volume;          // target
    const double kv        = generalParams.volume_spring_constant;   // stiffness

    double signOmega = 0.0;
    if (Omega_s > 0.0)      signOmega = 1.0;
    else if (Omega_s < 0.0) signOmega = -1.0;

    // Degenerate or uninitialized volume: skip forces
    if (Omega_abs < 1e-12 || Omega0 <= 0.0 || kv == 0.0 || signOmega == 0.0) {
        // std::cout << "VolumeSprings: degenerate volume or kv=0, skipping.\n";
        return;
    }

    // Force prefactor from E = k_v (Omega_abs - Omega0)^2:
    // F_n = -2 k_v (Omega_abs - Omega0) sign(Omega_s) ?Omega/?r_n
    const double prefactor = -2.0 * kv * (Omega_abs - Omega0) * signOmega;

    // ---- zip iterator over (node_id, dummy_int, Fx, Fy, Fz) ----
    // We use id_value_expanded twice; we ignore the second entry in the functor.
    auto begin = thrust::make_zip_iterator(
        thrust::make_tuple(
            auxVecs.id_value_expanded.begin(),   // node_id
            auxVecs.id_value_expanded.begin(),   // dummy / bucket
            coordInfoVecs.nodeForceX.begin(),    // forces to be updated
            coordInfoVecs.nodeForceY.begin(),
            coordInfoVecs.nodeForceZ.begin()));

    auto end = begin + auxVecs.id_value_expanded.size();

    // ---- construct functor ----
    VolumeSpringPrismFunctor functor(
        prefactor,
        N,
        Tper,
        L,
        numPrisms,
        thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()),
        thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()),
        thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()));

    // ---- apply volume forces ----
    thrust::transform(begin, end, begin, functor);
}
