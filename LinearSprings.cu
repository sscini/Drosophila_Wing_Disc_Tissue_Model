// LinearSprings.cu
// Computes linear spring forces using atomicAdd for direct force accumulation
// This approach is simpler and avoids sort_by_key synchronization issues

#include "System.h"
#include "SystemStructures.h"
#include "LinearSprings.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs)
{
    // Forces should already be zeroed in Solve_Forces()
    // If not, uncomment these lines:
    // thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    // thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    // thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

    // Create counting iterator for edge indices
    thrust::counting_iterator<int> edgeIdBegin(0);

    // Compute forces and energy using transform_reduce
    // The functor uses atomicAdd to accumulate forces directly into nodeForce arrays
    // The reduction sums up the energy from all edges
    linearSpringInfoVecs.linear_spring_energy = thrust::transform_reduce(
        thrust::make_zip_iterator(
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())) + coordInfoVecs.num_edges,
        LinearSpringFunctor(
            linearSpringInfoVecs.spring_constant,
            linearSpringInfoVecs.spring_constant_weak,
            linearSpringInfoVecs.spring_constant_vertical,
            thrust::raw_pointer_cast(linearSpringInfoVecs.edge_rest_length.data()),
            thrust::raw_pointer_cast(generalParams.edges_in_upperhem.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data()),
            generalParams.maxNodeCount),
        0.0,
        thrust::plus<double>());
    
    // Note: No sort_by_key or reduce_by_key needed!
    // The atomicAdd in the functor handles force accumulation directly.
    // This is both simpler and avoids potential synchronization issues.
}
