// LinearSprings.cu — per-edge spring constant array + atomicAdd

#include "System.h"
#include "SystemStructures.h"
#include "LinearSprings.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs)
{
    thrust::counting_iterator<int> edgeIdBegin(0);

    linearSpringInfoVecs.linear_spring_energy = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(
            edgeIdBegin,
            coordInfoVecs.edges2Nodes_1.begin(),
            coordInfoVecs.edges2Nodes_2.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(
            edgeIdBegin,
            coordInfoVecs.edges2Nodes_1.begin(),
            coordInfoVecs.edges2Nodes_2.begin())) + coordInfoVecs.num_edges,
        LinearSpringFunctor(
            thrust::raw_pointer_cast(linearSpringInfoVecs.edge_spring_k.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.edge_rest_length.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data()),
            generalParams.maxNodeCount),
        0.0,
        thrust::plus<double>());
}
