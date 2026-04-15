#include "System.h"
#include "SystemStructures.h"
#include "AreaTriangles.h"

// Function for computing the area-based forces on nodes.
void ComputeAreaTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

        // Create a counting iterator for element IDs.
        thrust::counting_iterator<int> elemId(0);
        
        // Initialize temporary vectors with zero values.
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
            areaTriangleInfoVecs.tempNodeForceXReduced.end(),
            0.0
            );
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
            areaTriangleInfoVecs.tempNodeForceYReduced.end(),
            0.0
            );
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceZReduced.begin(),
            areaTriangleInfoVecs.tempNodeForceZReduced.end(),
            0.0
            );
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
            areaTriangleInfoVecs.tempNodeForceXUnreduced.end(),
            0.0
            );
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
            areaTriangleInfoVecs.tempNodeForceYUnreduced.end(),
            0.0
            );
        thrust::fill(
            areaTriangleInfoVecs.tempNodeForceZUnreduced.begin(),
            areaTriangleInfoVecs.tempNodeForceZUnreduced.end(),
            0.0
            );
    
        
       // Compute the area-based triangle energy using transform_reduce.
       areaTriangleInfoVecs.area_triangle_energy = thrust::transform_reduce( 
			      thrust::make_zip_iterator(
				        thrust::make_tuple(
                    elemId,
          					coordInfoVecs.triangles2Nodes_1.begin(),
          					coordInfoVecs.triangles2Nodes_2.begin(),
                    coordInfoVecs.triangles2Nodes_3.begin(),
                    coordInfoVecs.triangles2Edges_1.begin(),
                    coordInfoVecs.triangles2Edges_2.begin(),
                    coordInfoVecs.triangles2Edges_3.begin()
                    )
               ),
      			thrust::make_zip_iterator(
      				    thrust::make_tuple(
                          elemId,
                					coordInfoVecs.triangles2Nodes_1.begin(),
                					coordInfoVecs.triangles2Nodes_2.begin(),
                          coordInfoVecs.triangles2Nodes_3.begin(),
                          coordInfoVecs.triangles2Edges_1.begin(),
                          coordInfoVecs.triangles2Edges_2.begin(),
                          coordInfoVecs.triangles2Edges_3.begin()
                          )
            ) + coordInfoVecs.num_triangles,
            AreaSpringFunctor( 
                    generalParams.SCALE_TYPE,
                    generalParams.nonuniform_wall_weakening_area,
                    generalParams.maxSpringScaler_area,
                    generalParams.scaling_pow,
                    generalParams.gausssigma,
                    generalParams.hilleqnconst,
                    generalParams.hilleqnpow,
                    thrust::raw_pointer_cast(coordInfoVecs.scaling_per_edge.data()),
                    areaTriangleInfoVecs.initial_area,
                    areaTriangleInfoVecs.spring_constant,
                    areaTriangleInfoVecs.spring_constant_weak,
                    thrust::raw_pointer_cast(generalParams.triangles_in_upperhem.data()),
                    thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                    thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                    thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
                    // NEW: pass node layer flags for cross-layer triangle filtering
                    thrust::raw_pointer_cast(generalParams.nodes_in_upperhem.data()),
                    thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data()),
                    thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceXUnreduced.data()),
                    thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceYUnreduced.data()),
                    thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceZUnreduced.data())
                    ),
           0.0, 
           thrust::plus<double>() 
           );
        
        
           // We now have unreducede forces. Sort these forces by node ID.
           // First key, then value. Vectors are returned sorted. 	
           thrust::sort_by_key(
                   areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
                   areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor*coordInfoVecs.num_triangles),
                   thrust::make_zip_iterator(
                       thrust::make_tuple(
      			                areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
      					            areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
      					            areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<int>());

            // Reduce the forces by node ID.
            int endKey = thrust::get<0>(
                  thrust::reduce_by_key(
                  areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
                  areaTriangleInfoVecs.tempNodeIdUnreduced.begin() + (areaTriangleInfoVecs.factor*coordInfoVecs.num_triangles),
                  thrust::make_zip_iterator(
                      thrust::make_tuple(
                          areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                          areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                          areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
                  areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                  thrust::make_zip_iterator(
                      thrust::make_tuple(
                          areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                          areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                          areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
                  thrust::equal_to<int>(), CVec3Add())) - areaTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
       
        // Apply the reduced forces to all nodes.
        thrust::for_each(
            thrust::make_zip_iterator(//1st begin
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::make_zip_iterator(//1st end
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
            AddForceFunctor (
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
            
};