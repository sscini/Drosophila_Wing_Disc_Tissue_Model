#include <stdio.h>
#include "System.h"
#include "SystemStructures.h"
//#include "AreaTriangles.h" // not needed
//#include "BendingTriangles.h"
//#include "MemRepulsionSprings_universal.h"
//#include "MemRepulsionSprings_local.h"
//#include "MemRepulsionEnergy.h"
#include "LinearSprings.h"
#include "NodeAdvance.h"
#include "Storage.h"
#include "Utilities.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"
#include <bits/stdc++.h> // apparently very bad and should be removed. 
//#include "LineTensionSprings.h" // only needed for either boundary nodes or for septin ring. Can be removed. 
#include <math.h>
#include <list>
//#include "TurgorForce.h"
#include "LinearSpringsEnergy.h"
#include "StrainTensor.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include "gradientRelax.h"

// somehow the gradient is not being set in my version - Kevin

// Helper function to count elements greater than or equal to zero in a vector.
int count_bigger(const std::vector<int> &elems)
{
    return std::count_if(elems.begin(), elems.end(), [](int c)
                         { return c >= 0; });
}

// Constructor for the System class.
System::System() {};

// Print net force on nodes along a radial line (?  0) from disc center to boundary
void System::PrintForce() {
    // Copy device forces to host
    thrust::host_vector<double> h_fx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> h_fy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> h_fz = coordInfoVecs.nodeForceZ;
    thrust::host_vector<double> h_x  = coordInfoVecs.nodeLocX;
    thrust::host_vector<double> h_y  = coordInfoVecs.nodeLocY;

    const double desiredTheta = 0.0;    // along +x axis
    const double eps = 0.01;            // angular tolerance (rad)
    std::vector<std::pair<double,int>> picks;
    int N = static_cast<int>(h_x.size());
    for (int i = 0; i < N; ++i) {
        double r = std::hypot(h_x[i], h_y[i]);
        double theta = std::atan2(h_y[i], h_x[i]);
        double diff = std::fabs(theta - desiredTheta);
        if (diff > M_PI) diff = 2*M_PI - diff;
        if (diff < eps) picks.emplace_back(r, i);
    }
    std::sort(picks.begin(), picks.end());

    std::printf("   r      Fx       Fy       Fz\n");
    for (auto &pr : picks) {
        int idx = pr.second;
        std::printf("%6.3f  %7.3e  %7.3e  %7.3e\n",
                    pr.first,
                    h_fx[idx], h_fy[idx], h_fz[idx]);
    }
}

// Function to solve the forces in the system.
void System::Solve_Forces()
{

    // Reset all forces to zero.
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

    // Compute forces and energy due to linear springs.
    ComputeLinearSprings(
        generalParams,
        coordInfoVecs,
        linearSpringInfoVecs,
        ljInfoVecs);

    // Compute forces and energy due to area springs. Nav commented out to test Active shape programming mesh type 02/27/2025  . Put back in 03/23/25
//      	ComputeAreaTriangleSprings(
//      		generalParams,
//      		coordInfoVecs,
//      		areaTriangleInfoVecs);

    // Compute forces and energy due to turgor pressure springs. (nav - commenting these out for now for flat surface 5/29/24) nav reintroducing the turgor pressure because the eversion wing does have turgor pressure. 8/17/2024
    // ComputeTurgorSprings(
    // generalParams,
    // coordInfoVecs,
    // areaTriangleInfoVecs
    //);

    // Compute forces and energy due to bending springs. Turn this off 10/10/24
//      	ComputeCosTriangleSprings(
//      		generalParams,
//      		coordInfoVecs,
//      		bendingTriangleInfoVecs);

    // Compute forces and energy due to membrane repulsion springs.// Nav commented out to test Active shape programming mesh type 02/27/2025. PUt back in 03/23/25
//      	ComputeMemRepulsionSprings_local(
//      		coordInfoVecs,
//      		linearSpringInfoVecs,
//      		capsidInfoVecs,
//      		generalParams,
//      		auxVecs);

        ComputeVolume(
          generalParams,
          coordInfoVecs,
          linearSpringInfoVecs,
          ljInfoVecs,
          prismInfoVecs);
    // Compute forces and energy due to volume springs. //(nav - commenting these out for now for flat surface 5/29/24) Nav had uncommented but she's bringing the comment back because testing out Active shape mesh 02/27/25
      	ComputeVolumeSprings(
      		coordInfoVecs,
          linearSpringInfoVecs,
          capsidInfoVecs,
          generalParams,
          auxVecs);
    
    // Now print forces along the radial line
   // PrintForce();
};

// Function to solve the entire system.
void System::solveSystem()
{
    // 0. Basic timestep
    generalParams.dt = 0.01;

    // 1. Compute initial center (apical mean)
    double cx = 0, cy = 0, cz = 0;
    int count = 0;
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        cx += coordInfoVecs.nodeLocX[i];
        cy += coordInfoVecs.nodeLocY[i];
        cz += coordInfoVecs.nodeLocZ[i];
        count++;
    }
    generalParams.centerX = cx / count;
    generalParams.centerY = cy / count;
    generalParams.centerZ = cz / count;

    // 2. Initial force check
    Solve_Forces();

    // 3. Print initial VTK
    storage->print_VTK_File();

    // ============================================
    //              STRAIN TENSOR STAGES
    // ============================================

    int stages = generalParams.Tf;
    double frac = 1.0;   // full-field application per stage

    for (int stage = 0; stage < stages; stage++)
    {
        // -- load lambda values for this stage --
        generalParams.lambda_iso_center_outDV = generalParams.lambda_iso_center_outDV_v[stage];
        generalParams.lambda_iso_edge_outDV   = generalParams.lambda_iso_edge_outDV_v[stage];
        generalParams.lambda_aniso_center_outDV = generalParams.lambda_aniso_center_outDV_v[stage];
        generalParams.lambda_aniso_edge_outDV   = generalParams.lambda_aniso_edge_outDV_v[stage];

        generalParams.lambda_iso_center_inDV = generalParams.lambda_iso_center_inDV_v[stage];
        generalParams.lambda_iso_edge_inDV   = generalParams.lambda_iso_edge_inDV_v[stage];
        generalParams.lambda_aniso_center_inDV = generalParams.lambda_aniso_center_inDV_v[stage];
        generalParams.lambda_aniso_edge_inDV   = generalParams.lambda_aniso_edge_inDV_v[stage];

        // Build vertex-level lambda
        LambdaField lambda;
        StrainTensorGPU::buildVertexLambda(generalParams, coordInfoVecs, lambda, frac);

        // Update rest lengths (initial_length → final_length)
        int layerflag = 0;   // 0 = whole tissue
        StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs, generalParams, lambda, linearSpringInfoVecs, layerflag);

        // Relaxation parameters
        generalParams.tol = 1e-5;
        int Nsteps = 10000;

        // ============================================
        //     GRADIENT RELAXATION LOOP
        // ============================================

        for (int iter = 0; iter < Nsteps; iter++)
        {
            // linearly increment rest lengths if doing time sweep
            for (int e = 0; e < coordInfoVecs.num_edges; e++) {
                double dl = (linearSpringInfoVecs.edge_final_length[e] -
                             linearSpringInfoVecs.edge_initial_length[e]) / double(Nsteps);
                linearSpringInfoVecs.edge_rest_length[e] += dl;
            }

            int k = relaxUntilConverged(*this);

            double E = linearSpringInfoVecs.linear_spring_energy;
            std::cout << "Relax iter " << iter 
                      << " Stage " << stage 
                      << " | E = " << E 
                      << " | Mov = " << generalParams.dx 
                      << " | Steps = " << k << std::endl;

            if (iter % 10 == 0)
                storage->print_VTK_File();
        }

        storage->print_VTK_File();
    }
}




// Function to assign the shared pointer to storage.
void System::assignStorage(std::shared_ptr<Storage> _storage)
{
    storage = _storage;
};

// Function to set the weak pointer to the SystemBuilder.
void System::set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr)
{
    weak_bld_ptr = _weak_bld_ptr;
};

// Function to initialize memory for thrust vectors and set coordInfoVecs values from input.
void System::initializeSystem(HostSetInfoVecs & hostSetInfoVecs)
{
    std::cout << "Initializing" << std::endl;

    // Set the max node count, edge count and triangle count.
    generalParams.maxNodeCount = hostSetInfoVecs.nodeLocX.size();
    coordInfoVecs.num_edges = hostSetInfoVecs.edges2Nodes_1.size();
    coordInfoVecs.num_triangles = hostSetInfoVecs.triangles2Nodes_1.size();

    std::cout << "num nodes: " << generalParams.maxNodeCount << std::endl;
    std::cout << "num edges: " << coordInfoVecs.num_edges << std::endl;
    std::cout << "num elems: " << coordInfoVecs.num_triangles << std::endl;
    // Allocate memory for various vectors using preallocated memory size.
    int mem_prealloc = 1;
    coordInfoVecs.scaling_per_edge.resize(mem_prealloc * coordInfoVecs.num_edges);
    // simplest safe default: all ones
    thrust::fill(coordInfoVecs.scaling_per_edge.begin(),
                 coordInfoVecs.scaling_per_edge.end(), 1.0);

   
    coordInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), false);
    coordInfoVecs.prevNodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    coordInfoVecs.prevNodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    coordInfoVecs.nodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.nodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.nodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    // coordInfoVecs.nodeVelX.resize(mem_prealloc*hostSetInfoVecs.nodeVelX.size(), 0.0);
    // coordInfoVecs.nodeVelY.resize(mem_prealloc*hostSetInfoVecs.nodeVelY.size(), 0.0);
    // coordInfoVecs.nodeVelZ.resize(mem_prealloc*hostSetInfoVecs.nodeVelZ.size(), 0.0);

    coordInfoVecs.nodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
    coordInfoVecs.nodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
    coordInfoVecs.nodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);

    coordInfoVecs.triangles2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Nodes_3.resize(mem_prealloc * coordInfoVecs.num_triangles);

    coordInfoVecs.triangles2Edges_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Edges_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Edges_3.resize(mem_prealloc * coordInfoVecs.num_triangles);

//    coordInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    coordInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    coordInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);

    hostSetInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    hostSetInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    hostSetInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);

    coordInfoVecs.edges2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_edges);
    coordInfoVecs.edges2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_edges);

    coordInfoVecs.edges2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_edges);
    coordInfoVecs.edges2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_edges);

    coordInfoVecs.SurfaceNormalX.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.SurfaceNormalY.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.SurfaceNormalZ.resize(mem_prealloc * generalParams.maxNodeCount);

    generalParams.nodes_in_upperhem.resize(mem_prealloc * generalParams.maxNodeCount);
    generalParams.triangles_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_triangles);
    generalParams.edges_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges);
    generalParams.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
    generalParams.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
    generalParams.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);

    hostSetInfoVecs.nodes_in_upperhem.resize(generalParams.nodes_in_upperhem.size());
    generalParams.nodes_in_upperhem = hostSetInfoVecs.nodes_in_upperhem;
    hostSetInfoVecs.triangles_in_upperhem.resize(generalParams.triangles_in_upperhem.size());
    hostSetInfoVecs.edges_in_upperhem.resize(generalParams.edges_in_upperhem.size());
    generalParams.edges_in_upperhem = hostSetInfoVecs.edges_in_upperhem;
    hostSetInfoVecs.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
    hostSetInfoVecs.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
    hostSetInfoVecs.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);



  
    // copy info to GPU
    std::cout << "Copying" << std::endl;
    thrust::copy(hostSetInfoVecs.isNodeFixed.begin(), hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());

    // Print information about fixed nodes in hostSetInfoVecs and coordInfoVecs.
    std::cout << "fixed_node_in_host: " << std::endl;
    for (int k = 0; k < hostSetInfoVecs.isNodeFixed.size(); k++)
    {
    }
    std::cout << "end_of_fixed_node_host_printout" << std::endl;
    std::cout << "fixed_node_in_device: " << std::endl;
    for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++)
    {
    }
    std::cout << "end_of_fixed_node_device_printout" << std::endl;
    std::cout << "size of host fixed " << hostSetInfoVecs.isNodeFixed.size() << std::endl;
    std::cout << "size of device fixed " << coordInfoVecs.isNodeFixed.size() << std::endl;

    // initialize various vectors with zeros or values from hostSetInfoVecs.
    //  Fill operations for other nodeForce and prevNodeForce vectors.
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
    
    // ==================== initialize lambda stage vectors ====================

    // outDV
    generalParams.lambda_iso_center_outDV_v.resize(3);
    generalParams.lambda_iso_edge_outDV_v.resize(3);
    generalParams.lambda_aniso_center_outDV_v.resize(3);
    generalParams.lambda_aniso_edge_outDV_v.resize(3);
    
    // inDV
    generalParams.lambda_iso_center_inDV_v.resize(3);
    generalParams.lambda_iso_edge_inDV_v.resize(3);
    generalParams.lambda_aniso_center_inDV_v.resize(3);
    generalParams.lambda_aniso_edge_inDV_v.resize(3);
    
    // temporary: no growth (all ones)
    thrust::fill(generalParams.lambda_iso_center_outDV_v.begin(),
                 generalParams.lambda_iso_center_outDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_iso_edge_outDV_v.begin(),
                 generalParams.lambda_iso_edge_outDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_aniso_center_outDV_v.begin(),
                 generalParams.lambda_aniso_center_outDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_aniso_edge_outDV_v.begin(),
                 generalParams.lambda_aniso_edge_outDV_v.end(), 1.0);
    
    thrust::fill(generalParams.lambda_iso_center_inDV_v.begin(),
                 generalParams.lambda_iso_center_inDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_iso_edge_inDV_v.begin(),
                 generalParams.lambda_iso_edge_inDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_aniso_center_inDV_v.begin(),
                 generalParams.lambda_aniso_center_inDV_v.end(), 1.0);
    thrust::fill(generalParams.lambda_aniso_edge_inDV_v.begin(),
                 generalParams.lambda_aniso_edge_inDV_v.end(), 1.0);
    
    // ========================================================================

    thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);

    // Copy node locations and other related information from hostSetInfoVecs to coordInfoVecs and other copy operations for triangles, edges and other related vectors.
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.prevNodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.prevNodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.nodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(), hostSetInfoVecs.triangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(), hostSetInfoVecs.triangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(), hostSetInfoVecs.triangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin());

    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(), hostSetInfoVecs.triangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(), hostSetInfoVecs.triangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(), hostSetInfoVecs.triangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin());

    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(), hostSetInfoVecs.edges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(), hostSetInfoVecs.edges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin());

    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(), hostSetInfoVecs.edges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(), hostSetInfoVecs.edges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin());

    thrust::copy(hostSetInfoVecs.nndata1.begin(), hostSetInfoVecs.nndata1.end(), coordInfoVecs.nndata1.begin());
    thrust::copy(hostSetInfoVecs.nndata2.begin(), hostSetInfoVecs.nndata2.end(), coordInfoVecs.nndata2.begin());
    thrust::copy(hostSetInfoVecs.nndata3.begin(), hostSetInfoVecs.nndata3.end(), coordInfoVecs.nndata3.begin());
    thrust::copy(hostSetInfoVecs.nndata4.begin(), hostSetInfoVecs.nndata4.end(), coordInfoVecs.nndata4.begin());
    thrust::copy(hostSetInfoVecs.nndata5.begin(), hostSetInfoVecs.nndata5.end(), coordInfoVecs.nndata5.begin());
    thrust::copy(hostSetInfoVecs.nndata6.begin(), hostSetInfoVecs.nndata6.end(), coordInfoVecs.nndata6.begin());
    thrust::copy(hostSetInfoVecs.nndata7.begin(), hostSetInfoVecs.nndata7.end(), coordInfoVecs.nndata7.begin());
    thrust::copy(hostSetInfoVecs.nndata8.begin(), hostSetInfoVecs.nndata8.end(), coordInfoVecs.nndata8.begin());
    thrust::copy(hostSetInfoVecs.nndata9.begin(), hostSetInfoVecs.nndata9.end(), coordInfoVecs.nndata9.begin());

    // Resize and initialize the 'u' vector.
    coordInfoVecs.u.resize(mem_prealloc * coordInfoVecs.num_triangles);

    // Part 18

    // Allocate memory for additiional data structures.

    // Area triangle info vec.
    // Number of area springs is the number of triangles
    std::cout << "Mem" << std::endl;
        
    // linear springs info vectors.
    //  Allocate memory for temporary node information in unreduced form for linear springs.
    linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);

    // Allocate memory for temporary node information in reduced form for bending springs.
    linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);

    // Clear edge_initial_length vector for linear springs.
    // linearSpringInfoVecs.edge_initial_length.clear();
    // linearSpringInfoVecs.edge_rest_length.clear();

    // linearSpringInfoVecs.edge_rest_length.resize(hostSetInfoVecs.edge_rest_length.size());
    linearSpringInfoVecs.edge_final_length.resize(hostSetInfoVecs.edge_initial_length.size());
    linearSpringInfoVecs.edge_initial_length = hostSetInfoVecs.edge_initial_length;
    linearSpringInfoVecs.edge_final_length = hostSetInfoVecs.edge_initial_length;
    
    
    std::cout << "host edge_initial_length size = " << hostSetInfoVecs.edge_initial_length.size() << std::endl;
    std::cout << "device edge_initial_length size = " << linearSpringInfoVecs.edge_initial_length.size() << std::endl;

    //  for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
    //    int i = coordInfoVecs.edges2Nodes_1[e];
    //    int j = coordInfoVecs.edges2Nodes_2[e];
    //    double dx = hostSetInfoVecs.nodeLocX[j] - hostSetInfoVecs.nodeLocX[i];
    //    double dy = hostSetInfoVecs.nodeLocY[j] - hostSetInfoVecs.nodeLocY[i];
    //    double dz = hostSetInfoVecs.nodeLocZ[j] - hostSetInfoVecs.nodeLocZ[i];
    //    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    //    //hostSetInfoVecs.edge_initial_length.push_back(dist);    // already done for initial
    //    hostSetInfoVecs.edge_rest_length.push_back(dist);
    //    //std::cout<< "edge_rest_length = " << hostSetInfoVecs.edge_initial_length[i]<<std::endl; in the current data structure it gave me 1021 edges. That's good. Now that they have been initialized I should start changing the rest lengths.
    //    }

    // thrust::copy(linearSpringInfoVecs.edge_rest_length.begin(),
    //              linearSpringInfoVecs.edge_rest_length.end(),
    //              hostSetInfoVecs.edge_rest_length.begin());

    linearSpringInfoVecs.edge_rest_length = hostSetInfoVecs.edge_initial_length;

    //linearSpringInfoVecs.edge_rest_length.resize(1/(generalParams.dt*generalParams.tol)) //= hostSetInfoVecs.edge_rest_length;

    std::cout << "host edge_rest_length size = " << hostSetInfoVecs.edge_rest_length.size() << std::endl;
    std::cout << "device edge_rest_length size = " << linearSpringInfoVecs.edge_rest_length.size() << std::endl;

    //linearSpringInfoVecs.edge_final_length.resize(coordInfoVecs.num_edges);

    //linearSpringInfoVecs.edge_final_length = linearSpringInfoVecs.edge_initial_length;

    // This loop is to test out whether the edge_rest_length is being initialized properly.
    // for (int i = 0; i < hostSetInfoVecs.edge_rest_length.size(); i++) {
    //    std::cout<< "edge_rest_length # "<< i <<" = "<< hostSetInfoVecs.edge_rest_length[i]<<std::endl;
    //    std::cout<< "edge_initial_length # "<< i <<" = "<< hostSetInfoVecs.edge_initial_length[i]<<std::endl;
    //}

    // Resize the hostSetInfoVecs for data transfer between host and device.
    hostSetInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    hostSetInfoVecs.nodeLocX.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeLocY.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeLocZ.resize(coordInfoVecs.nodeLocX.size());
    std::cout << "Host_nodeLocX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;

    // hostSetInfoVecs.nodeVelX.resize(coordInfoVecs.nodeVelX.size());
    // hostSetInfoVecs.nodeVelY.resize(coordInfoVecs.nodeVelY.size());
    // hostSetInfoVecs.nodeVelZ.resize(coordInfoVecs.nodeVelZ.size());

    hostSetInfoVecs.nodeForceX.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeForceY.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeForceZ.resize(coordInfoVecs.nodeLocX.size());
    std::cout << "Host_nodeForceX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;

    hostSetInfoVecs.triangles2Nodes_1.resize(coordInfoVecs.triangles2Nodes_1.size());
    hostSetInfoVecs.triangles2Nodes_2.resize(coordInfoVecs.triangles2Nodes_2.size());
    hostSetInfoVecs.triangles2Nodes_3.resize(coordInfoVecs.triangles2Nodes_3.size());
    std::cout << "Host_triangles2Nodes size = " << hostSetInfoVecs.triangles2Nodes_1.size() << std::endl;

    hostSetInfoVecs.triangles2Edges_1.resize(coordInfoVecs.triangles2Edges_1.size());
    hostSetInfoVecs.triangles2Edges_2.resize(coordInfoVecs.triangles2Edges_2.size());
    hostSetInfoVecs.triangles2Edges_3.resize(coordInfoVecs.triangles2Edges_3.size());
    std::cout << "Host_triangles2Edges size = " << hostSetInfoVecs.triangles2Edges_1.size() << std::endl;

    hostSetInfoVecs.edges2Nodes_1.resize(coordInfoVecs.edges2Nodes_1.size());
    hostSetInfoVecs.edges2Nodes_2.resize(coordInfoVecs.edges2Nodes_2.size());
    std::cout << "Host_edges2Nodes size = " << hostSetInfoVecs.edges2Nodes_1.size() << std::endl;

    hostSetInfoVecs.edges2Triangles_1.resize(coordInfoVecs.edges2Triangles_1.size());
    hostSetInfoVecs.edges2Triangles_2.resize(coordInfoVecs.edges2Triangles_2.size());
    std::cout << "Host_edges2Triangles size = " << hostSetInfoVecs.edges2Triangles_1.size() << std::endl;

// Print message indicating the system is ready.
    std::cout << "System Ready" << std::endl;

    // Allocate memory for buckets.
    auxVecs.id_bucket.resize(generalParams.maxNodeCount);
    auxVecs.id_value.resize(generalParams.maxNodeCount);
    auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
    auxVecs.id_value_expanded.resize(27 * (generalParams.maxNodeCount));
};






