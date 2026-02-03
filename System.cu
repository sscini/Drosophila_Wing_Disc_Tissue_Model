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

#include <thrust/host_vector.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>

// ============================================================================
// Helper function: sixV_tet_host
// Computes 6 * signed volume of tetrahedron with vertices i,j,k,l
// This should already exist in your System.cu, included here for completeness
// ============================================================================

static inline double sixV_tet_host(
    const std::vector<double>& X,
    const std::vector<double>& Y,
    const std::vector<double>& Z,
    int i, int j, int k, int l)
{
    const double ux = X[i] - X[l];
    const double uy = Y[i] - Y[l];
    const double uz = Z[i] - Z[l];

    const double vx = X[j] - X[l];
    const double vy = Y[j] - Y[l];
    const double vz = Z[j] - Z[l];

    const double wx = X[k] - X[l];
    const double wy = Y[k] - Y[l];
    const double wz = Z[k] - Z[l];

    // Cross product v × w
    const double cx = vy * wz - vz * wy;
    const double cy = vz * wx - vx * wz;
    const double cz = vx * wy - vy * wx;

    // Dot product u · (v × w) = 6 * signed volume
    return ux * cx + uy * cy + uz * cz;
}

// Compute volume of a prism using 3-tet decomposition
// Returns the 6*V value (divide by 6 to get actual volume)
static double computePrism6V(
    const std::vector<double>& X,
    const std::vector<double>& Y,
    const std::vector<double>& Z,
    int a, int b, int c,  // Top triangle (apical)
    int A, int B, int C)  // Bottom triangle (basal) - A below a, B below b, C below c
{
    // Standard 3-tet decomposition with apex at 'a':
    // Tet 1: (b, c, A, a)
    // Tet 2: (b, A, C, a)  
    // Tet 3: (A, B, C, a)
    double s1 = sixV_tet_host(X, Y, Z, b, c, A, a);
    double s2 = sixV_tet_host(X, Y, Z, b, A, C, a);
    double s3 = sixV_tet_host(X, Y, Z, A, B, C, a);
    
    return s1 + s2 + s3;
}

// ============================================================================
// FIXED BuildPrismsFromLayerFlags for Dome/Curved Geometry
// ============================================================================
// 
// Problem: The previous version used geometric filtering (XY vs Z displacement)
// to identify "straight-down" edges. This FAILS for dome geometry where the
// vertical edges are radial, not Cartesian-vertical.
//
// Solution: Use INDEX-BASED CORRESPONDENCE. Your mesh is constructed such that:
//   - Apical nodes are indices 0 to (A-1)
//   - Basal nodes are indices A to (2A-1)  [for 2-layer case]
//   - The straight-down edge from apical node i goes to basal node (i + A)
//
// This works because XMLParser.cu adds nodes in order:
//   1. Apical nodes (with layer flag = maxLayer)
//   2. Body nodes if any (with layer flags 1 to maxLayer-1)
//   3. Basal nodes (with layer flag = 0)
//
// ============================================================================

static void BuildPrismsFromLayerFlags(
    const GeneralParams& gp,
    const CoordInfoVecs& coord,
    PrismInfoVecs& prism)
{
    thrust::host_vector<int> hLayer = gp.nodes_in_upperhem;
    thrust::host_vector<int> t1 = coord.triangles2Nodes_1;
    thrust::host_vector<int> t2 = coord.triangles2Nodes_2;
    thrust::host_vector<int> t3 = coord.triangles2Nodes_3;

    const int Ntot = (int)hLayer.size();
    if (Ntot == 0) {
        std::cout << "ERROR: nodes_in_upperhem is empty.\n";
        prism.num_prisms = 0;
        return;
    }

    // Copy coordinates to host for volume calculations
    const int Nnodes = (int)coord.nodeLocX.size();
    std::vector<double> X(Nnodes), Y(Nnodes), Z(Nnodes);
    for (int i = 0; i < Nnodes; ++i) {
        X[i] = coord.nodeLocX[i];
        Y[i] = coord.nodeLocY[i];
        Z[i] = coord.nodeLocZ[i];
    }

    // ========================================================================
    // STEP 1: Determine layer structure
    // ========================================================================
    
    int maxFlag = -1;
    int minFlag = INT_MAX;
    for (int i = 0; i < Ntot; ++i) {
        if (hLayer[i] >= 0) {
            maxFlag = std::max(maxFlag, hLayer[i]);
            minFlag = std::min(minFlag, hLayer[i]);
        }
    }

    if (maxFlag < 1) {
        std::cout << "ERROR: Need at least two layers (maxFlag >= 1).\n";
        prism.num_prisms = 0;
        return;
    }

    std::cout << "BuildPrisms: Layer flags range from " << minFlag << " to " << maxFlag << std::endl;
    std::cout << "BuildPrisms: Apical layer = " << maxFlag << ", Basal layer = " << minFlag << std::endl;

    // ========================================================================
    // STEP 2: Build layer membership lists
    // Collect which nodes belong to which layer
    // ========================================================================
    
    std::vector<std::vector<int>> layers(maxFlag + 1);
    for (int i = 0; i < Ntot; ++i) {
        int ell = hLayer[i];
        if (ell >= 0 && ell <= maxFlag) {
            layers[ell].push_back(i);
        }
    }

    // Print layer sizes for debugging
    std::cout << "BuildPrisms: Layer sizes:" << std::endl;
    for (int ell = 0; ell <= maxFlag; ++ell) {
        std::cout << "  Layer " << ell << ": " << layers[ell].size() << " nodes" << std::endl;
    }

    // ========================================================================
    // STEP 3: Validate that adjacent layers have same node count
    // This is required for index-based correspondence
    // ========================================================================
    
    for (int ell = 1; ell <= maxFlag; ++ell) {
        if (layers[ell].empty() || layers[ell-1].empty()) {
            std::cout << "ERROR: Layer " << ell << " or " << (ell-1) << " is empty.\n";
            prism.num_prisms = 0;
            return;
        }
        if (layers[ell].size() != layers[ell-1].size()) {
            std::cout << "ERROR: Layer " << ell << " has " << layers[ell].size()
                      << " nodes but layer " << (ell-1) << " has " << layers[ell-1].size()
                      << ". Cannot use index-based correspondence.\n";
            prism.num_prisms = 0;
            return;
        }
    }

    // ========================================================================
    // STEP 4: Build downMap using INDEX-BASED CORRESPONDENCE
    // 
    // Key insight: Nodes are added in consistent order by XMLParser:
    //   - First all apical nodes (layer maxFlag)
    //   - Then body nodes if any (layers 1 to maxFlag-1)
    //   - Then basal nodes (layer 0)
    //
    // Within each layer, nodes maintain the same local ordering.
    // So node k in layer ell corresponds to node k in layer ell-1.
    // ========================================================================
    
    std::vector<int> downMap(Nnodes, -1);
    
    // Build local index within each layer
    std::vector<int> localIndex(Nnodes, -1);
    for (int ell = 0; ell <= maxFlag; ++ell) {
        for (int k = 0; k < (int)layers[ell].size(); ++k) {
            int globalNodeId = layers[ell][k];
            localIndex[globalNodeId] = k;
        }
    }

    // Now build downMap: for each node in layer ell > 0, 
    // find the node in layer ell-1 with the same local index
    int mappedCount = 0;
    for (int ell = 1; ell <= maxFlag; ++ell) {
        for (int k = 0; k < (int)layers[ell].size(); ++k) {
            int upperNode = layers[ell][k];
            int lowerNode = layers[ell-1][k];  // Same local index k
            downMap[upperNode] = lowerNode;
            mappedCount++;
        }
    }

    std::cout << "BuildPrisms: Created " << mappedCount 
              << " node mappings using index-based correspondence" << std::endl;

    // Debug: Print first few mappings
    std::cout << "BuildPrisms: Sample mappings (upper -> lower):" << std::endl;
    int printCount = 0;
    for (int i = 0; i < Nnodes && printCount < 10; ++i) {
        if (downMap[i] >= 0) {
            std::cout << "  Node " << i << " (layer " << hLayer[i] << ") -> Node " 
                      << downMap[i] << " (layer " << hLayer[downMap[i]] << ")" << std::endl;
            printCount++;
        }
    }

    // ========================================================================
    // STEP 5: Build prisms from triangles in upper layers
    // ========================================================================
    
    auto downNode = [&](int node) -> int {
        if (node < 0 || node >= Nnodes) return -1;
        return downMap[node];
    };

    std::vector<int> P1, P2, P3, P4, P5, P6;
    const int T = coord.num_triangles;
    P1.reserve(T); P2.reserve(T); P3.reserve(T);
    P4.reserve(T); P5.reserve(T); P6.reserve(T);

    int flipped = 0;
    int skipped_no_mapping = 0;
    int skipped_degenerate = 0;
    int skipped_wrong_layer = 0;

    for (int ti = 0; ti < T; ++ti) {
        int a = t1[ti], b = t2[ti], c = t3[ti];
        
        // Validate node indices
        if (a < 0 || b < 0 || c < 0) continue;
        if (a >= Nnodes || b >= Nnodes || c >= Nnodes) continue;

        int la = hLayer[a], lb = hLayer[b], lc = hLayer[c];

        // Only process triangles fully in one layer that has a layer below
        if (la < 1) {
            // Layer 0 (basal) has no lower neighbor - skip
            continue;
        }
        
        if (la != lb || la != lc) {
            // Triangle spans multiple layers - not a sheet triangle
            skipped_wrong_layer++;
            continue;
        }

        // Get corresponding nodes in the layer below
        int A = downNode(a);
        int B = downNode(b);
        int C = downNode(c);

        if (A < 0 || B < 0 || C < 0) {
            skipped_no_mapping++;
            continue;
        }

        // Verify the lower nodes are actually in layer la-1
        if (hLayer[A] != la-1 || hLayer[B] != la-1 || hLayer[C] != la-1) {
            std::cout << "WARNING: Prism mapping error for triangle " << ti 
                      << ": lower nodes not in expected layer" << std::endl;
            skipped_no_mapping++;
            continue;
        }

        // ====================================================================
        // Compute signed volume using 3-tetrahedron decomposition
        // Same decomposition as in VolumeComp.cu:
        //   s1 = 6V_tet(P2,P3,P4,P1) = 6V_tet(b,c,A,a)
        //   s2 = 6V_tet(P2,P4,P6,P1) = 6V_tet(b,A,C,a)
        //   s3 = 6V_tet(P4,P5,P6,P1) = 6V_tet(A,B,C,a)
        // ====================================================================
        
        double s1 = sixV_tet_host(X, Y, Z, b, c, A, a);
        double s2 = sixV_tet_host(X, Y, Z, b, A, C, a);
        double s3 = sixV_tet_host(X, Y, Z, A, B, C, a);
        double sum6V = s1 + s2 + s3;

        // Skip degenerate prisms (near-zero volume)
        if (std::fabs(sum6V) < 1e-12) {
            skipped_degenerate++;
            continue;
        }

        // Ensure consistent orientation - make all prisms have positive volume
        if (sum6V < 0.0) {
            std::swap(b, c);
            std::swap(B, C);
            flipped++;
        }

        // Add the prism
        P1.push_back(a); P2.push_back(b); P3.push_back(c);
        P4.push_back(A); P5.push_back(B); P6.push_back(C);
    }

    // Store results
    prism.num_prisms = (int)P1.size();
    prism.P1 = P1; prism.P2 = P2; prism.P3 = P3;
    prism.P4 = P4; prism.P5 = P5; prism.P6 = P6;

    // ========================================================================
    // STEP 6: Report results and validate
    // ========================================================================
    
    std::cout << "BuildPrisms: Created " << prism.num_prisms << " prisms" << std::endl;
    std::cout << "  - Flipped for consistent orientation: " << flipped << std::endl;
    std::cout << "  - Skipped (no mapping): " << skipped_no_mapping << std::endl;
    std::cout << "  - Skipped (degenerate): " << skipped_degenerate << std::endl;
    std::cout << "  - Skipped (wrong layer): " << skipped_wrong_layer << std::endl;

    // Validate: Check for negative volumes in final prisms
    int negativeCount = 0;
    double totalVolume = 0.0;
    
    for (int p = 0; p < prism.num_prisms; ++p) {
        double s1 = sixV_tet_host(X, Y, Z, P2[p], P3[p], P4[p], P1[p]);
        double s2 = sixV_tet_host(X, Y, Z, P2[p], P4[p], P6[p], P1[p]);
        double s3 = sixV_tet_host(X, Y, Z, P4[p], P5[p], P6[p], P1[p]);
        double vol = (s1 + s2 + s3) / 6.0;
        totalVolume += vol;

        if (vol < 0) {
            if (negativeCount < 10) {
                std::cout << "  WARNING: Prism " << p << " has negative volume = " 
                          << vol << std::endl;
            }
            negativeCount++;
        }
    }

    if (negativeCount > 0) {
        std::cout << "  WARNING: " << negativeCount << " prisms have negative volume!" << std::endl;
    } else {
        std::cout << "  SUCCESS: All prisms have positive volume" << std::endl;
    }
    
    std::cout << "  Total prism volume: " << totalVolume << std::endl;
}




// Diagnostic code to add to your System.cu after BuildPrismsFromLayerFlags
// This will print detailed information about each prism to help debug

void DiagnosePrisms(
    const GeneralParams& gp,
    const CoordInfoVecs& coord,
    const PrismInfoVecs& prism)
{
    const int P = prism.num_prisms;
    if (P <= 0) {
        std::cout << "DiagnosePrisms: No prisms!" << std::endl;
        return;
    }

    const int N = (int)coord.nodeLocX.size();
    
    // Copy to host
    thrust::host_vector<int> hP1 = prism.P1;
    thrust::host_vector<int> hP2 = prism.P2;
    thrust::host_vector<int> hP3 = prism.P3;
    thrust::host_vector<int> hP4 = prism.P4;
    thrust::host_vector<int> hP5 = prism.P5;
    thrust::host_vector<int> hP6 = prism.P6;
    thrust::host_vector<int> hLayer = gp.nodes_in_upperhem;

    std::vector<double> X(N), Y(N), Z(N);
    for (int i = 0; i < N; ++i) {
        X[i] = coord.nodeLocX[i];
        Y[i] = coord.nodeLocY[i];
        Z[i] = coord.nodeLocZ[i];
    }

    std::cout << "\n========== PRISM DIAGNOSTIC ==========" << std::endl;
    std::cout << "Total prisms: " << P << std::endl;
    std::cout << "Total nodes: " << N << std::endl;

    // Print all node positions and layers
    std::cout << "\nNode positions:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "  Node " << i << ": (" << X[i] << ", " << Y[i] << ", " << Z[i] 
                  << ") layer=" << hLayer[i] << std::endl;
    }

    std::cout << "\nPrism details:" << std::endl;
    
    double totalVol = 0.0;
    
    for (int p = 0; p < P; ++p) {
        int a = hP1[p], b = hP2[p], c = hP3[p];
        int A = hP4[p], B = hP5[p], C = hP6[p];

        std::cout << "\nPrism " << p << ":" << std::endl;
        std::cout << "  Top triangle: " << a << ", " << b << ", " << c << std::endl;
        std::cout << "  Bottom triangle: " << A << ", " << B << ", " << C << std::endl;
        
        std::cout << "  Top node positions:" << std::endl;
        std::cout << "    " << a << ": (" << X[a] << ", " << Y[a] << ", " << Z[a] << ") L=" << hLayer[a] << std::endl;
        std::cout << "    " << b << ": (" << X[b] << ", " << Y[b] << ", " << Z[b] << ") L=" << hLayer[b] << std::endl;
        std::cout << "    " << c << ": (" << X[c] << ", " << Y[c] << ", " << Z[c] << ") L=" << hLayer[c] << std::endl;
        
        std::cout << "  Bottom node positions:" << std::endl;
        std::cout << "    " << A << ": (" << X[A] << ", " << Y[A] << ", " << Z[A] << ") L=" << hLayer[A] << std::endl;
        std::cout << "    " << B << ": (" << X[B] << ", " << Y[B] << ", " << Z[B] << ") L=" << hLayer[B] << std::endl;
        std::cout << "    " << C << ": (" << X[C] << ", " << Y[C] << ", " << Z[C] << ") L=" << hLayer[C] << std::endl;

        // Check vertical alignment
        double dist_aA = sqrt(pow(X[a]-X[A],2) + pow(Y[a]-Y[A],2));
        double dist_bB = sqrt(pow(X[b]-X[B],2) + pow(Y[b]-Y[B],2));
        double dist_cC = sqrt(pow(X[c]-X[C],2) + pow(Y[c]-Y[C],2));
        
        std::cout << "  XY distances (should be ~0 for vertical alignment):" << std::endl;
        std::cout << "    a-A: " << dist_aA << std::endl;
        std::cout << "    b-B: " << dist_bB << std::endl;
        std::cout << "    c-C: " << dist_cC << std::endl;

        // Compute tetrahedra volumes
        auto sixV = [&](int i, int j, int k, int l) {
            double ux = X[i] - X[l], uy = Y[i] - Y[l], uz = Z[i] - Z[l];
            double vx = X[j] - X[l], vy = Y[j] - Y[l], vz = Z[j] - Z[l];
            double wx = X[k] - X[l], wy = Y[k] - Y[l], wz = Z[k] - Z[l];
            double cx = vy * wz - vz * wy;
            double cy = vz * wx - vx * wz;
            double cz = vx * wy - vy * wx;
            return ux * cx + uy * cy + uz * cz;
        };

        double s1 = sixV(b, c, A, a);
        double s2 = sixV(b, A, C, a);
        double s3 = sixV(A, B, C, a);
        double prismVol = (s1 + s2 + s3) / 6.0;

        std::cout << "  Tetrahedra 6V values:" << std::endl;
        std::cout << "    Tet1 (b,c,A,a): " << s1 << " -> V=" << s1/6.0 << std::endl;
        std::cout << "    Tet2 (b,A,C,a): " << s2 << " -> V=" << s2/6.0 << std::endl;
        std::cout << "    Tet3 (A,B,C,a): " << s3 << " -> V=" << s3/6.0 << std::endl;
        std::cout << "  Prism volume: " << prismVol << std::endl;

        totalVol += prismVol;
    }

    std::cout << "\nTotal volume from prisms: " << totalVol << std::endl;
    std::cout << "========================================\n" << std::endl;
}

// Call this right after BuildPrismsFromLayerFlags:
// DiagnosePrisms(generalParams, coordInfoVecs, prismInfoVecs);

static inline int clampN(int n, int maxn) { return std::min(std::max(n, 0), maxn); }

void DebugPrintMeshConnectivity(System& sys,
                                int maxTriPrint  = 50,
                                int maxEdgePrint = 100,
                                bool printAll    = false) {
    auto& coord  = sys.coordInfoVecs;
    auto& gp     = sys.generalParams;
    auto& spring = sys.linearSpringInfoVecs;

    // -------- triangles (elems) ----------
    const int T = coord.num_triangles;
    thrust::host_vector<int> h_t1 = coord.triangles2Nodes_1;
    thrust::host_vector<int> h_t2 = coord.triangles2Nodes_2;
    thrust::host_vector<int> h_t3 = coord.triangles2Nodes_3;

    std::cout << "\n===== TRIANGLES (" << T << ") =====\n";
    int triN = printAll ? T : clampN(maxTriPrint, T);
    for (int t = 0; t < triN; ++t) {
        std::cout << "T" << t << ": " << h_t1[t] << " " << h_t2[t] << " " << h_t3[t] << "\n";
    }

    // -------- edges (edgeinfos) ----------
    const int E = coord.num_edges;
    thrust::host_vector<int> h_e1 = coord.edges2Nodes_1;
    thrust::host_vector<int> h_e2 = coord.edges2Nodes_2;

    // Edge layer flags (you store these in generalParams.edges_in_upperhem)
    thrust::host_vector<int> h_edgeLayer;
    if ((int)gp.edges_in_upperhem.size() == E) {
        h_edgeLayer = gp.edges_in_upperhem;
    }

    // Node layer flags (optional, useful to sanity-check vertical edges)
    thrust::host_vector<int> h_nodeLayer;
    if ((int)gp.nodes_in_upperhem.size() == (int)coord.nodeLocX.size()) {
        h_nodeLayer = gp.nodes_in_upperhem;
    }

    // Edge rest length (if present)
    thrust::host_vector<double> h_rest;
    bool hasRest = ((int)spring.edge_rest_length.size() == E);
    if (hasRest) h_rest = spring.edge_rest_length;

    std::cout << "\n===== EDGES (" << E << ") =====\n";
    int edgeN = printAll ? E : clampN(maxEdgePrint, E);
    for (int e = 0; e < edgeN; ++e) {
        int a = h_e1[e];
        int b = h_e2[e];

        std::cout << "E" << e << ": " << a << " " << b;

        if (hasRest) std::cout << "   rest=" << h_rest[e];

        if ((int)h_edgeLayer.size() == E) std::cout << "   edgeLayer=" << h_edgeLayer[e];

        // Bonus: also print node layers so you can spot “vertical” edges instantly
        if ((int)h_nodeLayer.size() > std::max(a, b)) {
            std::cout << "   (nodeLayers " << h_nodeLayer[a] << "," << h_nodeLayer[b] << ")";
        }

        std::cout << "\n";
    }

    std::cout << std::endl;
}


// somehow the gradient is not being set in my version - Kevin

// Helper function to count elements greater than or equal to zero in a vector.
int count_bigger(const std::vector<int> &elems)
{
    return std::count_if(elems.begin(), elems.end(), [](int c)
                         { return c >= 0; });
}


// Constructor for the System class.
System::System() {};

void System::printEdges() {
    const auto& n1 = coordInfoVecs.edges2Nodes_1;
    const auto& n2 = coordInfoVecs.edges2Nodes_2;
    const auto& rest = linearSpringInfoVecs.edge_rest_length;
    const auto& layer = generalParams.edges_in_upperhem;

    int E = coordInfoVecs.num_edges;
    std::cout << "\n===== EDGES (" << E << ") =====\n";

    for (int i = 0; i < E; i++) {
        std::cout << "E" << i << ": "
                  << n1[i] << " "
                  << n2[i] << "   "
                  << rest[i] << "   "
                  << layer[i] << "\n";
    }
}

void System::printTriangles() {
    const auto& tri1 = coordInfoVecs.triangles2Nodes_1;
    const auto& tri2 = coordInfoVecs.triangles2Nodes_2;
    const auto& tri3 = coordInfoVecs.triangles2Nodes_3;

    int T = coordInfoVecs.num_triangles;
    std::cout << "\n===== TRIANGLES (" << T << ") =====\n";

    for (int i = 0; i < T; i++) {
        std::cout << "T" << i << ": "
                  << tri1[i] << " "
                  << tri2[i] << " "
                  << tri3[i] << "\n";
    }
}

// ============================================================================
// FORCE DECOMPOSITION DIAGNOSTIC
//
// This code verifies that: Total Force = Linear Spring Force + Volume Force
//
// Add this code inside your relaxation loop, after Solve_Forces() is called,
// or add it as a separate diagnostic function.
//
// The idea is to:
// 1. Save the total forces (which are computed by Solve_Forces)
// 2. Compute linear forces only and save them
// 3. Compute volume forces only and save them  
// 4. Verify that linear + volume = total
// ============================================================================

// OPTION 1: Add this as a diagnostic function in System.cu
// Call it once before the main relaxation loop to verify force computation


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

void DebugPrintLayerFlags(const GeneralParams& gp, const CoordInfoVecs& coord) {
    thrust::host_vector<int> hLayer = gp.nodes_in_upperhem;
    thrust::host_vector<int> eLayer = gp.edges_in_upperhem;
    
    std::map<int, int> nodeCounts;
    std::map<int, int> edgeCounts;
    
    for (int i = 0; i < hLayer.size(); ++i) {
        nodeCounts[hLayer[i]]++;
    }
    for (int i = 0; i < eLayer.size(); ++i) {
        edgeCounts[eLayer[i]]++;
    }
    
    std::cout << "=== LAYER FLAG DISTRIBUTION ===" << std::endl;
    std::cout << "Nodes by layer:" << std::endl;
    for (auto& kv : nodeCounts) {
        std::cout << "  Layer " << kv.first << ": " << kv.second << " nodes" << std::endl;
    }
    std::cout << "Edges by layer:" << std::endl;
    for (auto& kv : edgeCounts) {
        std::cout << "  Layer " << kv.first << ": " << kv.second << " edges" << std::endl;
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
        linearSpringInfoVecs);

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
          
          
//        {
//            thrust::host_vector<bool> h_fixed = coordInfoVecs.isNodeFixed;
//            int fixed = thrust::count(h_fixed.begin(), h_fixed.end(), true);
//            int total = h_fixed.size();
//            std::cout << "isNodeFixed: " << fixed << "/" << total << " nodes fixed ("
//                      << (100.0 * fixed / total) << "%)" << std::endl;
//            if (fixed == total) {
//                std::cout << "*** ERROR: ALL NODES ARE FIXED! ***" << std::endl;
//            }
//        }

          // Add this diagnostic code in Solve_Forces(), right after ComputeVolume():

        static int debug_call_count = 0;
        if (debug_call_count < 10 && generalParams.current_total_volume < 0) {
            std::cout << "\n=== NEGATIVE VOLUME DEBUG (call " << debug_call_count << ") ===" << std::endl;
            
            const int P = prismInfoVecs.num_prisms;
            const int N = (int)coordInfoVecs.nodeLocX.size();
            
            // Copy prism connectivity to host
            thrust::host_vector<int> hP1 = prismInfoVecs.P1;
            thrust::host_vector<int> hP2 = prismInfoVecs.P2;
            thrust::host_vector<int> hP3 = prismInfoVecs.P3;
            thrust::host_vector<int> hP4 = prismInfoVecs.P4;
            thrust::host_vector<int> hP5 = prismInfoVecs.P5;
            thrust::host_vector<int> hP6 = prismInfoVecs.P6;
            
            // Copy coordinates to host
            std::vector<double> X(N), Y(N), Z(N);
            for (int i = 0; i < N; ++i) {
                X[i] = coordInfoVecs.nodeLocX[i];
                Y[i] = coordInfoVecs.nodeLocY[i];
                Z[i] = coordInfoVecs.nodeLocZ[i];
            }
            
            // Lambda to compute 6*volume of tetrahedron
            auto sixV_tet = [&](int i, int j, int k, int l) -> double {
                double ux = X[i] - X[l], uy = Y[i] - Y[l], uz = Z[i] - Z[l];
                double vx = X[j] - X[l], vy = Y[j] - Y[l], vz = Z[j] - Z[l];
                double wx = X[k] - X[l], wy = Y[k] - Y[l], wz = Z[k] - Z[l];
                double cx = vy * wz - vz * wy;
                double cy = vz * wx - vx * wz;
                double cz = vx * wy - vy * wx;
                return ux * cx + uy * cy + uz * cz;
            };
            
            int neg_count = 0;
            double pos_sum = 0.0, neg_sum = 0.0;
            
            // Track first few negative prisms for detailed output
            std::vector<int> first_negative_prisms;
            
            for (int p = 0; p < P; ++p) {
                int a = hP1[p], b = hP2[p], c = hP3[p];
                int A = hP4[p], B = hP5[p], C = hP6[p];
                
                // Same decomposition as VolumeComp.cu
                double s1 = sixV_tet(b, c, A, a);
                double s2 = sixV_tet(b, A, C, a);
                double s3 = sixV_tet(A, B, C, a);
                double vol = (s1 + s2 + s3) / 6.0;
                
                if (vol < 0) {
                    neg_count++;
                    neg_sum += vol;
                    if (first_negative_prisms.size() < 5) {
                        first_negative_prisms.push_back(p);
                    }
                } else {
                    pos_sum += vol;
                }
            }
            
            std::cout << "Total prisms: " << P << std::endl;
            std::cout << "Positive prisms: " << (P - neg_count) << ", sum = " << pos_sum << std::endl;
            std::cout << "Negative prisms: " << neg_count << ", sum = " << neg_sum << std::endl;
            std::cout << "Net volume: " << (pos_sum + neg_sum) << std::endl;
            std::cout << "Reported volume: " << generalParams.current_total_volume << std::endl;
            
            // Print details of first few negative prisms
            if (!first_negative_prisms.empty()) {
                std::cout << "\nFirst negative prisms:" << std::endl;
                thrust::host_vector<int> hLayer = generalParams.nodes_in_upperhem;
                
                for (int p : first_negative_prisms) {
                    int a = hP1[p], b = hP2[p], c = hP3[p];
                    int A = hP4[p], B = hP5[p], C = hP6[p];
                    
                    double s1 = sixV_tet(b, c, A, a);
                    double s2 = sixV_tet(b, A, C, a);
                    double s3 = sixV_tet(A, B, C, a);
                    double vol = (s1 + s2 + s3) / 6.0;
                    
                    std::cout << "  Prism " << p << ": vol=" << vol << std::endl;
                    std::cout << "    Top (a,b,c): " << a << "," << b << "," << c 
                              << " layers: " << hLayer[a] << "," << hLayer[b] << "," << hLayer[c] << std::endl;
                    std::cout << "    Bot (A,B,C): " << A << "," << B << "," << C 
                              << " layers: " << hLayer[A] << "," << hLayer[B] << "," << hLayer[C] << std::endl;
                    std::cout << "    Tet volumes: s1=" << s1/6.0 << " s2=" << s2/6.0 << " s3=" << s3/6.0 << std::endl;
                    
                    // Print node positions
                    std::cout << "    Positions:" << std::endl;
                    std::cout << "      a(" << a << "): " << X[a] << "," << Y[a] << "," << Z[a] << std::endl;
                    std::cout << "      A(" << A << "): " << X[A] << "," << Y[A] << "," << Z[A] << std::endl;
                }
            }
            
            std::cout << "=== END DEBUG ===" << std::endl << std::endl;
            debug_call_count++;
        }
    // Compute forces and energy due to volume springs. //(nav - commenting these out for now for flat surface 5/29/24) Nav had uncommented but she's bringing the comment back because testing out Active shape mesh 02/27/25
      	ComputeVolumeSprings(
      		coordInfoVecs,
          linearSpringInfoVecs,
          capsidInfoVecs,
          generalParams,
          //auxVecs,
          prismInfoVecs);
    
    // Now print forces along the radial line
   // PrintForce();
};


void System::verifyForceDecomposition() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "      FORCE DECOMPOSITION VERIFICATION" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int N = generalParams.maxNodeCount;
    
    // ========================================================================
    // Step 1: Compute TOTAL forces (linear + volume)
    // ========================================================================
    
    // Zero out force arrays
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
    
    // Compute all forces
    Solve_Forces();  // This should compute both linear and volume forces
    
    // Copy total forces to host
    thrust::host_vector<double> h_totalFx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> h_totalFy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> h_totalFz = coordInfoVecs.nodeForceZ;
    
    // Compute max total force
    double maxTotalF = 0.0;
    for (int i = 0; i < N; ++i) {
        double F = sqrt(h_totalFx[i]*h_totalFx[i] + h_totalFy[i]*h_totalFy[i] + h_totalFz[i]*h_totalFz[i]);
        maxTotalF = std::max(maxTotalF, F);
    }
    
    std::cout << "\n--- Total Forces (from Solve_Forces) ---" << std::endl;
    std::cout << "  Max |F_total| = " << maxTotalF << std::endl;
    
    // ========================================================================
    // Step 2: Compute LINEAR forces only
    // ========================================================================
    
    // Zero out force arrays
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
    
    // Compute linear spring forces only
    ComputeLinearSprings(
        generalParams,
        coordInfoVecs,
        linearSpringInfoVecs);
    
    // Copy linear forces to host
    thrust::host_vector<double> h_linearFx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> h_linearFy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> h_linearFz = coordInfoVecs.nodeForceZ;
    
    // Compute max linear force
    double maxLinearF = 0.0;
    double sumLinearF = 0.0;
    for (int i = 0; i < N; ++i) {
        double F = sqrt(h_linearFx[i]*h_linearFx[i] + h_linearFy[i]*h_linearFy[i] + h_linearFz[i]*h_linearFz[i]);
        maxLinearF = std::max(maxLinearF, F);
        sumLinearF += F;
    }
    
    std::cout << "\n--- Linear Spring Forces ---" << std::endl;
    std::cout << "  Max |F_linear| = " << maxLinearF << std::endl;
    std::cout << "  Avg |F_linear| = " << sumLinearF / N << std::endl;
    
    // ========================================================================
    // Step 3: Compute VOLUME forces only
    // ========================================================================
    
    // Zero out force arrays
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
    
    // Compute volume forces only
    ComputeVolumeSprings(
      		coordInfoVecs,
          linearSpringInfoVecs,
          capsidInfoVecs,
          generalParams,
          //auxVecs,
          prismInfoVecs);
    
    // Copy volume forces to host
    thrust::host_vector<double> h_volumeFx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> h_volumeFy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> h_volumeFz = coordInfoVecs.nodeForceZ;
    
    // Compute max volume force
    double maxVolumeF = 0.0;
    double sumVolumeF = 0.0;
    for (int i = 0; i < N; ++i) {
        double F = sqrt(h_volumeFx[i]*h_volumeFx[i] + h_volumeFy[i]*h_volumeFy[i] + h_volumeFz[i]*h_volumeFz[i]);
        maxVolumeF = std::max(maxVolumeF, F);
        sumVolumeF += F;
    }
    
    std::cout << "\n--- Volume Forces ---" << std::endl;
    std::cout << "  Max |F_volume| = " << maxVolumeF << std::endl;
    std::cout << "  Avg |F_volume| = " << sumVolumeF / N << std::endl;
    
    // ========================================================================
    // Step 4: Verify F_total = F_linear + F_volume
    // ========================================================================
    
    std::cout << "\n--- Verification: F_total == F_linear + F_volume ---" << std::endl;
    
    double maxError = 0.0;
    double sumError = 0.0;
    int errorCount = 0;
    
    for (int i = 0; i < N; ++i) {
        // Compute expected total = linear + volume
        double expectedFx = h_linearFx[i] + h_volumeFx[i];
        double expectedFy = h_linearFy[i] + h_volumeFy[i];
        double expectedFz = h_linearFz[i] + h_volumeFz[i];
        
        // Compute error
        double errX = fabs(h_totalFx[i] - expectedFx);
        double errY = fabs(h_totalFy[i] - expectedFy);
        double errZ = fabs(h_totalFz[i] - expectedFz);
        double err = sqrt(errX*errX + errY*errY + errZ*errZ);
        
        maxError = std::max(maxError, err);
        sumError += err;
        
        if (err > 1e-10) {
            errorCount++;
            if (errorCount <= 5) {  // Print first 5 errors
                std::cout << "  Node " << i << ": error = " << err << std::endl;
                std::cout << "    Total:    (" << h_totalFx[i] << ", " << h_totalFy[i] << ", " << h_totalFz[i] << ")" << std::endl;
                std::cout << "    Linear:   (" << h_linearFx[i] << ", " << h_linearFy[i] << ", " << h_linearFz[i] << ")" << std::endl;
                std::cout << "    Volume:   (" << h_volumeFx[i] << ", " << h_volumeFy[i] << ", " << h_volumeFz[i] << ")" << std::endl;
                std::cout << "    Expected: (" << expectedFx << ", " << expectedFy << ", " << expectedFz << ")" << std::endl;
            }
        }
    }
    
    std::cout << "\n  Max error:   " << maxError << std::endl;
    std::cout << "  Avg error:   " << sumError / N << std::endl;
    std::cout << "  Nodes with error > 1e-10: " << errorCount << std::endl;
    
    if (maxError < 1e-10) {
        std::cout << "\n  ✓ VERIFIED: F_total = F_linear + F_volume" << std::endl;
    } else {
        std::cout << "\n  ✗ MISMATCH: Forces don't add up correctly!" << std::endl;
        std::cout << "    Check if Solve_Forces() includes other force terms." << std::endl;
    }
    
    // ========================================================================
    // Step 5: Print sample forces for first few nodes
    // ========================================================================
    
    std::cout << "\n--- Sample Forces (first 5 nodes) ---" << std::endl;
    for (int i = 0; i < std::min(5, N); ++i) {
        double totalMag = sqrt(h_totalFx[i]*h_totalFx[i] + h_totalFy[i]*h_totalFy[i] + h_totalFz[i]*h_totalFz[i]);
        double linearMag = sqrt(h_linearFx[i]*h_linearFx[i] + h_linearFy[i]*h_linearFy[i] + h_linearFz[i]*h_linearFz[i]);
        double volumeMag = sqrt(h_volumeFx[i]*h_volumeFx[i] + h_volumeFy[i]*h_volumeFy[i] + h_volumeFz[i]*h_volumeFz[i]);
        
        std::cout << "  Node " << i << ":" << std::endl;
        std::cout << "    |F_total|  = " << totalMag << std::endl;
        std::cout << "    |F_linear| = " << linearMag << std::endl;
        std::cout << "    |F_volume| = " << volumeMag << std::endl;
    }
    
    std::cout << std::string(60, '=') << std::endl << std::endl;
    
    // Reset forces to total for continued simulation
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
    Solve_Forces();
}


// ============================================================================
// OPTION 2: Simpler inline diagnostic (add inside relaxation loop)
//
// Add this code after Solve_Forces() to print force breakdown each iteration
// ============================================================================

/*
// After Solve_Forces():
{
    // Get current total forces
    thrust::host_vector<double> hfx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> hfy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> hfz = coordInfoVecs.nodeForceZ;
    
    double maxF = 0.0;
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        double F = sqrt(hfx[i]*hfx[i] + hfy[i]*hfy[i] + hfz[i]*hfz[i]);
        maxF = std::max(maxF, F);
    }
    
    // Get energies (these are computed in Solve_Forces or can be computed separately)
    double linearE = ComputeLinearSpringEnergy(...);  // If you have this function
    double volumeE = ComputeVolumeSpringEnergy(...);  // If you have this function
    
    std::cout << "Max|F|=" << maxF 
              << " LinearE=" << linearE 
              << " VolumeE=" << volumeE << std::endl;
}
*/


// ============================================================================
// FIXED solveSystem() function for System.cu
// 
// This version works with StrainTensor_FIXED.cu which:
// 1. Uses Gram-Schmidt orthonormalization for averaged basis vectors
// 2. Applies λ_hh = 1/λ_iso² to vertical edges for volume conservation
// 3. Properly handles both horizontal and vertical edges
//
// Replace your existing solveSystem() in System.cu with this version.
// ============================================================================
__global__ void k_interpolateRestLength(
    int E,
    const double* L_init,
    const double* L_final,
    double* L_rest,
    double t)  // t in [0, 1]
{
    int eid = blockIdx.x * blockDim.x + threadIdx.x;
    if (eid >= E) return;
    L_rest[eid] = (1.0 - t) * L_init[eid] + t * L_final[eid];
}

// ============================================================================
// solveSystem() - Main simulation driver
// 
// Following Python reference logic:
// 1. Compute target rest lengths for ALL stages ONCE at the beginning
// 2. Stages are quasi-static relaxations that interpolate between targets
// ============================================================================

void System::solveSystem()
{
    // ========================================================================
    //                          INITIALIZATION
    // ========================================================================
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "      DROSOPHILA WING DISC EVERSION SIMULATION" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Invert layer flags if geometry is inverted
    // (i.e., if higher layer numbers have lower z-coordinates)
    {
        thrust::host_vector<int> hLayer = generalParams.nodes_in_upperhem;
        int maxFlag = *std::max_element(hLayer.begin(), hLayer.end());
        
        // Check if layer 4 has lower z than expected
        double avg_z_top = 0, avg_z_bot = 0;
        int cnt_top = 0, cnt_bot = 0;
        for (int i = 0; i < hLayer.size(); ++i) {
            if (hLayer[i] == maxFlag) {
                avg_z_top += coordInfoVecs.nodeLocZ[i];
                cnt_top++;
            } else if (hLayer[i] == 0) {
                avg_z_bot += coordInfoVecs.nodeLocZ[i];
                cnt_bot++;
            }
        }
        avg_z_top /= cnt_top;
        avg_z_bot /= cnt_bot;
        
        if (avg_z_top < avg_z_bot) {
            std::cout << "WARNING: Layer flags appear inverted. Fixing..." << std::endl;
            for (int i = 0; i < hLayer.size(); ++i) {
                hLayer[i] = maxFlag - hLayer[i];
            }
            generalParams.nodes_in_upperhem = hLayer;
        }
    }
    // Build prisms for volume computation
    BuildPrismsFromLayerFlags(generalParams, coordInfoVecs, prismInfoVecs);
    
    // Setting boundary nodes
    
     generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);

    std::cout << "boundaries in upperhem = " << generalParams.boundaries_in_upperhem.size() << std::endl;

    std::vector<int> boundary_edge_list;
    std::vector<int> boundary_node_list;

    std::cout << "edges2Triangles_1 = " << coordInfoVecs.edges2Triangles_1.size() << std::endl;
    std::cout << "edges2Triangles_2 = " << coordInfoVecs.edges2Triangles_2.size() << std::endl;

   
    // Initialize all nodes as FREE
    thrust::fill(coordInfoVecs.isNodeFixed.begin(), coordInfoVecs.isNodeFixed.end(), false);
    
    // Mark ONLY boundary nodes as fixed
    for (int i = 0; i < coordInfoVecs.num_edges; i++) {
        int T1 = static_cast<int>(coordInfoVecs.edges2Triangles_1[i]);
        int T2 = static_cast<int>(coordInfoVecs.edges2Triangles_2[i]);
    
        if (T1 < 0 || T2 < 0 || T1 >= (INT_MAX - 1000) || T2 >= (INT_MAX - 1000))
            continue;
        
        if (T1 == T2) {
            // Boundary edge - fix these nodes
            generalParams.boundaries_in_lowerhem[i] = 1;
            boundary_edge_list.push_back(i);
            
            int bdry_node1 = static_cast<int>(coordInfoVecs.edges2Nodes_1[i]);
            int bdry_node2 = static_cast<int>(coordInfoVecs.edges2Nodes_2[i]);
            boundary_node_list.push_back(bdry_node1);
            boundary_node_list.push_back(bdry_node2);
            
            coordInfoVecs.isNodeFixed[bdry_node1] = true;
            coordInfoVecs.isNodeFixed[bdry_node2] = true;
        } else {
            // Interior edge - just mark it, DON'T touch isNodeFixed
            generalParams.boundaries_in_upperhem[i] = -1;
        }
    }
    
    // =============================================================================
// DIAGNOSTIC: Check edges2Triangles values
// =============================================================================
// 
// Add this diagnostic code RIGHT AFTER the boundary detection loop (after line 1179)
// to understand why all nodes are being marked as fixed.
//
// The likely cause: edges2Triangles_1[i] == edges2Triangles_2[i] for ALL edges,
// either because:
// 1. The data wasn't loaded correctly from the XML
// 2. Both arrays have the same default value (0 or -1)
// 3. The edge-to-triangle mapping is incorrect
// =============================================================================

// ADD THIS CODE after line 1179 (after the boundary detection loop closes)

    // ========== DIAGNOSTIC: Check boundary detection results ==========
    {
        std::cout << "\n=== BOUNDARY DETECTION DIAGNOSTIC ===" << std::endl;
        
        // Count how many edges have T1 == T2
        int t1_eq_t2_count = 0;
        int t1_neq_t2_count = 0;
        int invalid_count = 0;
        
        // Sample some edge values
        std::cout << "Sample edges2Triangles values:" << std::endl;
        for (int i = 0; i < std::min(20, (int)coordInfoVecs.num_edges); i++) {
            int T1 = static_cast<int>(coordInfoVecs.edges2Triangles_1[i]);
            int T2 = static_cast<int>(coordInfoVecs.edges2Triangles_2[i]);
            std::cout << "  Edge " << i << ": T1=" << T1 << ", T2=" << T2;
            if (T1 == T2) std::cout << " [BOUNDARY]";
            std::cout << std::endl;
        }
        
        // Count all edges
        for (int i = 0; i < coordInfoVecs.num_edges; i++) {
            int T1 = static_cast<int>(coordInfoVecs.edges2Triangles_1[i]);
            int T2 = static_cast<int>(coordInfoVecs.edges2Triangles_2[i]);
            
            if (T1 < 0 || T2 < 0 || T1 >= (INT_MAX - 1000) || T2 >= (INT_MAX - 1000)) {
                invalid_count++;
            } else if (T1 == T2) {
                t1_eq_t2_count++;
            } else {
                t1_neq_t2_count++;
            }
        }
        
        std::cout << "\nEdge classification:" << std::endl;
        std::cout << "  Total edges: " << coordInfoVecs.num_edges << std::endl;
        std::cout << "  Boundary edges (T1 == T2): " << t1_eq_t2_count << std::endl;
        std::cout << "  Interior edges (T1 != T2): " << t1_neq_t2_count << std::endl;
        std::cout << "  Invalid edges (skipped): " << invalid_count << std::endl;
        
        // Count fixed nodes
        thrust::host_vector<bool> h_fixed = coordInfoVecs.isNodeFixed;
        int fixed_count = thrust::count(h_fixed.begin(), h_fixed.end(), true);
        std::cout << "\nNode status after boundary detection:" << std::endl;
        std::cout << "  Fixed nodes: " << fixed_count << "/" << generalParams.maxNodeCount << std::endl;
        std::cout << "  Free nodes: " << (generalParams.maxNodeCount - fixed_count) << std::endl;
        
        if (t1_eq_t2_count == coordInfoVecs.num_edges) {
            std::cout << "\n*** WARNING: ALL edges have T1 == T2! ***" << std::endl;
            std::cout << "    This means edges2Triangles data is likely incorrect." << std::endl;
            std::cout << "    Check how edges2Triangles_1 and edges2Triangles_2 are loaded." << std::endl;
        }
        
        if (t1_neq_t2_count == 0 && invalid_count == coordInfoVecs.num_edges) {
            std::cout << "\n*** WARNING: ALL edges are invalid! ***" << std::endl;
            std::cout << "    edges2Triangles arrays may contain invalid values." << std::endl;
        }
        
        std::cout << "=== END DIAGNOSTIC ===" << std::endl << std::endl;
    }
    // ========== END DIAGNOSTIC ==========


// =============================================================================
// ALTERNATIVE: If edges2Triangles data is wrong, DON'T fix any nodes
// =============================================================================
// 
// If you want to test without any fixed nodes (let the whole mesh move freely),
// replace the entire boundary detection section with just:

/*
    // Temporarily disable boundary fixing for testing
    thrust::fill(coordInfoVecs.isNodeFixed.begin(), coordInfoVecs.isNodeFixed.end(), false);
    std::cout << "TEST MODE: All nodes set to FREE (no boundary fixing)" << std::endl;
*/

// This will let the mesh deform without any fixed boundaries.
// The mesh might drift in space, but you'll see if the deformation code works.
    // Use a small timestep for stability
    generalParams.dt = 0.01;

//    // Compute mesh center
//    double cx = 0.0, cy = 0.0, cz = 0.0;
//    int count = 0;
//    for (int i = 0; i < generalParams.maxNodeCount; i++) {
//        cx += coordInfoVecs.nodeLocX[i];
//        cy += coordInfoVecs.nodeLocY[i];
//        cz += coordInfoVecs.nodeLocZ[i];
//        count++;
//    }
//    if (count > 0) {
//        generalParams.centerX = cx / count;
//        generalParams.centerY = cy / count;
//        generalParams.centerZ = cz / count;
//    }
//    std::cout << "Mesh center: (" << generalParams.centerX << ", " 
//              << generalParams.centerY << ", " << generalParams.centerZ << ")" << std::endl;
              
    // ============================================================================
    // Computing Disc Radius for a Dome Mesh
    // 
    // The dome is a spherical cap. We need the PLANAR radius (in XY plane),
    // not the spherical radius.
    //
    //                    * <- apex of dome
    //                   /|\
    //                  / | \
    //                 /  |  \    <- spherical surface
    //                /   |   \
    //               /    |    \
    //              *-----|-----*  <- base circle (disc radius = this distance)
    //                    |
    //              disc_radius
    //
    // Method: Find the maximum XY distance from the mesh center (in XY only)
    // ============================================================================
    
    // PASTE THIS CODE in System.cu where you compute the mesh center
    // (around line 959-973 in your current code)
    
    // ============================================================================
    // Compute mesh center and disc radius
    // ============================================================================
    
    // First, compute the centroid of all nodes
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        cx += coordInfoVecs.nodeLocX[i];
        cy += coordInfoVecs.nodeLocY[i];
        cz += coordInfoVecs.nodeLocZ[i];
    }
    cx /= generalParams.maxNodeCount;
    cy /= generalParams.maxNodeCount;
    cz /= generalParams.maxNodeCount;
    
    generalParams.centerX = cx;
    generalParams.centerY = cy;
    generalParams.centerZ = cz;
    
    std::cout << "Mesh centroid: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
    
    // ============================================================================
    // Compute disc_radius = maximum PLANAR (XY) distance from center
    // This is the radius of the dome's circular base
    // ============================================================================
    
//    double max_planar_r = 0.0;
//    double max_spherical_r = 0.0;
//    
//    for (int i = 0; i < generalParams.maxNodeCount; i++) {
//        double dx = coordInfoVecs.nodeLocX[i] - cx;
//        double dy = coordInfoVecs.nodeLocY[i] - cy;
//        double dz = coordInfoVecs.nodeLocZ[i] - cz;
//        
//        // Planar distance (XY only) - this is what we want for disc_radius
//        double planar_r = sqrt(dx*dx + dy*dy);
//        max_planar_r = std::max(max_planar_r, planar_r);
//        
//        // Spherical distance (for comparison)
//        double spherical_r = sqrt(dx*dx + dy*dy + dz*dz);
//        max_spherical_r = std::max(max_spherical_r, spherical_r);
//    }
//    
//    // Set the disc radius to the maximum planar distance
//    generalParams.disc_radius = max_planar_r;
//    
//    std::cout << "Disc radius (planar XY): " << generalParams.disc_radius << std::endl;
//    std::cout << "Sphere radius (3D):      " << max_spherical_r << std::endl;
//    
    
    // ============================================================================
    // ALTERNATIVE: Compute from boundary nodes only (more accurate for dome)
    // 
    // If your mesh has boundary nodes marked, you can compute disc_radius
    // from those specifically, which might be more accurate.
    // ============================================================================
    
    /*
    double boundary_max_r = 0.0;
    int boundary_count = 0;
    
    for (int i = 0; i < coordInfoVecs.num_edges; i++) {
        // Check if this is a boundary edge (T1 == T2 in your code)
        int T1 = coordInfoVecs.edges2Triangles_1[i];
        int T2 = coordInfoVecs.edges2Triangles_2[i];
        
        if (T1 == T2) {
            // This is a boundary edge - get its nodes
            int n1 = coordInfoVecs.edges2Nodes_1[i];
            int n2 = coordInfoVecs.edges2Nodes_2[i];
            
            // Compute planar distance for both nodes
            double dx1 = coordInfoVecs.nodeLocX[n1] - cx;
            double dy1 = coordInfoVecs.nodeLocY[n1] - cy;
            double r1 = sqrt(dx1*dx1 + dy1*dy1);
            
            double dx2 = coordInfoVecs.nodeLocX[n2] - cx;
            double dy2 = coordInfoVecs.nodeLocY[n2] - cy;
            double r2 = sqrt(dx2*dx2 + dy2*dy2);
            
            boundary_max_r = std::max(boundary_max_r, std::max(r1, r2));
            boundary_count++;
        }
    }
    
    if (boundary_count > 0) {
        std::cout << "Disc radius from boundary nodes: " << boundary_max_r << std::endl;
        generalParams.disc_radius = boundary_max_r;
    }
    
    */
    
    // ============================================================================
    // ALTERNATIVE 2: Use nodes in a specific layer (e.g., apical layer only)
    //
    // Since your dome has multiple layers, you might want to compute
    // disc_radius from just one layer (e.g., the apical/top layer)
    // ============================================================================
    
    
    thrust::host_vector<int> h_layer = generalParams.nodes_in_upperhem;
    int max_layer = *std::max_element(h_layer.begin(), h_layer.end());
    
    double apical_max_r = 0.0;
    int apical_count = 0;
    
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        if (h_layer[i] == max_layer) {  // Apical layer
            double dx = coordInfoVecs.nodeLocX[i] - cx;
            double dy = coordInfoVecs.nodeLocY[i] - cy;
            double planar_r = sqrt(dx*dx + dy*dy);
            apical_max_r = std::max(apical_max_r, planar_r);
            apical_count++;
        }
    }
    
    std::cout << "Disc radius from apical layer (" << apical_count << " nodes): " 
              << apical_max_r << std::endl;
    generalParams.disc_radius = apical_max_r;
    

    // Compute initial volume
    ComputeVolume(generalParams, coordInfoVecs, linearSpringInfoVecs, 
                  ljInfoVecs, prismInfoVecs);
    generalParams.eq_total_volume = generalParams.current_total_volume;
    std::cout << "Initial volume: " << generalParams.eq_total_volume << std::endl;

    // ========================================================================
    //                       INITIAL RELAXATION
    // ========================================================================
    
    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << "PHASE 1: Initial Relaxation (no strain)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // Set rest lengths to initial lengths
    thrust::copy(linearSpringInfoVecs.edge_initial_length.begin(),
                 linearSpringInfoVecs.edge_initial_length.end(),
                 linearSpringInfoVecs.edge_rest_length.begin());

    generalParams.tol = 1e-4;
    int initial_relax_iters = relaxUntilConvergedWithParams(
    *this,
    0.05,      // force_tolerance - stop when max|F| < 0.1
    1e-3,     // displacement_tolerance (secondary)
    10000,   // max_iterations
    100);     // print every 1000 iterations
//    storage->print_VTK_File();
//    
//    // After initial_relax_iters
//    coordInfoVecs.nodeLocX[0] += 1.0;  // Perturb node 0
//    cudaDeviceSynchronize();
//    storage->print_VTK_File();
//    int initial_relax_iters_1 = relaxUntilConverged(*this); 
//    storage->print_VTK_File();
    double E_initial = linearSpringInfoVecs.linear_spring_energy;
    std::cout << "Initial relaxation before node 0 pulling : " << initial_relax_iters << " iterations" << std::endl;
   // std::cout << "Initial relaxation after node 0 pulling : " << initial_relax_iters_1 << " iterations" << std::endl;
    std::cout << "Initial energy: " << E_initial << std::endl;

    storage->print_VTK_File();


    // ========================================================================
    //            COMPUTE BASIS VECTORS WITH DV SEPARATION
    // ========================================================================
    
    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << "PHASE 2: Computing Basis Vectors with DV Separation" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // DV boundary parameters
    double theta_DV = 0.1931;  // DV boundary angle (radians) - ~11 degrees
    
    // Auto-detect mesh radius from coordinates
    double R = 0.0;
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        double x = coordInfoVecs.nodeLocX[i];
        double y = coordInfoVecs.nodeLocY[i];
        double z = coordInfoVecs.nodeLocZ[i];
        double r = sqrt(x*x + y*y + z*z);
        R = std::max(R, r);
    }
    if (R < 1e-10) R = 1.0;
    std::cout << "Detected mesh radius R = " << R << std::endl;
    
    
    // Compute basis vectors with DV separation
    // This will:
    // 1. Classify nodes into DV, dorsal, and ventral regions
    // 2. Compute region-specific origins (OV-OD line, OD, OV)
    // 3. Compute radial basis vectors from each origin
    // 4. Set up pathlengths for lambda interpolation
    StrainTensorGPU::computeBasisVectorsWithDVSeparation(
        generalParams, coordInfoVecs, theta_DV, R);
      // Build vertex-level lambda
        LambdaField lambda;
        
        // ====================================================================
        // FIX: Copy basis vectors from CoordInfoVecs to LambdaField
        // The basis vectors were computed by computeBasisVectorsWithDVSeparation()
        // and stored in coordInfoVecs, but buildVertexLambda() reads from lambda.
        // ====================================================================
        {
            int N = generalParams.maxNodeCount;
            
            // Resize the lambda field
            lambda.resize(N);
            
            // Copy basis vectors from CoordInfoVecs to LambdaField
            thrust::host_vector<double> h_eR_x = coordInfoVecs.e_R_x;
            thrust::host_vector<double> h_eR_y = coordInfoVecs.e_R_y;
            thrust::host_vector<double> h_eR_z = coordInfoVecs.e_R_z;
            
            thrust::host_vector<double> h_ePhi_x = coordInfoVecs.e_phi_x;
            thrust::host_vector<double> h_ePhi_y = coordInfoVecs.e_phi_y;
            thrust::host_vector<double> h_ePhi_z = coordInfoVecs.e_phi_z;
            
            thrust::host_vector<double> h_eH_x = coordInfoVecs.e_h_x;
            thrust::host_vector<double> h_eH_y = coordInfoVecs.e_h_y;
            thrust::host_vector<double> h_eH_z = coordInfoVecs.e_h_z;
            
            // Create host vectors for CVec3 basis vectors
            thrust::host_vector<CVec3> h_e_R(N);
            thrust::host_vector<CVec3> h_e_phi(N);
            thrust::host_vector<CVec3> h_e_h(N);
            
            for (int i = 0; i < N; ++i) {
                h_e_R[i] = thrust::make_tuple(h_eR_x[i], h_eR_y[i], h_eR_z[i]);
                h_e_phi[i] = thrust::make_tuple(h_ePhi_x[i], h_ePhi_y[i], h_ePhi_z[i]);
                h_e_h[i] = thrust::make_tuple(h_eH_x[i], h_eH_y[i], h_eH_z[i]);
            }
            
            // Copy to device (LambdaField)
            lambda.e_R = h_e_R;
            lambda.e_phi = h_e_phi;
            lambda.e_h = h_e_h;
            
            std::cout << "Basis vectors copied to LambdaField. N = " << N << std::endl;
        }
        
        // FIX: Ensure lambda_aniso_edge_outDV = 1.0 (your XML might have 0.5)
        generalParams.lambda_aniso_edge_outDV = 1.0;
        double frac = 1.0;   // full-field application per stage
        
        // Strain field (lambda) values in inDV and outDV regions at different stages of eversion. 
        generalParams.lambda_iso_center_outDV = 1.0;
        generalParams.lambda_iso_edge_outDV = 1.0;
        generalParams.lambda_aniso_center_outDV = 1.0;
        generalParams.lambda_aniso_edge_outDV = 1.0;
        
        generalParams.lambda_iso_center_inDV = 1.5; //    wl3-0hAPF (-0.09848994) | wl3-2hAPF (-0.11692544) | wl3-4hAPF (-0.06151876)
        generalParams.lambda_iso_edge_inDV = 1.0;//       wl3-0hAPF ( 1.18401136) | wl3-2hAPF ( 1.21007540) | wl3-4hAPF ( 1.47472744)
        generalParams.lambda_aniso_center_inDV = 1.0; //  wl3-0hAPF (-0.12904887) | wl3-2hAPF (-0.21271504) | wl3-4hAPF (-0.30567972)
        generalParams.lambda_aniso_edge_inDV = 1.5; //    wl3-0hAPF ( 1.03128453) | wl3-2hAPF ( 1.24178074) | wl3-4hAPF ( 1.29370391)
        
        
        // NOW build the lambda field (it will use the basis vectors we just copied)
        StrainTensorGPU::buildVertexLambda(generalParams, coordInfoVecs, lambda, frac);
        
        // ============================================================================
        // DIAGNOSTIC: Check if basis vectors are orthonormal and tensor is identity
        //
        // Add this code AFTER buildVertexLambda() and BEFORE updateEdgeRestLengths()
        // to verify the strain tensor is correct.
        // ============================================================================
        
        {
            std::cout << "\n" << std::string(60, '=') << std::endl;
            std::cout << "      BASIS VECTOR AND TENSOR DIAGNOSTICS" << std::endl;
            std::cout << std::string(60, '=') << std::endl;
            
            int N = generalParams.maxNodeCount;
            
            // Copy basis vectors from LambdaField to host
            thrust::host_vector<CVec3> h_e_R = lambda.e_R;
            thrust::host_vector<CVec3> h_e_phi = lambda.e_phi;
            thrust::host_vector<CVec3> h_e_h = lambda.e_h;
            thrust::host_vector<Mat_3x3> h_lam_alpha = lambda.lam_alpha;
            
            // Check first few nodes
            int num_to_check = std::min(5, N);
            
            std::cout << "\n--- Checking first " << num_to_check << " nodes ---" << std::endl;
            
            double max_non_orthonormal = 0.0;
            double max_tensor_error = 0.0;
            
            for (int i = 0; i < num_to_check; ++i) {
                double eR_x = thrust::get<0>(h_e_R[i]);
                double eR_y = thrust::get<1>(h_e_R[i]);
                double eR_z = thrust::get<2>(h_e_R[i]);
                
                double ePhi_x = thrust::get<0>(h_e_phi[i]);
                double ePhi_y = thrust::get<1>(h_e_phi[i]);
                double ePhi_z = thrust::get<2>(h_e_phi[i]);
                
                double eH_x = thrust::get<0>(h_e_h[i]);
                double eH_y = thrust::get<1>(h_e_h[i]);
                double eH_z = thrust::get<2>(h_e_h[i]);
                
                // Check magnitudes (should be 1)
                double mag_R = sqrt(eR_x*eR_x + eR_y*eR_y + eR_z*eR_z);
                double mag_phi = sqrt(ePhi_x*ePhi_x + ePhi_y*ePhi_y + ePhi_z*ePhi_z);
                double mag_h = sqrt(eH_x*eH_x + eH_y*eH_y + eH_z*eH_z);
                
                // Check orthogonality (dot products should be 0)
                double dot_R_phi = eR_x*ePhi_x + eR_y*ePhi_y + eR_z*ePhi_z;
                double dot_R_h = eR_x*eH_x + eR_y*eH_y + eR_z*eH_z;
                double dot_phi_h = ePhi_x*eH_x + ePhi_y*eH_y + ePhi_z*eH_z;
                
                std::cout << "\nNode " << i << ":" << std::endl;
                std::cout << "  e_R   = (" << eR_x << ", " << eR_y << ", " << eR_z << "), |e_R| = " << mag_R << std::endl;
                std::cout << "  e_phi = (" << ePhi_x << ", " << ePhi_y << ", " << ePhi_z << "), |e_phi| = " << mag_phi << std::endl;
                std::cout << "  e_h   = (" << eH_x << ", " << eH_y << ", " << eH_z << "), |e_h| = " << mag_h << std::endl;
                std::cout << "  Orthogonality: e_R·e_phi=" << dot_R_phi << ", e_R·e_h=" << dot_R_h << ", e_phi·e_h=" << dot_phi_h << std::endl;
                
                // Check tensor (should be identity when lambda=1)
                Mat_3x3& L = h_lam_alpha[i];
                CVec3 row0 = thrust::get<0>(L);
                CVec3 row1 = thrust::get<1>(L);
                CVec3 row2 = thrust::get<2>(L);
                
                std::cout << "  Tensor L = [" << std::endl;
                std::cout << "    [" << thrust::get<0>(row0) << ", " << thrust::get<1>(row0) << ", " << thrust::get<2>(row0) << "]" << std::endl;
                std::cout << "    [" << thrust::get<0>(row1) << ", " << thrust::get<1>(row1) << ", " << thrust::get<2>(row1) << "]" << std::endl;
                std::cout << "    [" << thrust::get<0>(row2) << ", " << thrust::get<1>(row2) << ", " << thrust::get<2>(row2) << "]" << std::endl;
                std::cout << "  ]" << std::endl;
                
                // Error from identity
                double err_00 = fabs(thrust::get<0>(row0) - 1.0);
                double err_11 = fabs(thrust::get<1>(row1) - 1.0);
                double err_22 = fabs(thrust::get<2>(row2) - 1.0);
                double err_01 = fabs(thrust::get<1>(row0));
                double err_02 = fabs(thrust::get<2>(row0));
                double err_10 = fabs(thrust::get<0>(row1));
                double err_12 = fabs(thrust::get<2>(row1));
                double err_20 = fabs(thrust::get<0>(row2));
                double err_21 = fabs(thrust::get<1>(row2));
                
                double max_diag_err = std::max({err_00, err_11, err_22});
                double max_offdiag_err = std::max({err_01, err_02, err_10, err_12, err_20, err_21});
                
                std::cout << "  Error from Identity: diag_err=" << max_diag_err << ", offdiag_err=" << max_offdiag_err << std::endl;
                
                max_tensor_error = std::max(max_tensor_error, std::max(max_diag_err, max_offdiag_err));
                max_non_orthonormal = std::max({max_non_orthonormal, fabs(mag_R-1), fabs(mag_phi-1), fabs(mag_h-1),
                                                fabs(dot_R_phi), fabs(dot_R_h), fabs(dot_phi_h)});
            }
            
            // Check ALL nodes for overall statistics
            std::cout << "\n--- Checking ALL " << N << " nodes ---" << std::endl;
            
            double sum_mag_R = 0, sum_mag_phi = 0, sum_mag_h = 0;
            double min_mag_R = 1e10, min_mag_phi = 1e10, min_mag_h = 1e10;
            double max_mag_R = 0, max_mag_phi = 0, max_mag_h = 0;
            int count_zero_basis = 0;
            
            for (int i = 0; i < N; ++i) {
                double eR_x = thrust::get<0>(h_e_R[i]);
                double eR_y = thrust::get<1>(h_e_R[i]);
                double eR_z = thrust::get<2>(h_e_R[i]);
                
                double ePhi_x = thrust::get<0>(h_e_phi[i]);
                double ePhi_y = thrust::get<1>(h_e_phi[i]);
                double ePhi_z = thrust::get<2>(h_e_phi[i]);
                
                double eH_x = thrust::get<0>(h_e_h[i]);
                double eH_y = thrust::get<1>(h_e_h[i]);
                double eH_z = thrust::get<2>(h_e_h[i]);
                
                double mag_R = sqrt(eR_x*eR_x + eR_y*eR_y + eR_z*eR_z);
                double mag_phi = sqrt(ePhi_x*ePhi_x + ePhi_y*ePhi_y + ePhi_z*ePhi_z);
                double mag_h = sqrt(eH_x*eH_x + eH_y*eH_y + eH_z*eH_z);
                
                sum_mag_R += mag_R;
                sum_mag_phi += mag_phi;
                sum_mag_h += mag_h;
                
                min_mag_R = std::min(min_mag_R, mag_R);
                min_mag_phi = std::min(min_mag_phi, mag_phi);
                min_mag_h = std::min(min_mag_h, mag_h);
                
                max_mag_R = std::max(max_mag_R, mag_R);
                max_mag_phi = std::max(max_mag_phi, mag_phi);
                max_mag_h = std::max(max_mag_h, mag_h);
                
                if (mag_R < 1e-10 || mag_phi < 1e-10 || mag_h < 1e-10) {
                    count_zero_basis++;
                }
            }
            
            std::cout << "  |e_R|:   min=" << min_mag_R << ", max=" << max_mag_R << ", avg=" << sum_mag_R/N << std::endl;
            std::cout << "  |e_phi|: min=" << min_mag_phi << ", max=" << max_mag_phi << ", avg=" << sum_mag_phi/N << std::endl;
            std::cout << "  |e_h|:   min=" << min_mag_h << ", max=" << max_mag_h << ", avg=" << sum_mag_h/N << std::endl;
            std::cout << "  Nodes with zero-magnitude basis vectors: " << count_zero_basis << std::endl;
            
            if (count_zero_basis > 0) {
                std::cout << "\n*** WARNING: " << count_zero_basis << " nodes have zero basis vectors! ***" << std::endl;
                std::cout << "    This will cause incorrect strain calculation!" << std::endl;
            }
            
            if (max_tensor_error > 0.01) {
                std::cout << "\n*** WARNING: Tensor deviates from identity by " << max_tensor_error << " ***" << std::endl;
                std::cout << "    With lambda=1, tensor should be identity!" << std::endl;
            }
            
            std::cout << std::string(60, '=') << std::endl << std::endl;
        }

        // Update rest lengths (initial_length → final_length)
        //int layerflag = 0;
        //StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs, generalParams, lambda, linearSpringInfoVecs, layerflag);

//    // Copy DV classification to GeneralParams for use in buildVertexLambda
//    generalParams.nodes_in_DV.resize(hostSetInfoVecs.nodes_in_DV.size());
//    thrust::copy(hostSetInfoVecs.nodes_in_DV.begin(), 
//                 hostSetInfoVecs.nodes_in_DV.end(),
//                 generalParams.nodes_in_DV.begin());
    // Compute basis vectors ONCE using initial geometry
    //LambdaField lambda;
   // double theta_DV = 0.1931;  // DV boundary angle (radians)
    //double R = 1.0;            // Will be auto-detected from mesh
    
   // StrainTensorGPU::computeBasisVectorsAndPathlength(
     //   generalParams, coordInfoVecs, lambda, theta_DV, R);
// ============================================
    //              STRAIN TENSOR STAGES
    // ============================================

    int stages = 2;//generalParams.Tf;
    
    linearSpringInfoVecs.edge_rest_length = linearSpringInfoVecs.edge_initial_length; // copy initial lengths to the rest length vector. 

    for (int stage = 0; stage < stages; stage++)
    {
        
        
        // -- load lambda values for this stage --
//        generalParams.lambda_iso_center_outDV = generalParams.lambda_iso_center_outDV_v[stage];
//        generalParams.lambda_iso_edge_outDV   = generalParams.lambda_iso_edge_outDV_v[stage];
//        generalParams.lambda_aniso_center_outDV = generalParams.lambda_aniso_center_outDV_v[stage];
//        generalParams.lambda_aniso_edge_outDV   = generalParams.lambda_aniso_edge_outDV_v[stage];
//
//        generalParams.lambda_iso_center_inDV = generalParams.lambda_iso_center_inDV_v[stage];
//        generalParams.lambda_iso_edge_inDV   = generalParams.lambda_iso_edge_inDV_v[stage];
//        generalParams.lambda_aniso_center_inDV = generalParams.lambda_aniso_center_inDV_v[stage];
//        generalParams.lambda_aniso_edge_inDV   = generalParams.lambda_aniso_edge_inDV_v[stage];

        // Build vertex-level lambda
        //LambdaField lambda;
        //StrainTensorGPU::buildVertexLambda(generalParams, coordInfoVecs, lambda, frac);

        // Update rest lengths (initial_length → final_length)
        int layerflag = 0;   // 0 = basal layer, 
                             // 1 - N = Body layers,
                             // N+1 = Apical layer,
                             // -1 = vertical layer
        // ============================================================================
        // DIAGNOSTIC CODE: Print Lambda Values and Strain Statistics
        // 
        // Paste this code BEFORE the call to StrainTensorGPU::updateEdgeRestLengths()
        // in your System.cu solveSystem() function.
        //
        // This will print:
        //   1. Lambda parameter values being used
        //   2. Strain statistics (min, max, average) for the computed rest lengths
        // ============================================================================
        
        
        {
            std::cout << "\n" << std::string(60, '=') << std::endl;
            std::cout << "      LAMBDA VALUES AND STRAIN DIAGNOSTICS" << std::endl;
            std::cout << std::string(60, '=') << std::endl;
            
           // generalParams.lambda_aniso_edge_outDV = 1.2;
           
          
            
            // Print lambda parameter values
            std::cout << "\n--- Lambda Parameters (from GeneralParams) ---" << std::endl;
            std::cout << "Outside DV Region:" << std::endl;
            std::cout << "  lambda_iso_center_outDV   = " << generalParams.lambda_iso_center_outDV << std::endl;
            std::cout << "  lambda_iso_edge_outDV     = " << generalParams.lambda_iso_edge_outDV << std::endl;
            std::cout << "  lambda_aniso_center_outDV = " << generalParams.lambda_aniso_center_outDV << std::endl;
            std::cout << "  lambda_aniso_edge_outDV   = " << generalParams.lambda_aniso_edge_outDV << std::endl;
            
            std::cout << "Inside DV Region:" << std::endl;
            std::cout << "  lambda_iso_center_inDV    = " << generalParams.lambda_iso_center_inDV << std::endl;
            std::cout << "  lambda_iso_edge_inDV      = " << generalParams.lambda_iso_edge_inDV << std::endl;
            std::cout << "  lambda_aniso_center_inDV  = " << generalParams.lambda_aniso_center_inDV << std::endl;
            std::cout << "  lambda_aniso_edge_inDV    = " << generalParams.lambda_aniso_edge_inDV << std::endl;
            
            std::cout << "\nGeometry Parameters:" << std::endl;
            std::cout << "  disc_radius = " << generalParams.disc_radius << std::endl;
            std::cout << "  centerX     = " << generalParams.centerX << std::endl;
            std::cout << "  centerY     = " << generalParams.centerY << std::endl;
            std::cout << "  centerZ     = " << generalParams.centerZ << std::endl;
            
            // Print per-vertex lambda field statistics (if lambda field has been built)
            if (lambda.lam_rr.size() > 0) {
                thrust::host_vector<double> h_lam_rr = lambda.lam_rr;
                thrust::host_vector<double> h_lam_pp = lambda.lam_pp;
                thrust::host_vector<double> h_lam_ss = lambda.lam_ss;
                thrust::host_vector<double> h_rho = generalParams.rho;
                
                double min_rr = 1e10, max_rr = -1e10, sum_rr = 0;
                double min_pp = 1e10, max_pp = -1e10, sum_pp = 0;
                double min_ss = 1e10, max_ss = -1e10, sum_ss = 0;
                double min_rho = 1e10, max_rho = -1e10, sum_rho = 0;
                
                int N = h_lam_rr.size();
                for (int i = 0; i < N; i++) {
                    min_rr = std::min(min_rr, h_lam_rr[i]);
                    max_rr = std::max(max_rr, h_lam_rr[i]);
                    sum_rr += h_lam_rr[i];
                    
                    min_pp = std::min(min_pp, h_lam_pp[i]);
                    max_pp = std::max(max_pp, h_lam_pp[i]);
                    sum_pp += h_lam_pp[i];
                    
                    min_ss = std::min(min_ss, h_lam_ss[i]);
                    max_ss = std::max(max_ss, h_lam_ss[i]);
                    sum_ss += h_lam_ss[i];
                    
                    if (h_rho.size() > i) {
                        min_rho = std::min(min_rho, h_rho[i]);
                        max_rho = std::max(max_rho, h_rho[i]);
                        sum_rho += h_rho[i];
                    }
                }
                
                std::cout << "\n--- Per-Vertex Lambda Field Statistics (N = " << N << " vertices) ---" << std::endl;
                std::cout << "  lambda_rr (radial):       min=" << min_rr << ", max=" << max_rr << ", avg=" << sum_rr/N << std::endl;
                std::cout << "  lambda_pp (circumfer.):   min=" << min_pp << ", max=" << max_pp << ", avg=" << sum_pp/N << std::endl;
                std::cout << "  lambda_ss (thickness):    min=" << min_ss << ", max=" << max_ss << ", avg=" << sum_ss/N << std::endl;
                std::cout << "  rho (normalized radius):  min=" << min_rho << ", max=" << max_rho << ", avg=" << sum_rho/N << std::endl;
            }
            
            std::cout << std::string(60, '=') << std::endl;
        }
        
        // --- NOW CALL updateEdgeRestLengths() ---

          StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs, generalParams, lambda, linearSpringInfoVecs, layerflag);

          // --- PASTE THIS BLOCK AFTER updateEdgeRestLengths() ---
          
          {
              std::cout << "\n" << std::string(60, '=') << std::endl;
              std::cout << "      STRAIN STATISTICS (After Rest Length Update)" << std::endl;
              std::cout << std::string(60, '=') << std::endl;
              
              // Copy device vectors to host
              thrust::host_vector<double> h_init = linearSpringInfoVecs.edge_initial_length;
              thrust::host_vector<double> h_final = linearSpringInfoVecs.edge_final_length;
              thrust::host_vector<double> h_rest = linearSpringInfoVecs.edge_rest_length;
              thrust::host_vector<int> h_layer = generalParams.edges_in_upperhem;
              
              int E = coordInfoVecs.num_edges;
              
              // Compute strain = (L_final - L_init) / L_init for each edge
              double min_strain = 1e10, max_strain = -1e10, sum_strain = 0;
              double min_ratio = 1e10, max_ratio = -1e10, sum_ratio = 0;
              int count_positive_strain = 0;
              int count_negative_strain = 0;
              int count_zero_strain = 0;
              int count_horizontal = 0;
              int count_vertical = 0;
              
              // Per-layer statistics
              std::map<int, std::vector<double>> strain_by_layer;
              
              for (int e = 0; e < E; e++) {
                  double L0 = h_init[e];
                  double Lf = h_final[e];
                  int layer = h_layer[e];
                  
                  if (L0 > 1e-10) {
                      double strain = (Lf - L0) / L0;  // Engineering strain
                      double ratio = Lf / L0;           // Stretch ratio
                      
                      min_strain = std::min(min_strain, strain);
                      max_strain = std::max(max_strain, strain);
                      sum_strain += strain;
                      
                      min_ratio = std::min(min_ratio, ratio);
                      max_ratio = std::max(max_ratio, ratio);
                      sum_ratio += ratio;
                      
                      if (strain > 1e-9) count_positive_strain++;
                      else if (strain < -1e-9) count_negative_strain++;
                      else count_zero_strain++;
                      
                      if (layer == -1) count_vertical++;
                      else count_horizontal++;
                      
                      strain_by_layer[layer].push_back(strain);
                  }
              }
              
              std::cout << "\n--- Overall Strain Statistics (E = " << E << " edges) ---" << std::endl;
              std::cout << "  Engineering Strain = (L_final - L_init) / L_init" << std::endl;
              std::cout << "  Min strain:  " << min_strain << " (" << min_strain*100 << "%)" << std::endl;
              std::cout << "  Max strain:  " << max_strain << " (" << max_strain*100 << "%)" << std::endl;
              std::cout << "  Avg strain:  " << sum_strain/E << " (" << (sum_strain/E)*100 << "%)" << std::endl;
              
              std::cout << "\n--- Stretch Ratio Statistics ---" << std::endl;
              std::cout << "  Stretch Ratio = L_final / L_init" << std::endl;
              std::cout << "  Min ratio:   " << min_ratio << std::endl;
              std::cout << "  Max ratio:   " << max_ratio << std::endl;
              std::cout << "  Avg ratio:   " << sum_ratio/E << std::endl;
              
              std::cout << "\n--- Edge Classification ---" << std::endl;
              std::cout << "  Edges with positive strain (tension):    " << count_positive_strain << std::endl;
              std::cout << "  Edges with negative strain (compress.):  " << count_negative_strain << std::endl;
              std::cout << "  Edges with zero strain:                  " << count_zero_strain << std::endl;
              std::cout << "  Horizontal edges (layer >= 0):           " << count_horizontal << std::endl;
              std::cout << "  Vertical edges (layer == -1):            " << count_vertical << std::endl;
              
              std::cout << "\n--- Strain by Layer ---" << std::endl;
              for (auto& kv : strain_by_layer) {
                  int layer = kv.first;
                  std::vector<double>& strains = kv.second;
                  
                  double layer_min = 1e10, layer_max = -1e10, layer_sum = 0;
                  for (double s : strains) {
                      layer_min = std::min(layer_min, s);
                      layer_max = std::max(layer_max, s);
                      layer_sum += s;
                  }
                  
                  std::string layer_name;
                  if (layer == -1) layer_name = "Vertical";
                  else if (layer == 0) layer_name = "Basal";
                  else if (layer == 4) layer_name = "Apical";  // Adjust if your max layer differs
                  else layer_name = "Body " + std::to_string(layer);
                  
                  std::cout << "  Layer " << layer << " (" << layer_name << ", n=" << strains.size() << "): "
                            << "min=" << layer_min*100 << "%, max=" << layer_max*100 << "%, avg=" << (layer_sum/strains.size())*100 << "%" 
                            << std::endl;
              }
              
              // Check if strain is actually being applied
              if (max_strain < 1e-9 && min_strain > -1e-9) {
                  std::cout << "\n*** WARNING: No strain detected! ***" << std::endl;
                  std::cout << "    Possible causes:" << std::endl;
                  std::cout << "    1. All lambda values = 1.0 (no growth)" << std::endl;
                  std::cout << "    2. Strain tensor not being computed correctly" << std::endl;
                  std::cout << "    3. Edge layer flags preventing strain application" << std::endl;
              }
              
              std::cout << std::string(60, '=') << std::endl << std::endl;
          }

        //StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs, generalParams, lambda, linearSpringInfoVecs, layerflag);

        // Relaxation parameters
        //generalParams.tol = 1e-4;
        int Nsteps = 1000; // this should be equal to the inverse of tolerance

        // ============================================
        //     GRADIENT RELAXATION LOOP
        // ============================================
        
        for (int iter = 0; iter < Nsteps; iter++)
        {
//            // linearly increment rest lengths if doing time sweep
//            for (int e = 0; e < coordInfoVecs.num_edges; e++) {
//                double dl = (linearSpringInfoVecs.edge_final_length[e] -
//                             linearSpringInfoVecs.edge_initial_length[e]) / double(Nsteps);
//                linearSpringInfoVecs.edge_rest_length[e] += dl;
//            }
            double t = (iter + 1.0) / double(Nsteps);
            dim3 grid((coordInfoVecs.num_edges + 255) / 256);
            k_interpolateRestLength<<<grid, 256>>>(
                coordInfoVecs.num_edges,
                thrust::raw_pointer_cast(linearSpringInfoVecs.edge_initial_length.data()),
                thrust::raw_pointer_cast(linearSpringInfoVecs.edge_final_length.data()),
                thrust::raw_pointer_cast(linearSpringInfoVecs.edge_rest_length.data()),
                t);
            cudaDeviceSynchronize();
            // In System.cu:
            int k = relaxUntilConvergedWithParams(
                *this,
                0.05,      // force_tolerance - stop when max|F| < 0.1
                1e-4,     // displacement_tolerance (secondary)
                10000,   // max_iterations
                10);    // print every 1000 iterations
            
            //verifyForceDecomposition();
            
            double E = linearSpringInfoVecs.linear_spring_energy+generalParams.volume_energy;
            std::cout << "Relax iter " << iter 
                      << " Stage " << stage 
                      << " | E = " << E 
                      << " | Linear E = " << linearSpringInfoVecs.linear_spring_energy
                      << " | Volume E = " << generalParams.volume_energy
                      << " | Mov = " << generalParams.dx 
                      << " | Volume = " << generalParams.current_total_volume
                      << " | Steps = " << k << std::endl;

            if (iter % 5 == 0)
                storage->print_VTK_File();
                
//            BuildPrismsFromLayerFlags(
//              generalParams,
//              coordInfoVecs,
//              prismInfoVecs);
        }
//        std::cout //"Relax iter " << iter 
//                      << " Stage " << stage 
//                      //<< " | E = " << E 
//                      << " | Mov = " << generalParams.dx 
//                      << " | Volume = " << generalParams.current_total_volume << std::endl;

        storage->print_VTK_File();
        linearSpringInfoVecs.edge_initial_length = linearSpringInfoVecs.edge_rest_length;
    }
    std::cout << std::string(60, '=') << std::endl;
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
    
    coordInfoVecs.pathlength_scaled.resize(mem_prealloc * coordInfoVecs.pathlength_scaled.size(), 0.0);
    
    coordInfoVecs.e_R_x.resize(mem_prealloc * coordInfoVecs.e_R_x.size(), 0.0);
    coordInfoVecs.e_R_x.resize(mem_prealloc * coordInfoVecs.e_R_y.size(), 0.0);
    coordInfoVecs.e_R_x.resize(mem_prealloc * coordInfoVecs.e_R_z.size(), 0.0);
    
    coordInfoVecs.e_phi_x.resize(mem_prealloc * coordInfoVecs.e_phi_x.size(), 0.0);
    coordInfoVecs.e_phi_y.resize(mem_prealloc * coordInfoVecs.e_phi_y.size(), 0.0);
    coordInfoVecs.e_phi_z.resize(mem_prealloc * coordInfoVecs.e_phi_z.size(), 0.0);
    
    coordInfoVecs.e_h_x.resize(mem_prealloc * coordInfoVecs.e_h_x.size(), 0.0);
    coordInfoVecs.e_h_y.resize(mem_prealloc * coordInfoVecs.e_h_y.size(), 0.0);
    coordInfoVecs.e_h_z.resize(mem_prealloc * coordInfoVecs.e_h_z.size(), 0.0);
    

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

    int invalidEdges = 0;
    for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
        int n1 = coordInfoVecs.edges2Nodes_1[e];
        int n2 = coordInfoVecs.edges2Nodes_2[e];
        
        if (n1 < 0 || n1 >= generalParams.maxNodeCount ||
            n2 < 0 || n2 >= generalParams.maxNodeCount ||
            n1 == INT_MAX || n2 == INT_MAX) {
            invalidEdges++;
            coordInfoVecs.edges2Nodes_1[e] = 0;  // Fix to valid
            coordInfoVecs.edges2Nodes_2[e] = 0;
        }
    }
    std::cout << "Invalid edges found and fixed: " << invalidEdges << std::endl;

  
    // copy info to GPU
    std::cout << "Copying" << std::endl;
    //thrust::copy(hostSetInfoVecs.isNodeFixed.begin(), hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());

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
    
    prismInfoVecs.num_prisms = 0;
    prismInfoVecs.P1.clear(); prismInfoVecs.P2.clear(); prismInfoVecs.P3.clear();
    prismInfoVecs.P4.clear(); prismInfoVecs.P5.clear(); prismInfoVecs.P6.clear();

    
    // ==================== initialize lambda stage vectors ====================

//    // outDV
//    generalParams.lambda_iso_center_outDV_v.resize(3);
//    generalParams.lambda_iso_edge_outDV_v.resize(3);
//    generalParams.lambda_aniso_center_outDV_v.resize(3);
//    generalParams.lambda_aniso_edge_outDV_v.resize(3);
//    
//    // inDV
//    generalParams.lambda_iso_center_inDV_v.resize(3);
//    generalParams.lambda_iso_edge_inDV_v.resize(3);
//    generalParams.lambda_aniso_center_inDV_v.resize(3);
//    generalParams.lambda_aniso_edge_inDV_v.resize(3);
    
    // temporary: no growth (all ones)
//    thrust::fill(generalParams.lambda_iso_center_outDV_v.begin(),
//                 generalParams.lambda_iso_center_outDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_iso_edge_outDV_v.begin(),
//                 generalParams.lambda_iso_edge_outDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_aniso_center_outDV_v.begin(),
//                 generalParams.lambda_aniso_center_outDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_aniso_edge_outDV_v.begin(),
//                 generalParams.lambda_aniso_edge_outDV_v.end(), 1.0);
//    
//    thrust::fill(generalParams.lambda_iso_center_inDV_v.begin(),
//                 generalParams.lambda_iso_center_inDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_iso_edge_inDV_v.begin(),
//                 generalParams.lambda_iso_edge_inDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_aniso_center_inDV_v.begin(),
//                 generalParams.lambda_aniso_center_inDV_v.end(), 1.0);
//    thrust::fill(generalParams.lambda_aniso_edge_inDV_v.begin(),
//                 generalParams.lambda_aniso_edge_inDV_v.end(), 1.0);
    
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

    const int E = coordInfoVecs.num_edges;

    // These MUST be size E
    linearSpringInfoVecs.edge_initial_length.resize(E);
    linearSpringInfoVecs.edge_rest_length.resize(E);
    linearSpringInfoVecs.edge_final_length.resize(E);
    
    // Copy from host (host vectors must also be size E)
    linearSpringInfoVecs.edge_initial_length = hostSetInfoVecs.edge_initial_length;
    linearSpringInfoVecs.edge_rest_length    = hostSetInfoVecs.edge_initial_length;
    linearSpringInfoVecs.edge_final_length   = hostSetInfoVecs.edge_initial_length;

//    // linearSpringInfoVecs.edge_rest_length.resize(hostSetInfoVecs.edge_rest_length.size());
//    linearSpringInfoVecs.edge_final_length.resize(hostSetInfoVecs.edge_initial_length.size());
//    linearSpringInfoVecs.edge_initial_length = hostSetInfoVecs.edge_initial_length;
//    linearSpringInfoVecs.edge_final_length = hostSetInfoVecs.edge_initial_length;
    
    
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
    
    if ((int)hostSetInfoVecs.edge_initial_length.size() != E) {
        std::cerr << "FATAL: host edge_initial_length size "
                  << hostSetInfoVecs.edge_initial_length.size()
                  << " != num_edges " << E << "\n";
        std::abort();
    }


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
    auxVecs.id_bucket_expanded.resize( (generalParams.maxNodeCount));
    auxVecs.id_value_expanded.resize( (generalParams.maxNodeCount));
};






