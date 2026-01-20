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

// ============================================================================
// FIXED solveSystem() function for System.cu
// 
// Key fixes:
// 1. Uses thrust operations for rest length interpolation instead of host loop
// 2. Adds stability checks and early termination
// 3. Reduces timestep for strain application
// 4. Uses more substeps for gentler strain application
// ============================================================================

// Add this functor near the top of System.cu (after includes):
struct InterpolateRestLengthFunctor {
    double t;  // interpolation parameter [0,1]
    
    __host__ __device__
    InterpolateRestLengthFunctor(double _t) : t(_t) {}
    
    __device__
    double operator()(const thrust::tuple<double, double>& lengths) {
        double start_len = thrust::get<0>(lengths);
        double final_len = thrust::get<1>(lengths);
        return (1.0 - t) * start_len + t * final_len;
    }
};

// Replace your solveSystem() function with this version:

void System::solveSystem()
{
    // ============================================
    //              INITIALIZATION
    // ============================================
    
    BuildPrismsFromLayerFlags(generalParams, coordInfoVecs, prismInfoVecs);

    // Use a VERY small timestep for stability
    generalParams.dt = 0.0000001;  // 1e-7

    // Compute initial center
    double cx = 0.0, cy = 0.0, cz = 0.0;
    int count = 0;
    for (int i = 0; i < generalParams.maxNodeCount; i++) {
        cx += coordInfoVecs.nodeLocX[i];
        cy += coordInfoVecs.nodeLocY[i];
        cz += coordInfoVecs.nodeLocZ[i];
        count++;
    }
    if (count > 0) {
        generalParams.centerX = cx / count;
        generalParams.centerY = cy / count;
        generalParams.centerZ = cz / count;
    }

    // Compute initial volume
    ComputeVolume(generalParams, coordInfoVecs, linearSpringInfoVecs, 
                  ljInfoVecs, prismInfoVecs);
    generalParams.eq_total_volume = generalParams.current_total_volume;
    std::cout << "Initial volume = " << generalParams.eq_total_volume << std::endl;

    // ============================================
    //     INITIAL RELAXATION
    // ============================================
    
    std::cout << "\n========== INITIAL RELAXATION ==========" << std::endl;
    
    thrust::copy(linearSpringInfoVecs.edge_initial_length.begin(),
                 linearSpringInfoVecs.edge_initial_length.end(),
                 linearSpringInfoVecs.edge_rest_length.begin());

    generalParams.tol = 1e-8;
    int initial_relax_iters = relaxUntilConverged(*this);
    
    std::cout << "Initial relaxation converged in " << initial_relax_iters 
              << " iterations, E = " << linearSpringInfoVecs.linear_spring_energy << std::endl;

    storage->print_VTK_File();

    // ============================================
    //              STRAIN APPLICATION
    // ============================================
    
    int num_stages = static_cast<int>(generalParams.Tf);
    if (num_stages < 1) num_stages = 1;
    
    std::cout << "\n========== APPLYING STRAIN IN " << num_stages << " STAGE(S) ==========" << std::endl;

    // Compute basis vectors ONCE at the start
    LambdaField lambda;
    double theta_DV = 0.1931;
    double R = 1.0;  // Will be auto-detected
    
    StrainTensorGPU::computeBasisVectorsAndPathlength(
        generalParams, coordInfoVecs, lambda, theta_DV, R);

    // Flag to track simulation stability
    bool simulation_stable = true;

    for (int stage = 0; stage < num_stages && simulation_stable; stage++)
    {
        std::cout << "\n----- Stage " << stage + 1 << " of " << num_stages << " -----" << std::endl;

        // Get lambda coefficients for this stage
        StrainTensorGPU::getLambdaCoeffsForStage(
            stage,
            generalParams.lambda_iso_center_outDV,
            generalParams.lambda_iso_edge_outDV,
            generalParams.lambda_aniso_center_outDV,
            generalParams.lambda_aniso_edge_outDV,
            generalParams.lambda_iso_center_inDV,
            generalParams.lambda_iso_edge_inDV,
            generalParams.lambda_aniso_center_inDV,
            generalParams.lambda_aniso_edge_inDV);

        // Build vertex lambda values
        double frac = 1.0;
        StrainTensorGPU::buildVertexLambda(generalParams, lambda, frac);

        // Compute target rest lengths
        int layerflag = -1;  // Apply to all horizontal layers, skip vertical
        StrainTensorGPU::updateEdgeRestLengths(
            coordInfoVecs, generalParams, lambda, 
            linearSpringInfoVecs, layerflag);

        // ========================================================================
        // FIXED: Quasi-static strain application using GPU operations
        // ========================================================================
        
        int num_substeps = 10;           // More substeps for gentler application
        int relax_iters_per_substep = 10; // Fewer iters per substep, more substeps
        
        std::cout << "Applying strain quasi-statically in " << num_substeps << " substeps..." << std::endl;
        
        // Store starting and final rest lengths ON DEVICE (no host copy!)
        thrust::device_vector<double> d_start_length = linearSpringInfoVecs.edge_rest_length;
        thrust::device_vector<double> d_final_length = linearSpringInfoVecs.edge_final_length;
        
        for (int sub = 1; sub <= num_substeps && simulation_stable; sub++)
        {
            double t = static_cast<double>(sub) / static_cast<double>(num_substeps);
            
            // ========================================================================
            // GPU-accelerated interpolation (FIXED: no host-device loop!)
            // ========================================================================
            thrust::transform(
                thrust::make_zip_iterator(
                    thrust::make_tuple(d_start_length.begin(), d_final_length.begin())),
                thrust::make_zip_iterator(
                    thrust::make_tuple(d_start_length.end(), d_final_length.end())),
                linearSpringInfoVecs.edge_rest_length.begin(),
                InterpolateRestLengthFunctor(t)
            );
            
            // Relaxation loop with stability checking
            for (int relax_iter = 0; relax_iter < relax_iters_per_substep; relax_iter++) {
                Solve_Forces();
                AdvancePositions(coordInfoVecs, generalParams, domainParams);
                
                // Check for NaN in energy periodically
                if (relax_iter % 10 == 0) {
                    if (std::isnan(linearSpringInfoVecs.linear_spring_energy)) {
                        std::cout << "ERROR: NaN energy at substep " << sub 
                                  << ", iter " << relax_iter << std::endl;
                        simulation_stable = false;
                        break;
                    }
                }
            }
            
            // Progress reporting and volume check
            if (sub % 50 == 0 || sub == num_substeps) {
                ComputeVolume(generalParams, coordInfoVecs, linearSpringInfoVecs, 
                              ljInfoVecs, prismInfoVecs);
                std::cout << "  Substep " << sub << "/" << num_substeps 
                          << ", E = " << linearSpringInfoVecs.linear_spring_energy
                          << ", V = " << generalParams.current_total_volume << std::endl;
                          
                // Check for numerical problems
                if (std::isnan(linearSpringInfoVecs.linear_spring_energy) ||
                    generalParams.current_total_volume < 0 ||
                    std::isnan(generalParams.current_total_volume)) {
                    std::cout << "  ERROR: Numerical instability detected!" << std::endl;
                    std::cout << "  Attempting recovery with smaller timestep..." << std::endl;
                    
                    // Try to recover by reducing timestep
                    generalParams.dt *= 0.1;
                    if (generalParams.dt < 1e-12) {
                        std::cout << "  FATAL: Timestep too small, giving up." << std::endl;
                        simulation_stable = false;
                    }
                }
            }
        }
        
        if (!simulation_stable) {
            std::cout << "Aborting simulation due to instability at stage " << stage + 1 << std::endl;
            break;
        }
        
        // Final relaxation for this stage
        std::cout << "Final relaxation for stage " << stage + 1 << "..." << std::endl;
        generalParams.tol = 1e-8;
        int final_iters = relaxUntilConverged(*this);
        std::cout << "Final relaxation: " << final_iters << " iterations, "
                  << "E = " << linearSpringInfoVecs.linear_spring_energy << std::endl;
        
        // Check stability after final relaxation
        if (std::isnan(linearSpringInfoVecs.linear_spring_energy)) {
            std::cout << "ERROR: NaN energy after final relaxation." << std::endl;
            simulation_stable = false;
        }
        
        // Update edge_rest_length from edge_final_length for next stage
        thrust::copy(d_final_length.begin(), d_final_length.end(), 
                     linearSpringInfoVecs.edge_rest_length.begin());
        
        storage->print_VTK_File();
        
    } // end stage loop
    
    if (simulation_stable) {
        std::cout << "\n========== SIMULATION COMPLETE ==========" << std::endl;
    } else {
        std::cout << "\n========== SIMULATION FAILED DUE TO INSTABILITY ==========" << std::endl;
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
    
    prismInfoVecs.num_prisms = 0;
    prismInfoVecs.P1.clear(); prismInfoVecs.P2.clear(); prismInfoVecs.P3.clear();
    prismInfoVecs.P4.clear(); prismInfoVecs.P5.clear(); prismInfoVecs.P6.clear();

    
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






