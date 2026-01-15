// BuildPrismsFromLayerFlags_fixed.cpp
// 
// This is a replacement for the BuildPrismsFromLayerFlags function in System.cu
// that properly handles hexagonal prism meshes where apical and basal layers
// have the same node count and are connected by vertical edges.
//
// COPY THIS INTO YOUR System.cu, replacing the existing BuildPrismsFromLayerFlags function

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>

// Helper: compute 6 * signed volume of tetrahedron (i,j,k,l)
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

    const double cx = vy * wz - vz * wy;
    const double cy = vz * wx - vx * wz;
    const double cz = vx * wy - vy * wx;

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

static void BuildPrismsFromLayerFlags(
    const GeneralParams& gp,
    const CoordInfoVecs& coord,
    PrismInfoVecs& prism)
{
    // Copy data to host
    thrust::host_vector<int> hLayer = gp.nodes_in_upperhem;
    thrust::host_vector<int> t1 = coord.triangles2Nodes_1;
    thrust::host_vector<int> t2 = coord.triangles2Nodes_2;
    thrust::host_vector<int> t3 = coord.triangles2Nodes_3;
    
    // Get edge connectivity to find vertical edges
    thrust::host_vector<int> e1 = coord.edges2Nodes_1;
    thrust::host_vector<int> e2 = coord.edges2Nodes_2;
    thrust::host_vector<int> edgeLayer = gp.edges_in_upperhem;

    const int Nnodes = (int)coord.nodeLocX.size();
    const int Nedges = (int)coord.num_edges;
    const int Ntriangles = coord.num_triangles;

    if (Nnodes == 0) {
        std::cout << "BuildPrisms: No nodes!" << std::endl;
        prism.num_prisms = 0;
        return;
    }

    // Copy coordinates
    std::vector<double> X(Nnodes), Y(Nnodes), Z(Nnodes);
    for (int i = 0; i < Nnodes; ++i) {
        X[i] = coord.nodeLocX[i];
        Y[i] = coord.nodeLocY[i];
        Z[i] = coord.nodeLocZ[i];
    }

    // Determine layer structure
    int maxFlag = -1;
    int minFlag = INT_MAX;
    for (int i = 0; i < Nnodes; ++i) {
        if (hLayer[i] >= 0) {
            maxFlag = std::max(maxFlag, (int)hLayer[i]);
            minFlag = std::min(minFlag, (int)hLayer[i]);
        }
    }
    
    std::cout << "BuildPrisms: Layer flags range from " << minFlag << " to " << maxFlag << std::endl;

    if (maxFlag < 1) {
        std::cout << "BuildPrisms: Need at least 2 layers!" << std::endl;
        prism.num_prisms = 0;
        return;
    }

    // Build map of vertical connections using ACTUAL vertical edges from the mesh
    // edgeLayer == -1 means vertical edge
    std::map<int, int> apicalToBasal;  // Maps apical node -> basal node
    std::map<int, int> basalToApical;  // Maps basal node -> apical node
    
    int apicalLayer = maxFlag;  // Typically 1 for 2-layer mesh
    int basalLayer = minFlag;   // Typically 0

    std::cout << "BuildPrisms: Apical layer = " << apicalLayer << ", Basal layer = " << basalLayer << std::endl;

    for (int e = 0; e < Nedges; ++e) {
        int n1 = e1[e];
        int n2 = e2[e];
        
        // Check if this is a vertical edge (connects different layers)
        if (n1 < 0 || n2 < 0 || n1 >= Nnodes || n2 >= Nnodes) continue;
        
        int L1 = hLayer[n1];
        int L2 = hLayer[n2];
        
        // Vertical edge if layers differ (or edgeLayer == -1)
        bool isVertical = (edgeLayer[e] == -1) || (L1 != L2);
        
        if (isVertical && L1 >= 0 && L2 >= 0) {
            // Determine which is apical and which is basal
            int apicalNode = -1, basalNode = -1;
            
            if (L1 == apicalLayer && L2 == basalLayer) {
                apicalNode = n1;
                basalNode = n2;
            } else if (L2 == apicalLayer && L1 == basalLayer) {
                apicalNode = n2;
                basalNode = n1;
            }
            
            if (apicalNode >= 0 && basalNode >= 0) {
                apicalToBasal[apicalNode] = basalNode;
                basalToApical[basalNode] = apicalNode;
            }
        }
    }

    std::cout << "BuildPrisms: Found " << apicalToBasal.size() << " vertical edge pairs" << std::endl;

    // If no vertical edges found, fall back to nearest-neighbor matching
    if (apicalToBasal.empty()) {
        std::cout << "BuildPrisms: No vertical edges found, using nearest-neighbor matching" << std::endl;
        
        // Collect apical and basal nodes
        std::vector<int> apicalNodes, basalNodes;
        for (int i = 0; i < Nnodes; ++i) {
            if (hLayer[i] == apicalLayer) apicalNodes.push_back(i);
            else if (hLayer[i] == basalLayer) basalNodes.push_back(i);
        }
        
        // Match by XY distance
        for (int ap : apicalNodes) {
            double bestDist = 1e30;
            int bestBasal = -1;
            for (int ba : basalNodes) {
                double dx = X[ap] - X[ba];
                double dy = Y[ap] - Y[ba];
                double dist = dx*dx + dy*dy;
                if (dist < bestDist) {
                    bestDist = dist;
                    bestBasal = ba;
                }
            }
            if (bestBasal >= 0) {
                apicalToBasal[ap] = bestBasal;
                basalToApical[bestBasal] = ap;
            }
        }
    }

    // Now build prisms from apical triangles
    std::vector<int> P1, P2, P3, P4, P5, P6;
    int skippedNoMapping = 0;
    int flippedCount = 0;

    for (int ti = 0; ti < Ntriangles; ++ti) {
        int a = t1[ti], b = t2[ti], c = t3[ti];
        if (a < 0 || b < 0 || c < 0) continue;
        if (a >= Nnodes || b >= Nnodes || c >= Nnodes) continue;

        // Check if this triangle is on the apical layer
        int La = hLayer[a], Lb = hLayer[b], Lc = hLayer[c];
        
        // Only use triangles where ALL vertices are on apical layer
        if (La != apicalLayer || Lb != apicalLayer || Lc != apicalLayer) {
            continue;
        }

        // Find corresponding basal nodes
        auto itA = apicalToBasal.find(a);
        auto itB = apicalToBasal.find(b);
        auto itC = apicalToBasal.find(c);
        
        if (itA == apicalToBasal.end() || itB == apicalToBasal.end() || itC == apicalToBasal.end()) {
            skippedNoMapping++;
            continue;
        }

        int A = itA->second;  // Basal node below a
        int B = itB->second;  // Basal node below b
        int C = itC->second;  // Basal node below c

        // Compute prism volume with current orientation
        double vol6 = computePrism6V(X, Y, Z, a, b, c, A, B, C);
        
        // If negative, flip the triangle winding (swap b and c, and B and C)
        if (vol6 < 0) {
            std::swap(b, c);
            std::swap(B, C);
            vol6 = computePrism6V(X, Y, Z, a, b, c, A, B, C);
            flippedCount++;
        }
        
        // Check individual tetrahedra
        double s1 = sixV_tet_host(X, Y, Z, b, c, A, a);
        double s2 = sixV_tet_host(X, Y, Z, b, A, C, a);
        double s3 = sixV_tet_host(X, Y, Z, A, B, C, a);
        
        // If still has negative tets, try alternative decomposition
        if (s1 < 0 || s2 < 0 || s3 < 0) {
            // Try using different apex point
            // Alternative: use basal triangle as base, apical node as apex
            double alt_s1 = sixV_tet_host(X, Y, Z, B, C, a, A);
            double alt_s2 = sixV_tet_host(X, Y, Z, B, a, c, A);
            double alt_s3 = sixV_tet_host(X, Y, Z, a, b, c, A);
            
            double alt_vol6 = alt_s1 + alt_s2 + alt_s3;
            
            if (alt_vol6 > 0 && alt_s1 >= 0 && alt_s2 >= 0 && alt_s3 >= 0) {
                // Use alternative decomposition (swap roles of apical/basal)
                std::swap(a, A);
                std::swap(b, B);
                std::swap(c, C);
            }
        }

        // Store the prism
        P1.push_back(a);
        P2.push_back(b);
        P3.push_back(c);
        P4.push_back(A);
        P5.push_back(B);
        P6.push_back(C);
    }

    prism.num_prisms = (int)P1.size();
    prism.P1 = P1;
    prism.P2 = P2;
    prism.P3 = P3;
    prism.P4 = P4;
    prism.P5 = P5;
    prism.P6 = P6;

    std::cout << "BuildPrisms: Created " << prism.num_prisms << " prisms" << std::endl;
    std::cout << "  - Flipped: " << flippedCount << std::endl;
    std::cout << "  - Skipped (no mapping): " << skippedNoMapping << std::endl;

    // Verify all prisms have positive volume
    int negativeCount = 0;
    for (int p = 0; p < prism.num_prisms; ++p) {
        double vol6 = computePrism6V(X, Y, Z, 
            P1[p], P2[p], P3[p], 
            P4[p], P5[p], P6[p]);
        if (vol6 < 0) {
            negativeCount++;
            std::cout << "  WARNING: Prism " << p << " has negative volume = " << vol6/6.0 << std::endl;
        }
    }
    
    if (negativeCount > 0) {
        std::cout << "  WARNING: " << negativeCount << " prisms have negative volume!" << std::endl;
    } else {
        std::cout << "  All prisms have positive volume." << std::endl;
    }
}