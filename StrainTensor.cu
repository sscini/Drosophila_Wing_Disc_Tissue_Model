// StrainTensor.cu
// IMPORTANT: Include SystemStructures.h FIRST to get full struct definitions
#include "SystemStructures.h"
#include "StrainTensor.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace StrainTensorGPU {

// ============================================================================
// Helper functions (file-local)
// ============================================================================

static inline void normalize3(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-10) {
        x /= len;
        y /= len;
        z /= len;
    }
}

static inline double arcDistOnSphere(double x1, double y1, double z1,
                                      double x2, double y2, double z2, double R) {
    double dot = (x1*x2 + y1*y2 + z1*z2) / (R*R);
    dot = std::max(-1.0, std::min(1.0, dot));
    return R * std::acos(dot);
}

// ============================================================================
// computeBasisVectorsAndPathlength
// ============================================================================
void computeBasisVectorsAndPathlength(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LambdaField& lambda,
    double theta_DV,
    double R_input)
{
    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
    
    std::cout << "Computing basis vectors and pathlength for " << N << " nodes..." << std::endl;
    std::cout << "  theta_DV = " << theta_DV << " rad" << std::endl;
    
    lambda.resize(N);
    
    // Copy node positions from device to host
    thrust::host_vector<double> h_nodeLocX(N), h_nodeLocY(N), h_nodeLocZ(N);
    thrust::copy(coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocX.end(), h_nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocY.end(), h_nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(), coordInfoVecs.nodeLocZ.end(), h_nodeLocZ.begin());
    
    // Auto-detect sphere radius from mesh
    double R = 0.0;
    double r_min = 1e20, r_max = 0.0;
    for (int i = 0; i < N; i++) {
        double r = std::sqrt(h_nodeLocX[i]*h_nodeLocX[i] + 
                             h_nodeLocY[i]*h_nodeLocY[i] + 
                             h_nodeLocZ[i]*h_nodeLocZ[i]);
        r_min = std::min(r_min, r);
        r_max = std::max(r_max, r);
    }
    R = r_max;
    
    std::cout << "  Auto-detected radii: min=" << r_min << ", max=" << r_max << std::endl;
    std::cout << "  Using R=" << R << " for pathlength calculation" << std::endl;
    
    if (R < 1e-10) {
        std::cerr << "ERROR: Sphere radius is essentially zero!" << std::endl;
        return;
    }
    
    double DV_boundary = R * std::sin(theta_DV / 2.0);
    std::cout << "  DV boundary (x-coordinate): |x| <= " << DV_boundary << std::endl;
    
    // First pass: compute pathlength and basis vectors
    double max_pathlength_inDV = 0.0;
    double max_pathlength_outDV = 0.0;
    int count_inDV = 0, count_outDV = 0;
    
    for (int i = 0; i < N; i++) {
        double x = h_nodeLocX[i];
        double y = h_nodeLocY[i];
        double z = h_nodeLocZ[i];
        double r = std::sqrt(x*x + y*y + z*z);
        
        if (r < 1e-10) r = 1e-10;
        
        // Surface normal (e_h)
        lambda.e_h_x[i] = x / r;
        lambda.e_h_y[i] = y / r;
        lambda.e_h_z[i] = z / r;
        
        // DV region check
        bool inDV = (std::abs(x) <= DV_boundary);
        lambda.nodes_in_DV[i] = inDV ? 1 : 0;
        
        if (inDV) count_inDV++;
        else count_outDV++;
        
        // Project to outer sphere for pathlength calculation
        double x_norm = x * R / r;
        double y_norm = y * R / r;
        double z_norm = z * R / r;
        
        // Center point
        double cx, cy, cz;
        if (inDV) {
            double temp = R*R - x_norm*x_norm;
            cx = x_norm;
            cy = 0.0;
            cz = (temp > 0) ? std::sqrt(temp) : R;
        } else {
            double sign_x = (x >= 0) ? 1.0 : -1.0;
            double theta_center = theta_DV / 2.0;
            cx = R * std::sin(theta_center) * sign_x;
            cy = 0.0;
            cz = R * std::cos(theta_center);
        }
        
        // Pathlength
        double pathlength = arcDistOnSphere(x_norm, y_norm, z_norm, cx, cy, cz, R);
        lambda.pathlength_scaled[i] = pathlength;
        
        if (inDV) {
            max_pathlength_inDV = std::max(max_pathlength_inDV, pathlength);
        } else {
            max_pathlength_outDV = std::max(max_pathlength_outDV, pathlength);
        }
        
        // e_R: radial direction in surface
        double oax = x_norm - cx;
        double oay = y_norm - cy;
        double oaz = z_norm - cz;
        normalize3(oax, oay, oaz);
        
        double eh_dot_eoa = lambda.e_h_x[i]*oax + lambda.e_h_y[i]*oay + lambda.e_h_z[i]*oaz;
        double er_x = oax - eh_dot_eoa * lambda.e_h_x[i];
        double er_y = oay - eh_dot_eoa * lambda.e_h_y[i];
        double er_z = oaz - eh_dot_eoa * lambda.e_h_z[i];
        normalize3(er_x, er_y, er_z);
        
        lambda.e_R_x[i] = er_x;
        lambda.e_R_y[i] = er_y;
        lambda.e_R_z[i] = er_z;
        
        // e_phi: circumferential direction = e_h × e_R
        lambda.e_phi_x[i] = lambda.e_h_y[i] * er_z - lambda.e_h_z[i] * er_y;
        lambda.e_phi_y[i] = lambda.e_h_z[i] * er_x - lambda.e_h_x[i] * er_z;
        lambda.e_phi_z[i] = lambda.e_h_x[i] * er_y - lambda.e_h_y[i] * er_x;
    }
    
    // Normalize pathlength
    std::cout << "  Max pathlength (before norm): inDV=" << max_pathlength_inDV 
              << ", outDV=" << max_pathlength_outDV << std::endl;
    
    for (int i = 0; i < N; i++) {
        bool inDV = (lambda.nodes_in_DV[i] == 1);
        if (inDV && max_pathlength_inDV > 1e-10) {
            lambda.pathlength_scaled[i] /= max_pathlength_inDV;
        } else if (!inDV && max_pathlength_outDV > 1e-10) {
            lambda.pathlength_scaled[i] /= max_pathlength_outDV;
        }
    }
    
    std::cout << "  Nodes in DV: " << count_inDV << ", nodes in outDV: " << count_outDV << std::endl;
}

// ============================================================================
// buildVertexLambda
// ============================================================================
void buildVertexLambda(
    GeneralParams& generalParams,
    LambdaField& lambda,
    double frac)
{
    const int N = static_cast<int>(lambda.pathlength_scaled.size());
    
    std::cout << "Building vertex lambda values for " << N << " nodes (frac = " << frac << ")..." << std::endl;
    
    int count_inDV = 0, count_outDV = 0;
    double avg_lambda_iso = 0.0, avg_lambda_aniso = 0.0;
    
    for (int i = 0; i < N; i++) {
        double p = lambda.pathlength_scaled[i];
        p = std::max(0.0, std::min(1.0, p));
        
        bool inDV = (lambda.nodes_in_DV[i] == 1);
        
        double lambda_iso, lambda_aniso;
        
        if (inDV) {
            lambda_iso = generalParams.lambda_iso_center_inDV + 
                        (generalParams.lambda_iso_edge_inDV - generalParams.lambda_iso_center_inDV) * p;
            lambda_aniso = generalParams.lambda_aniso_center_inDV + 
                          (generalParams.lambda_aniso_edge_inDV - generalParams.lambda_aniso_center_inDV) * p;
            count_inDV++;
        } else {
            lambda_iso = generalParams.lambda_iso_center_outDV + 
                        (generalParams.lambda_iso_edge_outDV - generalParams.lambda_iso_center_outDV) * p;
            lambda_aniso = generalParams.lambda_aniso_center_outDV + 
                          (generalParams.lambda_aniso_edge_outDV - generalParams.lambda_aniso_center_outDV) * p;
            count_outDV++;
        }
        
        // Fractional strain for quasi-static loading
        lambda_iso = 1.0 + frac * (lambda_iso - 1.0);
        lambda_aniso = 1.0 + frac * (lambda_aniso - 1.0);
        
        // ?_RR = ?_iso * ?_aniso, ?_ff = ?_iso / ?_aniso
        lambda.lambda_RR[i] = lambda_iso * lambda_aniso;
        lambda.lambda_phiphi[i] = lambda_iso / lambda_aniso;
        lambda.lambda_hh[i] = 1.0;
        
        avg_lambda_iso += lambda_iso;
        avg_lambda_aniso += lambda_aniso;
    }
    
    if (N > 0) {
        avg_lambda_iso /= N;
        avg_lambda_aniso /= N;
    }
    
    std::cout << "  Lambda stats: inDV=" << count_inDV << ", outDV=" << count_outDV << std::endl;
    std::cout << "  Avg lambda_iso=" << avg_lambda_iso << ", avg lambda_aniso=" << avg_lambda_aniso << std::endl;
}

// ============================================================================
// updateEdgeRestLengths
// ============================================================================
void updateEdgeRestLengths(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    LambdaField& lambda,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    int layerflag)
{
    const int num_edges = coordInfoVecs.num_edges;
    const int N = static_cast<int>(lambda.pathlength_scaled.size());
    
    std::cout << "Updating rest lengths for " << num_edges << " edges..." << std::endl;
    std::cout << "  layerflag=" << layerflag << " (>=0 means apply to that layer only, <0 means skip vertical)" << std::endl;
    
    // Copy data from device to host
    thrust::host_vector<double> h_nodeLocX(N), h_nodeLocY(N), h_nodeLocZ(N);
    thrust::copy(coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocX.end(), h_nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocY.end(), h_nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(), coordInfoVecs.nodeLocZ.end(), h_nodeLocZ.begin());
    
    thrust::host_vector<int> h_edges2Nodes_1(num_edges), h_edges2Nodes_2(num_edges);
    thrust::copy(coordInfoVecs.edges2Nodes_1.begin(), coordInfoVecs.edges2Nodes_1.end(), h_edges2Nodes_1.begin());
    thrust::copy(coordInfoVecs.edges2Nodes_2.begin(), coordInfoVecs.edges2Nodes_2.end(), h_edges2Nodes_2.begin());
    
    thrust::host_vector<double> h_initial_length(num_edges);
    thrust::copy(linearSpringInfoVecs.edge_initial_length.begin(), 
                 linearSpringInfoVecs.edge_initial_length.end(), h_initial_length.begin());
    
    // Copy edge layer flags
    thrust::host_vector<int> h_edges_layer(num_edges, 0);
    if (generalParams.edges_in_upperhem.size() >= static_cast<size_t>(num_edges)) {
        thrust::copy(generalParams.edges_in_upperhem.begin(), 
                     generalParams.edges_in_upperhem.end(), h_edges_layer.begin());
    }
    
    // Output vector
    thrust::host_vector<double> h_final_length(num_edges);
    
    int edges_modified = 0;
    int edges_skipped_vertical = 0;
    int edges_skipped_layer = 0;
    double max_strain = 0.0;
    double avg_strain = 0.0;
    
    for (int e = 0; e < num_edges; e++) {
        int v1 = h_edges2Nodes_1[e];
        int v2 = h_edges2Nodes_2[e];
        int edge_layer = h_edges_layer[e];
        
        // Bounds check
        if (v1 < 0 || v1 >= N || v2 < 0 || v2 >= N) {
            h_final_length[e] = h_initial_length[e];
            continue;
        }
        
        // =====================================================
        // CRITICAL: Skip VERTICAL edges (layer flag == -1)
        // Vertical edges connect different layers and should NOT
        // be transformed by the in-plane strain tensor!
        // =====================================================
        if (edge_layer == -1) {
            h_final_length[e] = h_initial_length[e];
            edges_skipped_vertical++;
            continue;
        }
        
        // If a specific layer is requested, skip other layers
        if (layerflag >= 0 && edge_layer != layerflag) {
            h_final_length[e] = h_initial_length[e];
            edges_skipped_layer++;
            continue;
        }
        
        // Get spring vector
        double dx = h_nodeLocX[v2] - h_nodeLocX[v1];
        double dy = h_nodeLocY[v2] - h_nodeLocY[v1];
        double dz = h_nodeLocZ[v2] - h_nodeLocZ[v1];
        
        // Average lambda values at endpoints
        double avg_RR = 0.5 * (lambda.lambda_RR[v1] + lambda.lambda_RR[v2]);
        double avg_phiphi = 0.5 * (lambda.lambda_phiphi[v1] + lambda.lambda_phiphi[v2]);
        double avg_hh = 0.5 * (lambda.lambda_hh[v1] + lambda.lambda_hh[v2]);
        
        // Average and normalize basis vectors
        double avg_eR_x = 0.5 * (lambda.e_R_x[v1] + lambda.e_R_x[v2]);
        double avg_eR_y = 0.5 * (lambda.e_R_y[v1] + lambda.e_R_y[v2]);
        double avg_eR_z = 0.5 * (lambda.e_R_z[v1] + lambda.e_R_z[v2]);
        normalize3(avg_eR_x, avg_eR_y, avg_eR_z);
        
        double avg_ephi_x = 0.5 * (lambda.e_phi_x[v1] + lambda.e_phi_x[v2]);
        double avg_ephi_y = 0.5 * (lambda.e_phi_y[v1] + lambda.e_phi_y[v2]);
        double avg_ephi_z = 0.5 * (lambda.e_phi_z[v1] + lambda.e_phi_z[v2]);
        normalize3(avg_ephi_x, avg_ephi_y, avg_ephi_z);
        
        double avg_eh_x = 0.5 * (lambda.e_h_x[v1] + lambda.e_h_x[v2]);
        double avg_eh_y = 0.5 * (lambda.e_h_y[v1] + lambda.e_h_y[v2]);
        double avg_eh_z = 0.5 * (lambda.e_h_z[v1] + lambda.e_h_z[v2]);
        normalize3(avg_eh_x, avg_eh_y, avg_eh_z);
        
        // Project spring vector onto basis
        double v_dot_eR = dx*avg_eR_x + dy*avg_eR_y + dz*avg_eR_z;
        double v_dot_ephi = dx*avg_ephi_x + dy*avg_ephi_y + dz*avg_ephi_z;
        double v_dot_eh = dx*avg_eh_x + dy*avg_eh_y + dz*avg_eh_z;
        
        // Transform: v' = ? · v
        double vx_p = avg_RR * v_dot_eR * avg_eR_x + 
                      avg_phiphi * v_dot_ephi * avg_ephi_x + 
                      avg_hh * v_dot_eh * avg_eh_x;
        double vy_p = avg_RR * v_dot_eR * avg_eR_y + 
                      avg_phiphi * v_dot_ephi * avg_ephi_y + 
                      avg_hh * v_dot_eh * avg_eh_y;
        double vz_p = avg_RR * v_dot_eR * avg_eR_z + 
                      avg_phiphi * v_dot_ephi * avg_ephi_z + 
                      avg_hh * v_dot_eh * avg_eh_z;
        
        // New rest length
        double new_length = std::sqrt(vx_p*vx_p + vy_p*vy_p + vz_p*vz_p);
        
        // Sanity check - cap maximum strain to prevent instability
        double initial_length = h_initial_length[e];
        if (initial_length > 1e-10) {
            double strain = (new_length - initial_length) / initial_length;
            
            // Cap strain at ±50% to prevent instability
            const double MAX_STRAIN = 0.5;
            if (strain > MAX_STRAIN) {
                new_length = initial_length * (1.0 + MAX_STRAIN);
            } else if (strain < -MAX_STRAIN) {
                new_length = initial_length * (1.0 - MAX_STRAIN);
            }
            
            // Recalculate strain after capping
            strain = (new_length - initial_length) / initial_length;
            avg_strain += std::abs(strain);
            max_strain = std::max(max_strain, std::abs(strain));
            edges_modified++;
        }
        
        h_final_length[e] = new_length;
    }
    
    // Copy results back to device
    thrust::copy(h_final_length.begin(), h_final_length.end(), 
                 linearSpringInfoVecs.edge_final_length.begin());
    
    if (edges_modified > 0) {
        avg_strain /= edges_modified;
    }
    
    std::cout << "  Edges modified: " << edges_modified << std::endl;
    std::cout << "  Edges skipped (vertical): " << edges_skipped_vertical << std::endl;
    std::cout << "  Edges skipped (wrong layer): " << edges_skipped_layer << std::endl;
    std::cout << "  Avg strain: " << avg_strain * 100 << "%, max strain: " << max_strain * 100 << "%" << std::endl;
}
// ============================================================================
// PROGRESSIVE STRAIN TESTING for getLambdaCoeffsForStage()
// 
// Use this to systematically find the maximum stable strain level
// ============================================================================

void StrainTensorGPU::getLambdaCoeffsForStage(
    int stage,
    double& iso_center_outDV, double& iso_edge_outDV,
    double& aniso_center_outDV, double& aniso_edge_outDV,
    double& iso_center_inDV, double& iso_edge_inDV,
    double& aniso_center_inDV, double& aniso_edge_inDV)
{
    // =====================================================================
    // STRAIN TESTING CONFIGURATION
    // Change TEST_LEVEL to test different strain amounts
    // =====================================================================
    
    // Test levels:
    // 0 = No strain (all 1.0) - VERIFIED WORKING
    // 1 = 2% uniform isotropic
    // 2 = 5% uniform isotropic  
    // 3 = 10% uniform isotropic
    // 4 = 15% uniform isotropic
    // 5 = 20% uniform isotropic (close to paper values)
    // 6 = Paper values with reduced magnitude (50%)
    // 7 = Full paper values
    
    const int TEST_LEVEL = 1;  // <-- CHANGE THIS TO TEST DIFFERENT LEVELS
    
    switch (TEST_LEVEL) {
        case 0:
            // No strain - verified stable
            iso_center_outDV = 1.0;
            iso_edge_outDV = 1.0;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.0;
            iso_edge_inDV = 1.0;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 0: No strain (all lambda = 1.0)" << std::endl;
            break;
            
        case 1:
            // 2% uniform isotropic growth
            iso_center_outDV = 1.02;
            iso_edge_outDV = 1.02;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.02;
            iso_edge_inDV = 1.02;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 1: 2% uniform isotropic" << std::endl;
            break;
            
        case 2:
            // 5% uniform isotropic growth
            iso_center_outDV = 1.05;
            iso_edge_outDV = 1.05;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.05;
            iso_edge_inDV = 1.05;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 2: 5% uniform isotropic" << std::endl;
            break;
            
        case 3:
            // 10% uniform isotropic growth
            iso_center_outDV = 1.10;
            iso_edge_outDV = 1.10;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.10;
            iso_edge_inDV = 1.10;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 3: 10% uniform isotropic" << std::endl;
            break;
            
        case 4:
            // 15% uniform isotropic growth
            iso_center_outDV = 1.15;
            iso_edge_outDV = 1.15;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.15;
            iso_edge_inDV = 1.15;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 4: 15% uniform isotropic" << std::endl;
            break;
            
        case 5:
            // 20% uniform isotropic growth (similar magnitude to paper)
            iso_center_outDV = 1.20;
            iso_edge_outDV = 1.20;
            aniso_center_outDV = 1.0;
            aniso_edge_outDV = 1.0;
            iso_center_inDV = 1.20;
            iso_edge_inDV = 1.20;
            aniso_center_inDV = 1.0;
            aniso_edge_inDV = 1.0;
            std::cout << "TEST_LEVEL 5: 20% uniform isotropic" << std::endl;
            break;
            
        case 6:
            // Paper values at 50% magnitude
            // Original stage 0 values, but (lambda - 1) * 0.5 + 1
            iso_center_outDV = 1.0 + (1.20789 - 1.0) * 0.5;   // = 1.104
            iso_edge_outDV = 1.0 + (1.08383 - 1.0) * 0.5;     // = 1.042
            aniso_center_outDV = 1.0 + (1.01054 - 1.0) * 0.5; // = 1.005
            aniso_edge_outDV = 1.0 + (1.07226 - 1.0) * 0.5;   // = 1.036
            iso_center_inDV = 1.0 + (1.18401 - 1.0) * 0.5;    // = 1.092
            iso_edge_inDV = 1.0 + (1.08552 - 1.0) * 0.5;      // = 1.043
            aniso_center_inDV = 1.0 + (1.03128 - 1.0) * 0.5;  // = 1.016
            aniso_edge_inDV = 1.0 + (0.90224 - 1.0) * 0.5;    // = 0.951
            std::cout << "TEST_LEVEL 6: Paper values at 50% magnitude" << std::endl;
            break;
            
        case 7:
        default:
            // Full paper values (original hardcoded values for stage 0)
            iso_center_outDV = 1.20789;
            iso_edge_outDV = 1.08383;
            aniso_center_outDV = 1.01054;
            aniso_edge_outDV = 1.07226;
            iso_center_inDV = 1.18401;
            iso_edge_inDV = 1.08552;
            aniso_center_inDV = 1.03128;
            aniso_edge_inDV = 0.90224;
            std::cout << "TEST_LEVEL 7: Full paper values (stage 0)" << std::endl;
            break;
    }
    
    // Print the actual values being used
    std::cout << "  outDV: iso_center=" << iso_center_outDV 
              << ", iso_edge=" << iso_edge_outDV
              << ", aniso_center=" << aniso_center_outDV 
              << ", aniso_edge=" << aniso_edge_outDV << std::endl;
    std::cout << "  inDV:  iso_center=" << iso_center_inDV 
              << ", iso_edge=" << iso_edge_inDV
              << ", aniso_center=" << aniso_center_inDV 
              << ", aniso_edge=" << aniso_edge_inDV << std::endl;
}


// ============================================================================
// EXPECTED RESULTS FOR EACH TEST LEVEL
// ============================================================================
//
// Level 0 (no strain):     E ˜ 0.007, V = 273962 ? VERIFIED
// Level 1 (2% iso):        E should increase slightly, V stable
// Level 2 (5% iso):        E increases more, V stable
// Level 3 (10% iso):       E increases significantly, V stable
// Level 4 (15% iso):       Getting close to instability threshold
// Level 5 (20% iso):       May need more substeps/smaller dt
// Level 6 (50% paper):     Tests spatial variation
// Level 7 (full paper):    Full strain - requires all stability fixes
//
// If a level fails (NaN, negative volume), go back one level and:
// 1. Increase num_substeps (try 100, 200, 500)
// 2. Decrease dt (try 1e-8, 1e-9)
// 3. Check if anisotropic strain causes issues (test iso-only first)
// ============================================================================
}

} // namespace StrainTensorGPU
