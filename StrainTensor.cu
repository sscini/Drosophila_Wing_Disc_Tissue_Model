#include "StrainTensor.h"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace StrainTensorGPU {

/**
 * Helper function: Convert Cartesian to spherical coordinates
 */
static void cart2sph(double x, double y, double z, double& r, double& theta, double& phi) {
    r = std::sqrt(x*x + y*y + z*z);
    if (r < 1e-10) {
        theta = 0.0;
        phi = 0.0;
        return;
    }
    theta = std::acos(z / r);  // polar angle from +z axis
    phi = std::atan2(y, x);     // azimuthal angle
}

/**
 * Helper function: Convert spherical to Cartesian coordinates
 */
static void sph2cart(double r, double theta, double phi, double& x, double& y, double& z) {
    x = r * std::sin(theta) * std::cos(phi);
    y = r * std::sin(theta) * std::sin(phi);
    z = r * std::cos(theta);
}

/**
 * Helper function: Compute arc distance on sphere between two points
 */
static double arcDistOnSphere(double x1, double y1, double z1, 
                               double x2, double y2, double z2, double R) {
    double dot = (x1*x2 + y1*y2 + z1*z2) / (R*R);
    // Clamp for numerical safety
    dot = std::max(-1.0, std::min(1.0, dot));
    return R * std::acos(dot);
}

/**
 * Helper function: Normalize a 3D vector in place
 */
static void normalize3(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-10) {
        x /= len;
        y /= len;
        z /= len;
    }
}

void computeBasisVectorsAndPathlength(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LambdaField& lambda,
    double theta_DV,
    double R)
{
    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
    
    std::cout << "Computing basis vectors and pathlength for " << N << " nodes..." << std::endl;
    std::cout << "  theta_DV = " << theta_DV << " rad, R = " << R << std::endl;
    
    // Resize lambda field
    lambda.resize(N);
    
    // First pass: compute pathlength (unnormalized) and basis vectors
    double max_pathlength_inDV = 0.0;
    double max_pathlength_outDV = 0.0;
    
    for (int i = 0; i < N; i++) {
        double x = coordInfoVecs.nodeLocX[i];
        double y = coordInfoVecs.nodeLocY[i];
        double z = coordInfoVecs.nodeLocZ[i];
        double r = std::sqrt(x*x + y*y + z*z);
        
        if (r < 1e-10) r = 1e-10;  // Avoid division by zero
        
        // ========================================
        // 1. Surface normal (e_h) - points radially outward for spherical cap
        // ========================================
        lambda.e_h_x[i] = x / r;
        lambda.e_h_y[i] = y / r;
        lambda.e_h_z[i] = z / r;
        
        // ========================================
        // 2. Determine if point is inside DV region
        // ========================================
        // DV boundary is defined by |x| <= R * sin(theta_DV / 2)
        double DV_boundary = R * std::sin(theta_DV / 2.0);
        bool inDV = (std::abs(x) <= DV_boundary);
        
        // Update nodes_in_DV if it exists and has correct size
        if (coordInfoVecs.nodes_in_DV.size() == static_cast<size_t>(N)) {
            coordInfoVecs.nodes_in_DV[i] = inDV ? 1 : 0;
        }
        
        // ========================================
        // 3. Compute center point for pathlength calculation
        // ========================================
        double cx, cy, cz;
        
        if (inDV) {
            // Inside DV: center is on the DV midline (y=0 plane)
            // Project point onto y=0 plane while staying on sphere
            cx = x;
            cy = 0.0;
            double temp = R*R - x*x;
            cz = (temp > 0) ? std::sqrt(temp) : 0.0;
        } else {
            // Outside DV: center is at the DV boundary edge
            // This is a point at angle theta_DV/2 from the pole
            double sign_x = (x >= 0) ? 1.0 : -1.0;
            double theta_center = theta_DV / 2.0;
            cx = R * std::sin(theta_center) * sign_x;
            cy = 0.0;
            cz = R * std::cos(theta_center);
        }
        
        // ========================================
        // 4. Compute pathlength: arc distance from center
        // ========================================
        double pathlength = arcDistOnSphere(x, y, z, cx, cy, cz, R);
        lambda.pathlength_scaled[i] = pathlength;  // Will normalize later
        
        // Track max pathlength for normalization
        if (inDV) {
            max_pathlength_inDV = std::max(max_pathlength_inDV, pathlength);
        } else {
            max_pathlength_outDV = std::max(max_pathlength_outDV, pathlength);
        }
        
        // ========================================
        // 5. Compute e_R: in-surface radial direction pointing away from center
        // ========================================
        // e_OA = (position - center) / |position - center|
        double oax = x - cx;
        double oay = y - cy;
        double oaz = z - cz;
        normalize3(oax, oay, oaz);
        
        // Project onto surface: e_R = e_OA - (e_h · e_OA) * e_h
        double eh_dot_eoa = lambda.e_h_x[i]*oax + lambda.e_h_y[i]*oay + lambda.e_h_z[i]*oaz;
        double er_x = oax - eh_dot_eoa * lambda.e_h_x[i];
        double er_y = oay - eh_dot_eoa * lambda.e_h_y[i];
        double er_z = oaz - eh_dot_eoa * lambda.e_h_z[i];
        normalize3(er_x, er_y, er_z);
        
        lambda.e_R_x[i] = er_x;
        lambda.e_R_y[i] = er_y;
        lambda.e_R_z[i] = er_z;
        
        // ========================================
        // 6. Compute e_phi: circumferential direction = e_h × e_R
        // ========================================
        lambda.e_phi_x[i] = lambda.e_h_y[i] * er_z - lambda.e_h_z[i] * er_y;
        lambda.e_phi_y[i] = lambda.e_h_z[i] * er_x - lambda.e_h_x[i] * er_z;
        lambda.e_phi_z[i] = lambda.e_h_x[i] * er_y - lambda.e_h_y[i] * er_x;
    }
    
    // ========================================
    // Second pass: normalize pathlength separately for DV and outDV
    // ========================================
    std::cout << "  Max pathlength inDV: " << max_pathlength_inDV 
              << ", outDV: " << max_pathlength_outDV << std::endl;
    
    int count_inDV = 0, count_outDV = 0;
    for (int i = 0; i < N; i++) {
        bool inDV = false;
        if (coordInfoVecs.nodes_in_DV.size() == static_cast<size_t>(N)) {
            inDV = (hostSetInfoVecs.nodes_in_DV[i] == 1);
        } else {
            // Fallback: check x coordinate
            double x = coordInfoVecs.nodeLocX[i];
            double DV_boundary = R * std::sin(theta_DV / 2.0);
            inDV = (std::abs(x) <= DV_boundary);
        }
        
        if (inDV && max_pathlength_inDV > 1e-10) {
            lambda.pathlength_scaled[i] /= max_pathlength_inDV;
            count_inDV++;
        } else if (!inDV && max_pathlength_outDV > 1e-10) {
            lambda.pathlength_scaled[i] /= max_pathlength_outDV;
            count_outDV++;
        }
    }
    
    std::cout << "  Nodes in DV: " << count_inDV << ", nodes in outDV: " << count_outDV << std::endl;
}

void buildVertexLambda(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LambdaField& lambda,
    double frac)
{
    const int N = static_cast<int>(coordInfoVecs.nodeLocX.size());
    
    std::cout << "Building vertex lambda values for " << N << " nodes (frac = " << frac << ")..." << std::endl;
    
    // Ensure lambda vectors are properly sized
    if (lambda.lambda_RR.size() != static_cast<size_t>(N)) {
        lambda.resize(N);
    }
    
    int count_inDV = 0, count_outDV = 0;
    double avg_lambda_iso = 0, avg_lambda_aniso = 0;
    
    for (int i = 0; i < N; i++) {
        // Get normalized pathlength (should be in [0, 1])
        double p = lambda.pathlength_scaled[i];
        p = std::max(0.0, std::min(1.0, p));  // Clamp just in case
        
        // Determine if in DV region
        bool inDV = false;
        if (coordInfoVecs.nodes_in_DV.size() == static_cast<size_t>(N)) {
            inDV = (hostSetInfoVecs.nodes_in_DV[i] == 1);
        }
        
        double lambda_iso, lambda_aniso;
        
        if (inDV) {
            // Inside DV: use inDV lambda parameters
            // Linear interpolation: ?(p) = center + (edge - center) * p
            lambda_iso = generalParams.lambda_iso_center_inDV + 
                        (generalParams.lambda_iso_edge_inDV - generalParams.lambda_iso_center_inDV) * p;
            lambda_aniso = generalParams.lambda_aniso_center_inDV + 
                          (generalParams.lambda_aniso_edge_inDV - generalParams.lambda_aniso_center_inDV) * p;
            count_inDV++;
        } else {
            // Outside DV: use outDV lambda parameters
            lambda_iso = generalParams.lambda_iso_center_outDV + 
                        (generalParams.lambda_iso_edge_outDV - generalParams.lambda_iso_center_outDV) * p;
            lambda_aniso = generalParams.lambda_aniso_center_outDV + 
                          (generalParams.lambda_aniso_edge_outDV - generalParams.lambda_aniso_center_outDV) * p;
            count_outDV++;
        }
        
        // Apply fractional strain for quasi-static loading
        // ?_applied = 1 + frac * (?_target - 1)
        lambda_iso = 1.0 + frac * (lambda_iso - 1.0);
        lambda_aniso = 1.0 + frac * (lambda_aniso - 1.0);
        
        // ========================================
        // Convert to directional components (PAPER'S FORMULA!)
        // ========================================
        // ?_RR = ?_iso * ?_aniso (stretch in radial direction)
        // ?_ff = ?_iso / ?_aniso (stretch in circumferential direction) - NOTE: DIVISION!
        // ?_hh = 1.0 (no stretch in height direction)
        
        lambda.lambda_RR[i] = lambda_iso * lambda_aniso;
        lambda.lambda_phiphi[i] = lambda_iso / lambda_aniso;  // CRITICAL: division, not multiplication!
        lambda.lambda_hh[i] = 1.0;
        
        avg_lambda_iso += lambda_iso;
        avg_lambda_aniso += lambda_aniso;
    }
    
    avg_lambda_iso /= N;
    avg_lambda_aniso /= N;
    
    std::cout << "  Lambda stats: inDV=" << count_inDV << ", outDV=" << count_outDV << std::endl;
    std::cout << "  Avg lambda_iso=" << avg_lambda_iso << ", avg lambda_aniso=" << avg_lambda_aniso << std::endl;
}

void updateEdgeRestLengths(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    LambdaField& lambda,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    int layerflag)
{
    const int num_edges = coordInfoVecs.num_edges;
    
    std::cout << "Updating rest lengths for " << num_edges << " edges (layerflag=" << layerflag << ")..." << std::endl;
    
    int edges_modified = 0;
    double max_strain = 0.0;
    double avg_strain = 0.0;
    
    for (int e = 0; e < num_edges; e++) {
        // Get edge endpoint indices
        int v1 = coordInfoVecs.edges2Nodes_1[e];
        int v2 = coordInfoVecs.edges2Nodes_2[e];
        
        // Skip edges not in target layer (if filtering is enabled)
        if (layerflag >= 0) {
            // Check if edge is in the specified layer using edges_in_upperhem
            if (generalParams.edges_in_upperhem.size() > static_cast<size_t>(e)) {
                int edge_layer = hostSetInfoVecs.edges_in_upperhem[e];
                if (edge_layer != layerflag) {
                    // Keep edge at current length
                    linearSpringInfoVecs.edge_final_length[e] = linearSpringInfoVecs.edge_initial_length[e];
                    continue;
                }
            }
        }
        
        // ========================================
        // 1. Get spring vector
        // ========================================
        double dx = coordInfoVecs.nodeLocX[v2] - coordInfoVecs.nodeLocX[v1];
        double dy = coordInfoVecs.nodeLocY[v2] - coordInfoVecs.nodeLocY[v1];
        double dz = coordInfoVecs.nodeLocZ[v2] - coordInfoVecs.nodeLocZ[v1];
        
        // ========================================
        // 2. Average lambda values at endpoints
        // ========================================
        double avg_RR = 0.5 * (lambda.lambda_RR[v1] + lambda.lambda_RR[v2]);
        double avg_phiphi = 0.5 * (lambda.lambda_phiphi[v1] + lambda.lambda_phiphi[v2]);
        double avg_hh = 0.5 * (lambda.lambda_hh[v1] + lambda.lambda_hh[v2]);
        
        // ========================================
        // 3. Average and renormalize basis vectors
        // ========================================
        // e_R
        double avg_eR_x = 0.5 * (lambda.e_R_x[v1] + lambda.e_R_x[v2]);
        double avg_eR_y = 0.5 * (lambda.e_R_y[v1] + lambda.e_R_y[v2]);
        double avg_eR_z = 0.5 * (lambda.e_R_z[v1] + lambda.e_R_z[v2]);
        normalize3(avg_eR_x, avg_eR_y, avg_eR_z);
        
        // e_phi
        double avg_ephi_x = 0.5 * (lambda.e_phi_x[v1] + lambda.e_phi_x[v2]);
        double avg_ephi_y = 0.5 * (lambda.e_phi_y[v1] + lambda.e_phi_y[v2]);
        double avg_ephi_z = 0.5 * (lambda.e_phi_z[v1] + lambda.e_phi_z[v2]);
        normalize3(avg_ephi_x, avg_ephi_y, avg_ephi_z);
        
        // e_h
        double avg_eh_x = 0.5 * (lambda.e_h_x[v1] + lambda.e_h_x[v2]);
        double avg_eh_y = 0.5 * (lambda.e_h_y[v1] + lambda.e_h_y[v2]);
        double avg_eh_z = 0.5 * (lambda.e_h_z[v1] + lambda.e_h_z[v2]);
        normalize3(avg_eh_x, avg_eh_y, avg_eh_z);
        
        // ========================================
        // 4. Project spring vector onto each basis direction
        // ========================================
        double v_dot_eR = dx*avg_eR_x + dy*avg_eR_y + dz*avg_eR_z;
        double v_dot_ephi = dx*avg_ephi_x + dy*avg_ephi_y + dz*avg_ephi_z;
        double v_dot_eh = dx*avg_eh_x + dy*avg_eh_y + dz*avg_eh_z;
        
        // ========================================
        // 5. Transform spring vector: v' = ? · v
        // ? = ?_RR * (e_R ? e_R) + ?_ff * (e_f ? e_f) + ?_hh * (e_h ? e_h)
        // v' = ?_RR * (v·e_R) * e_R + ?_ff * (v·e_f) * e_f + ?_hh * (v·e_h) * e_h
        // ========================================
        double vx_p = avg_RR * v_dot_eR * avg_eR_x + 
                      avg_phiphi * v_dot_ephi * avg_ephi_x + 
                      avg_hh * v_dot_eh * avg_eh_x;
        double vy_p = avg_RR * v_dot_eR * avg_eR_y + 
                      avg_phiphi * v_dot_ephi * avg_ephi_y + 
                      avg_hh * v_dot_eh * avg_eh_y;
        double vz_p = avg_RR * v_dot_eR * avg_eR_z + 
                      avg_phiphi * v_dot_ephi * avg_ephi_z + 
                      avg_hh * v_dot_eh * avg_eh_z;
        
        // ========================================
        // 6. New rest length = |v'|
        // ========================================
        double new_length = std::sqrt(vx_p*vx_p + vy_p*vy_p + vz_p*vz_p);
        linearSpringInfoVecs.edge_final_length[e] = new_length;
        
        // Track statistics
        double initial_length = linearSpringInfoVecs.edge_initial_length[e];
        if (initial_length > 1e-10) {
            double strain = (new_length - initial_length) / initial_length;
            avg_strain += std::abs(strain);
            max_strain = std::max(max_strain, std::abs(strain));
            edges_modified++;
        }
    }
    
    if (edges_modified > 0) {
        avg_strain /= edges_modified;
    }
    
    std::cout << "  Edges modified: " << edges_modified 
              << ", avg strain: " << avg_strain * 100 << "%" 
              << ", max strain: " << max_strain * 100 << "%" << std::endl;
}

} // namespace StrainTensorGPU
