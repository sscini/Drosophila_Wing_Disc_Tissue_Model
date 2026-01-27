// StrainTensor.cu  — DV-aware basis construction and rest-length projection
//
// Spontaneous strain tensor ? = ?_rr (e_R ? e_R) + ?_ff (e_f ? e_f) + ?_hh (e_h ? e_h)
// with ?_rr = ?_iso * ?_aniso,  ?_ff = ?_iso / ?_aniso,  ?_hh = 1.
//
// Inside the DV stripe:
//   - e_R is the *radial* direction away from the line segment joining O_D and O_V
//     (the DV boundary line).
//   - Only radial strain is applied: ?_rr as scheduled, ?_pp = 1, ?_ss = 1.
// Outside the DV stripe:
//   - e_R is constructed as before using OA from OD (dorsal) or OV (ventral),
//     with side decided by a fixed in-plane vector perpendicular to the DV axis.
//
// Vertical/pillar edges are flagged as -1 and are not altered in the rest-length update.

#include "StrainTensor.h"
#include "SystemStructures.h"
#include "System.h"

#include <cuda_runtime.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <cmath>

#define BLOCK_SZ 256

// ------------ tuple helpers -----------------
template<int I>  __host__ __device__
inline double c(const CVec3& v) { return thrust::get<I>(v); }

__host__ __device__
inline Mat_3x3 outer(const CVec3& a, const CVec3& b) {
    return Mat_3x3(
        CVec3( c<0>(a)*c<0>(b), c<0>(a)*c<1>(b), c<0>(a)*c<2>(b) ),
        CVec3( c<1>(a)*c<0>(b), c<1>(a)*c<1>(b), c<1>(a)*c<2>(b) ),
        CVec3( c<2>(a)*c<0>(b), c<2>(a)*c<1>(b), c<2>(a)*c<2>(b) )
    );
}

__host__ __device__
inline void axpy(double s, const Mat_3x3& A, Mat_3x3& C) {
    // row 0
    {
        CVec3 r = thrust::get<0>(C);
        CVec3 a = thrust::get<0>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<0>(C)   = r;
    }
    // row 1
    {
        CVec3 r = thrust::get<1>(C);
        CVec3 a = thrust::get<1>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<1>(C)   = r;
    }
    // row 2
    {
        CVec3 r = thrust::get<2>(C);
        CVec3 a = thrust::get<2>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<2>(C)   = r;
    }
}

__host__ __device__
inline CVec3 matVec(const Mat_3x3& M, const CVec3& v) {
    return CVec3(
        c<0>( thrust::get<0>(M) )*c<0>(v) + c<1>( thrust::get<0>(M) )*c<1>(v) + c<2>( thrust::get<0>(M) )*c<2>(v),
        c<0>( thrust::get<1>(M) )*c<0>(v) + c<1>( thrust::get<1>(M) )*c<1>(v) + c<2>( thrust::get<1>(M) )*c<2>(v),
        c<0>( thrust::get<2>(M) )*c<0>(v) + c<1>( thrust::get<2>(M) )*c<1>(v) + c<2>( thrust::get<2>(M) )*c<2>(v)
    );
}

__host__ __device__
inline double norm3(const CVec3& v) {
    return sqrt(thrust::get<0>(v)*thrust::get<0>(v) +
                thrust::get<1>(v)*thrust::get<1>(v) +
                thrust::get<2>(v)*thrust::get<2>(v)) + 1e-14;
}

__host__ __device__
inline CVec3 normalize(const CVec3& v) {
    double n = norm3(v);
    return CVec3( c<0>(v)/n, c<1>(v)/n, c<2>(v)/n );
}

__host__ __device__
inline CVec3 cross(const CVec3& a, const CVec3& b) {
    return CVec3(
        c<1>(a)*c<2>(b) - c<2>(a)*c<1>(b),
        c<2>(a)*c<0>(b) - c<0>(a)*c<2>(b),
        c<0>(a)*c<1>(b) - c<1>(a)*c<0>(b)
    );
}

__host__ __device__
inline double dot3(const CVec3& a, const CVec3& b) {
    return c<0>(a)*c<0>(b) + c<1>(a)*c<1>(b) + c<2>(a)*c<2>(b);
}


// ============================================================================
// Region Classification
// ============================================================================

enum RegionType {
    REGION_DORSAL = 0,
    REGION_DV = 1,
    REGION_VENTRAL = 2
};

static RegionType classifyNodeRegion(double y, double DV_half_width) {
    if (y < -DV_half_width) {
        return REGION_DORSAL;
    } else if (y > DV_half_width) {
        return REGION_VENTRAL;
    } else {
        return REGION_DV;
    }
}

// ============================================================================
// Structure to hold region-specific origins
// ============================================================================

struct RegionOrigins {
    double OD_x, OD_y, OD_z;
    double OV_x, OV_y, OV_z;
    double center_x, center_y, center_z;
    double DV_half_width;
    double max_r_dorsal;
    double max_r_ventral;
    double max_rho_DV;
};

// ============================================================================
// Compute Basis Vectors with Region-Specific Origins
// Works entirely with device vectors - no HostSetInfoVecs
// ============================================================================

void StrainTensorGPU::computeBasisVectorsWithDVSeparation(
    GeneralParams& params,
    CoordInfoVecs& coords,
    double theta_DV,
    double R)
{
    int N = params.maxNodeCount;
    
    // ========================================================================
    // STEP 1: Copy device vectors to local std::vectors for processing
    // ========================================================================
    
    std::vector<double> h_nodeLocX(N);
    std::vector<double> h_nodeLocY(N);
    std::vector<double> h_nodeLocZ(N);
    
    thrust::copy(coords.nodeLocX.begin(), coords.nodeLocX.begin() + N, h_nodeLocX.begin());
    thrust::copy(coords.nodeLocY.begin(), coords.nodeLocY.begin() + N, h_nodeLocY.begin());
    thrust::copy(coords.nodeLocZ.begin(), coords.nodeLocZ.begin() + N, h_nodeLocZ.begin());
    
    // ========================================================================
    // STEP 2: Compute region origins
    // ========================================================================
    
    RegionOrigins origins;
    origins.DV_half_width = R * sin(theta_DV);
    
    // Compute mesh center
    double sum_x = 0, sum_y = 0, sum_z = 0;
    for (int i = 0; i < N; ++i) {
        sum_x += h_nodeLocX[i];
        sum_y += h_nodeLocY[i];
        sum_z += h_nodeLocZ[i];
    }
    origins.center_x = sum_x / N;
    origins.center_y = sum_y / N;
    origins.center_z = sum_z / N;
    
    // Compute dorsal and ventral centroids
    double dorsal_sum_x = 0, dorsal_sum_y = 0, dorsal_sum_z = 0;
    double ventral_sum_x = 0, ventral_sum_y = 0, ventral_sum_z = 0;
    int dorsal_count = 0, ventral_count = 0;
    
    for (int i = 0; i < N; ++i) {
        double y = h_nodeLocY[i];
        
        if (y < -origins.DV_half_width) {
            dorsal_sum_x += h_nodeLocX[i];
            dorsal_sum_y += h_nodeLocY[i];
            dorsal_sum_z += h_nodeLocZ[i];
            dorsal_count++;
        } else if (y > origins.DV_half_width) {
            ventral_sum_x += h_nodeLocX[i];
            ventral_sum_y += h_nodeLocY[i];
            ventral_sum_z += h_nodeLocZ[i];
            ventral_count++;
        }
    }
    
    // Compute OD (dorsal origin)
    if (dorsal_count > 0) {
        origins.OD_x = dorsal_sum_x / dorsal_count;
        origins.OD_y = dorsal_sum_y / dorsal_count;
        origins.OD_z = dorsal_sum_z / dorsal_count;
    } else {
        origins.OD_x = origins.center_x;
        origins.OD_y = -origins.DV_half_width - R * 0.25;
        origins.OD_z = origins.center_z;
    }
    
    // Compute OV (ventral origin)
    if (ventral_count > 0) {
        origins.OV_x = ventral_sum_x / ventral_count;
        origins.OV_y = ventral_sum_y / ventral_count;
        origins.OV_z = ventral_sum_z / ventral_count;
    } else {
        origins.OV_x = origins.center_x;
        origins.OV_y = origins.DV_half_width + R * 0.25;
        origins.OV_z = origins.center_z;
    }
    
    // Compute maximum radial distances for normalization
    origins.max_r_dorsal = 0;
    origins.max_r_ventral = 0;
    origins.max_rho_DV = 0;
    
    for (int i = 0; i < N; ++i) {
        double x = h_nodeLocX[i];
        double y = h_nodeLocY[i];
        double z = h_nodeLocZ[i];
        
        RegionType region = classifyNodeRegion(y, origins.DV_half_width);
        
        if (region == REGION_DORSAL) {
            double dx = x - origins.OD_x;
            double dy = y - origins.OD_y;
            double dz = z - origins.OD_z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            origins.max_r_dorsal = std::max(origins.max_r_dorsal, r);
        } else if (region == REGION_VENTRAL) {
            double dx = x - origins.OV_x;
            double dy = y - origins.OV_y;
            double dz = z - origins.OV_z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            origins.max_r_ventral = std::max(origins.max_r_ventral, r);
        } else {
            double rho = fabs(y);
            origins.max_rho_DV = std::max(origins.max_rho_DV, rho);
        }
    }
    
    // Ensure non-zero max distances
    if (origins.max_r_dorsal < 1e-10) origins.max_r_dorsal = R;
    if (origins.max_r_ventral < 1e-10) origins.max_r_ventral = R;
    if (origins.max_rho_DV < 1e-10) origins.max_rho_DV = origins.DV_half_width;
    
    // Print region origins
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DV BOUNDARY SEPARATION - Region Origins" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "DV half-width: " << origins.DV_half_width << std::endl;
    std::cout << "Mesh center: (" << origins.center_x << ", " 
              << origins.center_y << ", " << origins.center_z << ")" << std::endl;
    std::cout << "OD (Dorsal origin): (" << origins.OD_x << ", " 
              << origins.OD_y << ", " << origins.OD_z << ")" << std::endl;
    std::cout << "OV (Ventral origin): (" << origins.OV_x << ", " 
              << origins.OV_y << ", " << origins.OV_z << ")" << std::endl;
    std::cout << "Max radial distances:" << std::endl;
    std::cout << "  Dorsal: " << origins.max_r_dorsal << std::endl;
    std::cout << "  Ventral: " << origins.max_r_ventral << std::endl;
    std::cout << "  DV (rho): " << origins.max_rho_DV << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // ========================================================================
    // STEP 3: Compute basis vectors and store in local vectors
    // ========================================================================
    
    std::vector<double> h_pathlength(N);
    std::vector<double> h_eR_x(N), h_eR_y(N), h_eR_z(N);
    std::vector<double> h_ePhi_x(N), h_ePhi_y(N), h_ePhi_z(N);
    std::vector<double> h_eH_x(N), h_eH_y(N), h_eH_z(N);
    std::vector<int> h_nodes_in_DV(N);
    
    int count_dorsal = 0, count_ventral = 0, count_DV = 0;
    
    for (int i = 0; i < N; ++i) {
        double x = h_nodeLocX[i];
        double y = h_nodeLocY[i];
        double z = h_nodeLocZ[i];
        
        RegionType region = classifyNodeRegion(y, origins.DV_half_width);
        
        double radial_x, radial_y, radial_z;
        double pathlength;
        
        if (region == REGION_DV) {
            h_nodes_in_DV[i] = 1;
            count_DV++;
            
            // Radial direction is purely in y (perpendicular to y=0)
            radial_x = 0;
            radial_y = (y >= 0) ? 1.0 : -1.0;
            radial_z = 0;
            
            double rho = fabs(y);
            pathlength = rho / origins.max_rho_DV;
            
        } else if (region == REGION_DORSAL) {
            h_nodes_in_DV[i] = 0;
            count_dorsal++;
            
            double dx = x - origins.OD_x;
            double dy = y - origins.OD_y;
            double dz = z - origins.OD_z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (r > 1e-10) {
                radial_x = dx / r;
                radial_y = dy / r;
                radial_z = dz / r;
            } else {
                radial_x = 0;
                radial_y = -1;
                radial_z = 0;
            }
            
            pathlength = r / origins.max_r_dorsal;
            
        } else {
            h_nodes_in_DV[i] = 0;
            count_ventral++;
            
            double dx = x - origins.OV_x;
            double dy = y - origins.OV_y;
            double dz = z - origins.OV_z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (r > 1e-10) {
                radial_x = dx / r;
                radial_y = dy / r;
                radial_z = dz / r;
            } else {
                radial_x = 0;
                radial_y = 1;
                radial_z = 0;
            }
            
            pathlength = r / origins.max_r_ventral;
        }
        
        // Store e_R
        h_eR_x[i] = radial_x;
        h_eR_y[i] = radial_y;
        h_eR_z[i] = radial_z;
        h_pathlength[i] = pathlength;
        
        // Compute surface normal (e_h)
        double nx = x - origins.center_x;
        double ny = y - origins.center_y;
        double nz = z - origins.center_z;
        double n_mag = sqrt(nx*nx + ny*ny + nz*nz);
        
        if (n_mag > 1e-10) {
            nx /= n_mag;
            ny /= n_mag;
            nz /= n_mag;
        } else {
            nx = 0; ny = 0; nz = 1;
        }
        
        h_eH_x[i] = nx;
        h_eH_y[i] = ny;
        h_eH_z[i] = nz;
        
        // ============================================================================
        // FIX: Orthonormalize Basis Vectors Using Gram-Schmidt
        //
        // The problem: e_R (radial from region origin) and e_h (surface normal from
        // mesh center) are NOT perpendicular on a spherical dome.
        //
        // The fix: Use Gram-Schmidt orthonormalization to ensure the basis is
        // orthonormal before assembling the strain tensor.
        //
        // In StrainTensor.cu, find the section in computeBasisVectorsWithDVSeparation()
        // where e_phi is computed (around lines 383-414), and REPLACE with this code:
        // ============================================================================
        
                // ====================================================================
                // GRAM-SCHMIDT ORTHONORMALIZATION
                // 
                // We want an orthonormal basis {e_R, e_phi, e_h} where:
                //   - e_R is the primary direction (radial from region origin)
                //   - e_h is perpendicular to e_R and roughly aligned with surface normal
                //   - e_phi is perpendicular to both (tangential/circumferential)
                //
                // Process:
                //   1. Start with e_R (already unit length)
                //   2. Orthogonalize e_h against e_R: e_h' = e_h - (e_h·e_R)*e_R
                //   3. Normalize e_h'
                //   4. Compute e_phi = e_h' × e_R (guaranteed perpendicular to both)
                // ====================================================================
                
                // Step 1: e_R is already computed and normalized (radial_x, radial_y, radial_z)
                // Store e_R
                h_eR_x[i] = radial_x;
                h_eR_y[i] = radial_y;
                h_eR_z[i] = radial_z;
                h_pathlength[i] = pathlength;
                
                // Step 2: Compute initial surface normal (e_h_raw)
                nx = x - origins.center_x;
                ny = y - origins.center_y;
                nz = z - origins.center_z;
                n_mag = sqrt(nx*nx + ny*ny + nz*nz);
                
                if (n_mag > 1e-10) {
                    nx /= n_mag;
                    ny /= n_mag;
                    nz /= n_mag;
                } else {
                    nx = 0; ny = 0; nz = 1;
                }
                
                // Step 3: Orthogonalize e_h against e_R using Gram-Schmidt
                // e_h_orth = e_h_raw - (e_h_raw · e_R) * e_R
                double dot_h_R = nx*radial_x + ny*radial_y + nz*radial_z;
                
                double hx = nx - dot_h_R * radial_x;
                double hy = ny - dot_h_R * radial_y;
                double hz = nz - dot_h_R * radial_z;
                
                // Normalize e_h
                double h_mag = sqrt(hx*hx + hy*hy + hz*hz);
                
                if (h_mag > 1e-10) {
                    hx /= h_mag;
                    hy /= h_mag;
                    hz /= h_mag;
                } else {
                    // e_h was parallel to e_R, need to pick an arbitrary perpendicular direction
                    // Find a vector not parallel to e_R
                    if (fabs(radial_x) < 0.9) {
                        hx = 1; hy = 0; hz = 0;
                    } else {
                        hx = 0; hy = 1; hz = 0;
                    }
                    // Orthogonalize
                    double dot = hx*radial_x + hy*radial_y + hz*radial_z;
                    hx -= dot * radial_x;
                    hy -= dot * radial_y;
                    hz -= dot * radial_z;
                    // Normalize
                    h_mag = sqrt(hx*hx + hy*hy + hz*hz);
                    if (h_mag > 1e-10) {
                        hx /= h_mag;
                        hy /= h_mag;
                        hz /= h_mag;
                    }
                }
                
                h_eH_x[i] = hx;
                h_eH_y[i] = hy;
                h_eH_z[i] = hz;
                
                // Step 4: Compute e_phi = e_h × e_R (cross product)
                // This is guaranteed perpendicular to both e_h and e_R
                double phi_x = hy * radial_z - hz * radial_y;
                double phi_y = hz * radial_x - hx * radial_z;
                double phi_z = hx * radial_y - hy * radial_x;
                
                // Normalize e_phi (should already be unit length if e_h and e_R are orthonormal)
                double phi_mag = sqrt(phi_x*phi_x + phi_y*phi_y + phi_z*phi_z);
                
                if (phi_mag > 1e-10) {
                    phi_x /= phi_mag;
                    phi_y /= phi_mag;
                    phi_z /= phi_mag;
                } else {
                    // Degenerate case - shouldn't happen with proper orthogonalization
                    phi_x = 0; phi_y = 0; phi_z = 1;
                }
                
                h_ePhi_x[i] = phi_x;
                h_ePhi_y[i] = phi_y;
                h_ePhi_z[i] = phi_z;
                
                // ====================================================================
                // VERIFICATION (optional, can remove after testing)
                // ====================================================================
                #ifdef DEBUG_BASIS_VECTORS
                // Check orthonormality
                double dot_R_phi = radial_x*phi_x + radial_y*phi_y + radial_z*phi_z;
                double dot_R_h = radial_x*hx + radial_y*hy + radial_z*hz;
                double dot_phi_h = phi_x*hx + phi_y*hy + phi_z*hz;
                
                if (fabs(dot_R_phi) > 1e-10 || fabs(dot_R_h) > 1e-10 || fabs(dot_phi_h) > 1e-10) {
                    std::cout << "WARNING: Node " << i << " basis not orthonormal!" << std::endl;
                    std::cout << "  e_R·e_phi = " << dot_R_phi << std::endl;
                    std::cout << "  e_R·e_h = " << dot_R_h << std::endl;
                    std::cout << "  e_phi·e_h = " << dot_phi_h << std::endl;
                }
                #endif
        
        
        // ============================================================================
        // COMPLETE REPLACEMENT for the basis vector computation loop in 
        // computeBasisVectorsWithDVSeparation() (lines ~294-415 in StrainTensor.cu)
        //
        // Replace the entire loop body with this:
        // ============================================================================
        
        /*
            for (int i = 0; i < N; ++i) {
                double x = h_nodeLocX[i];
                double y = h_nodeLocY[i];
                double z = h_nodeLocZ[i];
                
                RegionType region = classifyNodeRegion(y, origins.DV_half_width);
                
                double radial_x, radial_y, radial_z;
                double pathlength;
                
                if (region == REGION_DV) {
                    h_nodes_in_DV[i] = 1;
                    count_DV++;
                    
                    radial_x = 0;
                    radial_y = (y >= 0) ? 1.0 : -1.0;
                    radial_z = 0;
                    
                    double rho = fabs(y);
                    pathlength = rho / origins.max_rho_DV;
                    
                } else if (region == REGION_DORSAL) {
                    h_nodes_in_DV[i] = 0;
                    count_dorsal++;
                    
                    double dx = x - origins.OD_x;
                    double dy = y - origins.OD_y;
                    double dz = z - origins.OD_z;
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    if (r > 1e-10) {
                        radial_x = dx / r;
                        radial_y = dy / r;
                        radial_z = dz / r;
                    } else {
                        radial_x = 0;
                        radial_y = -1;
                        radial_z = 0;
                    }
                    
                    pathlength = r / origins.max_r_dorsal;
                    
                } else { // REGION_VENTRAL
                    h_nodes_in_DV[i] = 0;
                    count_ventral++;
                    
                    double dx = x - origins.OV_x;
                    double dy = y - origins.OV_y;
                    double dz = z - origins.OV_z;
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    
                    if (r > 1e-10) {
                        radial_x = dx / r;
                        radial_y = dy / r;
                        radial_z = dz / r;
                    } else {
                        radial_x = 0;
                        radial_y = 1;
                        radial_z = 0;
                    }
                    
                    pathlength = r / origins.max_r_ventral;
                }
                
                // Store e_R (radial direction)
                h_eR_x[i] = radial_x;
                h_eR_y[i] = radial_y;
                h_eR_z[i] = radial_z;
                h_pathlength[i] = pathlength;
                
                // Compute raw surface normal
                double nx = x - origins.center_x;
                double ny = y - origins.center_y;
                double nz = z - origins.center_z;
                double n_mag = sqrt(nx*nx + ny*ny + nz*nz);
                
                if (n_mag > 1e-10) {
                    nx /= n_mag;
                    ny /= n_mag;
                    nz /= n_mag;
                } else {
                    nx = 0; ny = 0; nz = 1;
                }
                
                // GRAM-SCHMIDT: Orthogonalize e_h against e_R
                double dot_h_R = nx*radial_x + ny*radial_y + nz*radial_z;
                
                double hx = nx - dot_h_R * radial_x;
                double hy = ny - dot_h_R * radial_y;
                double hz = nz - dot_h_R * radial_z;
                
                double h_mag = sqrt(hx*hx + hy*hy + hz*hz);
                
                if (h_mag > 1e-10) {
                    hx /= h_mag;
                    hy /= h_mag;
                    hz /= h_mag;
                } else {
                    // e_h parallel to e_R - pick arbitrary perpendicular
                    if (fabs(radial_x) < 0.9) {
                        hx = 1; hy = 0; hz = 0;
                    } else {
                        hx = 0; hy = 1; hz = 0;
                    }
                    double dot = hx*radial_x + hy*radial_y + hz*radial_z;
                    hx -= dot * radial_x;
                    hy -= dot * radial_y;
                    hz -= dot * radial_z;
                    h_mag = sqrt(hx*hx + hy*hy + hz*hz);
                    if (h_mag > 1e-10) {
                        hx /= h_mag; hy /= h_mag; hz /= h_mag;
                    }
                }
                
                h_eH_x[i] = hx;
                h_eH_y[i] = hy;
                h_eH_z[i] = hz;
                
                // e_phi = e_h × e_R (perpendicular to both)
                double phi_x = hy * radial_z - hz * radial_y;
                double phi_y = hz * radial_x - hx * radial_z;
                double phi_z = hx * radial_y - hy * radial_x;
                
                double phi_mag = sqrt(phi_x*phi_x + phi_y*phi_y + phi_z*phi_z);
                if (phi_mag > 1e-10) {
                    phi_x /= phi_mag;
                    phi_y /= phi_mag;
                    phi_z /= phi_mag;
                }
                
                h_ePhi_x[i] = phi_x;
                h_ePhi_y[i] = phi_y;
                h_ePhi_z[i] = phi_z;
            }
        */
    }
    
    // Print region statistics
    std::cout << "\nNode Region Distribution:" << std::endl;
    std::cout << "  Dorsal nodes: " << count_dorsal << std::endl;
    std::cout << "  DV boundary nodes: " << count_DV << std::endl;
    std::cout << "  Ventral nodes: " << count_ventral << std::endl;
    std::cout << "  Total: " << (count_dorsal + count_DV + count_ventral) << std::endl;
    
    // ========================================================================
    // STEP 4: Copy computed values to device vectors in CoordInfoVecs
    // ========================================================================
    
    // Resize device vectors
    coords.pathlength_scaled.resize(N);
    coords.e_R_x.resize(N);
    coords.e_R_y.resize(N);
    coords.e_R_z.resize(N);
    coords.e_phi_x.resize(N);
    coords.e_phi_y.resize(N);
    coords.e_phi_z.resize(N);
    coords.e_h_x.resize(N);
    coords.e_h_y.resize(N);
    coords.e_h_z.resize(N);
    
    // Copy to device
    thrust::copy(h_pathlength.begin(), h_pathlength.end(), coords.pathlength_scaled.begin());
    thrust::copy(h_eR_x.begin(), h_eR_x.end(), coords.e_R_x.begin());
    thrust::copy(h_eR_y.begin(), h_eR_y.end(), coords.e_R_y.begin());
    thrust::copy(h_eR_z.begin(), h_eR_z.end(), coords.e_R_z.begin());
    thrust::copy(h_ePhi_x.begin(), h_ePhi_x.end(), coords.e_phi_x.begin());
    thrust::copy(h_ePhi_y.begin(), h_ePhi_y.end(), coords.e_phi_y.begin());
    thrust::copy(h_ePhi_z.begin(), h_ePhi_z.end(), coords.e_phi_z.begin());
    thrust::copy(h_eH_x.begin(), h_eH_x.end(), coords.e_h_x.begin());
    thrust::copy(h_eH_y.begin(), h_eH_y.end(), coords.e_h_y.begin());
    thrust::copy(h_eH_z.begin(), h_eH_z.end(), coords.e_h_z.begin());
    
    // Copy nodes_in_DV to GeneralParams device vector
    params.nodes_in_DV.resize(N);
    thrust::copy(h_nodes_in_DV.begin(), h_nodes_in_DV.end(), params.nodes_in_DV.begin());
    
    std::cout << "Basis vectors computed and copied to device. N = " << N << std::endl;
}


//// ============================================================================ Removed on 01/25/26 by nav 
//// Mark DV stripe (independent of layer).
//// Stripe is |x - centerX| <= R * sin(theta_DV/2).
//// ============================================================================
//
//__global__
//void k_markDVstripe(int N,
//                    const double* x,
//                    double centerX, double R, double thetaDV,
//                    const int* /*isUpper, unused for gating*/,
//                    int* DVflag)
//{
//    int i = blockIdx.x*blockDim.x + threadIdx.x;
//    if (i>=N) return;
//
//    double halfw = R * sin(0.5*thetaDV);
//    DVflag[i] = (fabs(x[i] - centerX) <= halfw) ? 1 : 0;
//}

// ============================================================================
// Build local basis vectors (e_h, e_R, e_phi) at every vertex
// DV-aware variant:
//   - inside DV stripe: e_R is radial away from the O_D–O_V line segment
//   - outside stripe:   e_R from OD/OV edge-centers on the appropriate side
// ============================================================================

// ============================================================================
// Compute Basis Vectors with Region-Specific Origins
// ============================================================================


// ============================================================================
// Build ? at vertices
// ? := radial distance from disc centre in the plane, normalized by disc radius.
// ============================================================================

__global__
void k_buildLambda(int    N,
                   const double *x, const double *y, const double *z,
                   double cx, double cy, double cz,
                   // outDV schedule
                   double lam_iso_outDV_center,   double lam_iso_outDV_edge,
                   double lam_aniso_outDV_center, double lam_aniso_outDV_edge,
                   // inDV schedule
                   double lam_iso_inDV_center,    double lam_iso_inDV_edge,
                   double lam_aniso_inDV_center,  double lam_aniso_inDV_edge,
                   // geometry
                   double disc_radius,
                   // outputs / work
                   double *rho,
                   double *lam_rr, double *lam_pp, double *lam_ss,
                   const CVec3* e_R, const CVec3* e_phi, const CVec3* e_h,
                   Mat_3x3* lam_alpha,
                   const int * /*layerFlag*/,
                   const int * DVflag)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;

    // radial coordinate (flat projection from disc centre)
    double dx = x[tid] - cx;
    double dy = y[tid] - cy;
    (void)z; (void)cz;
    double r = sqrt(dx*dx + dy*dy) + 1e-14;
    double p = (disc_radius > 0.0) ? (r / disc_radius) : 0.0;
    p = (p < 0.0) ? 0.0 : ((p > 1.0) ? 1.0 : p);
    rho[tid] = p;

    // choose schedule
    const bool inDV = (DVflag[tid] != 0);

    double lamIso = inDV ?
        (lam_iso_inDV_center    + (lam_iso_inDV_edge    - lam_iso_inDV_center)   * (p*p)) :
        (lam_iso_outDV_center   + (lam_iso_outDV_edge   - lam_iso_outDV_center)  * (p*p));

    double lamAni = inDV ?
        (lam_aniso_inDV_center  + (lam_aniso_inDV_edge  - lam_aniso_inDV_center) * (p*p)) :
        (lam_aniso_outDV_center + (lam_aniso_outDV_edge - lam_aniso_outDV_center)* (p*p));

    lam_rr[tid] = lamIso * lamAni;
    lam_pp[tid] = lamIso / lamAni;
    lam_ss[tid] = 1.0; // no through-thickness growth

    // Inside DV stripe: only radial strain; kill tangential strain
    if (inDV) {
        lam_pp[tid] = 1.0;
    }

    // assemble tensor
    Mat_3x3 L = Mat_3x3{ CVec3(0,0,0), CVec3(0,0,0), CVec3(0,0,0) };
    axpy(lam_rr[tid], outer(e_R [tid], e_R [tid]), L);
    axpy(lam_pp[tid], outer(e_phi[tid], e_phi[tid]), L);
    axpy(lam_ss[tid], outer(e_h  [tid], e_h  [tid]), L);
    lam_alpha[tid] = L;
}

// ============================================================================
// Rest-length update by full projection of ? onto the edge direction.
// Vertical edges are flagged as -1 and are skipped.
// ============================================================================

__global__
void k_edgeRestProj(int    E,
                    const int    *e2n1,  const int *e2n2,
                    const double *x,     const double *y,   const double *z,
                    const Mat_3x3 *lam_alpha,
                    double *L0, double *Lstar,
                    const int *edgeLayerFlags) // -1 ? vertical/pillar
{
    int eid = blockIdx.x*blockDim.x + threadIdx.x;
    if (eid >= E) return;

    // Do not alter vertical/pillar edges. Layer flag = -1.
    if (edgeLayerFlags[eid] == -1) return;

    int a = e2n1[eid];
    int b = e2n2[eid];

    CVec3 dX = CVec3(x[a]-x[b], y[a]-y[b], z[a]-z[b]);
    L0[eid] = norm3(dX);

    // average tensor at edge endpoints
    Mat_3x3 La = lam_alpha[a], Lb = lam_alpha[b], Lp;

    thrust::get<0>(Lp) = CVec3(
        0.5*( c<0>(thrust::get<0>(La)) + c<0>(thrust::get<0>(Lb)) ),
        0.5*( c<1>(thrust::get<0>(La)) + c<1>(thrust::get<0>(Lb)) ),
        0.5*( c<2>(thrust::get<0>(La)) + c<2>(thrust::get<0>(Lb)) ) );

    thrust::get<1>(Lp) = CVec3(
        0.5*( c<0>(thrust::get<1>(La)) + c<0>(thrust::get<1>(Lb)) ),
        0.5*( c<1>(thrust::get<1>(La)) + c<1>(thrust::get<1>(Lb)) ),
        0.5*( c<2>(thrust::get<1>(La)) + c<2>(thrust::get<1>(Lb)) ) );

    thrust::get<2>(Lp) = CVec3(
        0.5*( c<0>(thrust::get<2>(La)) + c<0>(thrust::get<2>(Lb)) ),
        0.5*( c<1>(thrust::get<2>(La)) + c<1>(thrust::get<2>(Lb)) ),
        0.5*( c<2>(thrust::get<2>(La)) + c<2>(thrust::get<2>(Lb)) ) );

    CVec3 dX_stretch = matVec(Lp, dX);
    Lstar[eid] = norm3(dX_stretch);
}

// ============================================================================
// Public wrappers
// ============================================================================

namespace StrainTensorGPU {

void buildVertexLambda(GeneralParams& gp,
                       CoordInfoVecs& coord,
                       LambdaField&   field,
                       double         tFrac /*unused but kept for API parity*/)
{
    // size guards
    if (field.lam_rr.size() != coord.nodeLocX.size())
        field.resize(coord.nodeLocX.size());
    if (gp.rho.size() != coord.nodeLocX.size())
        gp.rho.resize(coord.nodeLocX.size());

    const int N = static_cast<int>(coord.nodeLocX.size());
    dim3 grid((N + BLOCK_SZ - 1) / BLOCK_SZ);



    // --- build ? field (uses DV mask & the basis we just built) ---
    k_buildLambda<<<grid,BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.nodeLocX.data()),
        thrust::raw_pointer_cast(coord.nodeLocY.data()),
        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        gp.centerX, gp.centerY, gp.centerZ,
        // outDV
        gp.lambda_iso_center_outDV,   gp.lambda_iso_edge_outDV,
        gp.lambda_aniso_center_outDV, gp.lambda_aniso_edge_outDV,
        // inDV
        gp.lambda_iso_center_inDV,    gp.lambda_iso_edge_inDV,
        gp.lambda_aniso_center_inDV,  gp.lambda_aniso_edge_inDV,
        // geometry
        gp.disc_radius,
        // outputs
        thrust::raw_pointer_cast(gp.rho.data()),
        thrust::raw_pointer_cast(field.lam_rr.data()),
        thrust::raw_pointer_cast(field.lam_pp.data()),
        thrust::raw_pointer_cast(field.lam_ss.data()),
        thrust::raw_pointer_cast(field.e_R.data()),
        thrust::raw_pointer_cast(field.e_phi.data()),
        thrust::raw_pointer_cast(field.e_h.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()),
        thrust::raw_pointer_cast(gp.nodes_in_upperhem.data()), // not used to gate DV
        thrust::raw_pointer_cast(gp.nodes_in_DV.data()) );
    cudaDeviceSynchronize();
}

void updateEdgeRestLengths(CoordInfoVecs&  coord,
                           GeneralParams&  gp,
                           const LambdaField& field,
                           LinearSpringInfoVecs& lsInfo,
                           int /*targetLayer, unused now*/)
{
    const int E = static_cast<int>(coord.num_edges);
    dim3 grid((E + BLOCK_SZ - 1) / BLOCK_SZ);

    k_edgeRestProj<<<grid,BLOCK_SZ>>>(
        E,
        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()),
        thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
        thrust::raw_pointer_cast(coord.nodeLocX.data()),
        thrust::raw_pointer_cast(coord.nodeLocY.data()),
        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()),
        thrust::raw_pointer_cast(lsInfo.edge_initial_length.data()),
        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
        thrust::raw_pointer_cast(gp.edges_in_upperhem.data())   // here: -1 denotes vertical edges
    );
    cudaDeviceSynchronize();
}

} // namespace StrainTensorGPU
