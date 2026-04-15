// ============================================================================
// StrainTensor.cu — Rewritten to match Fuhrmann et al. (Science Advances, 2024)
//
// KEY CHANGES from previous version:
//   1. ?_hh = 1.0 always (no volume conservation in growth tensor)
//   2. ALL edges processed in updateEdgeRestLengths (no vertical skip)
//   3. ? computed as geodesic arc length on sphere, not planar XY distance
//   4. Basis vectors: e_h = position/|position| (sphere normal),
//      e_R via Gram-Schmidt of (vertex-origin) against e_h,
//      e_phi = e_h × e_R
//   5. No sign convention correction — values are absolute stretch ratios
//   6. No DV strain redirection (inDV uses its own lambda values directly)
//   7. No strainMode gating — all edges always processed
//
// Author: Navaira Sherwani, 2025 (rewrite March 2026)
// ============================================================================

#include "StrainTensor.h"
#include "SystemStructures.h"
#include "System.h"

#include <cuda_runtime.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

#define BLOCK_SZ 256

// ============================================================================
// Tuple helpers (unchanged from before — these are correct)
// ============================================================================

template<int I> __host__ __device__
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
    {
        CVec3 r = thrust::get<0>(C);
        CVec3 a = thrust::get<0>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<0>(C) = r;
    }
    {
        CVec3 r = thrust::get<1>(C);
        CVec3 a = thrust::get<1>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<1>(C) = r;
    }
    {
        CVec3 r = thrust::get<2>(C);
        CVec3 a = thrust::get<2>(A);
        thrust::get<0>(r) += s * c<0>(a);
        thrust::get<1>(r) += s * c<1>(a);
        thrust::get<2>(r) += s * c<2>(a);
        thrust::get<2>(C) = r;
    }
}

__host__ __device__
inline CVec3 matVec(const Mat_3x3& M, const CVec3& v) {
    return CVec3(
        c<0>(thrust::get<0>(M))*c<0>(v) + c<1>(thrust::get<0>(M))*c<1>(v) + c<2>(thrust::get<0>(M))*c<2>(v),
        c<0>(thrust::get<1>(M))*c<0>(v) + c<1>(thrust::get<1>(M))*c<1>(v) + c<2>(thrust::get<1>(M))*c<2>(v),
        c<0>(thrust::get<2>(M))*c<0>(v) + c<1>(thrust::get<2>(M))*c<1>(v) + c<2>(thrust::get<2>(M))*c<2>(v)
    );
}

__host__ __device__
inline double norm3(const CVec3& v) {
    return sqrt(c<0>(v)*c<0>(v) + c<1>(v)*c<1>(v) + c<2>(v)*c<2>(v)) + 1e-14;
}

__host__ __device__
inline CVec3 normalize(const CVec3& v) {
    double n = norm3(v);
    return CVec3(c<0>(v)/n, c<1>(v)/n, c<2>(v)/n);
}

__host__ __device__
inline double dot3(const CVec3& a, const CVec3& b) {
    return c<0>(a)*c<0>(b) + c<1>(a)*c<1>(b) + c<2>(a)*c<2>(b);
}


// ============================================================================
// Region classification (unchanged — this is correct)
// ============================================================================

enum RegionType { REGION_DORSAL = 0, REGION_DV = 1, REGION_VENTRAL = 2 };

static RegionType classifyNodeRegion(double y, double DV_half_width) {
    if (y < -DV_half_width) return REGION_DORSAL;
    else if (y > DV_half_width) return REGION_VENTRAL;
    else return REGION_DV;
}

struct RegionOrigins {
    double OD_x, OD_y, OD_z;
    double OV_x, OV_y, OV_z;
    double center_x, center_y, center_z;
    double DV_half_width;
    double max_pathlength_dorsal;
    double max_pathlength_ventral;
    double max_pathlength_DV;
};


// ============================================================================
// computeBasisVectorsWithDVSeparation — HOST function
//
// Matching Fuhrmann's add_basis_vectors_to_Sph() exactly:
//   - e_h = position / |position|  (outward normal on a sphere centered at origin)
//   - e_OA = (vertex - origin) / |vertex - origin|
//   - e_R = e_OA - (e_h · e_OA) × e_h, then normalized  (Gram-Schmidt)
//   - e_phi = e_h × e_R
//   - ? = geodesic_distance(vertex, origin) / max_geodesic_distance_in_region
//
// Geodesic distance on sphere: d = R × arccos( (x_A · x_O) / R² )
// ============================================================================

void StrainTensorGPU::computeBasisVectorsWithDVSeparation(
    GeneralParams& params,
    CoordInfoVecs& coords,
    double theta_DV,
    double R)
{
    int N = params.maxNodeCount;

    // --- Copy positions to host ---
    std::vector<double> hX(N), hY(N), hZ(N);
    thrust::copy(coords.nodeLocX.begin(), coords.nodeLocX.begin() + N, hX.begin());
    thrust::copy(coords.nodeLocY.begin(), coords.nodeLocY.begin() + N, hY.begin());
    thrust::copy(coords.nodeLocZ.begin(), coords.nodeLocZ.begin() + N, hZ.begin());

    // --- Compute mesh center ---
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < N; ++i) { cx += hX[i]; cy += hY[i]; cz += hZ[i]; }
    cx /= N; cy /= N; cz /= N;

    // --- DV half-width ---
    double DV_hw = R * sin(theta_DV/2);

    // --- Compute region centroids for origin selection ---
    // Fuhrmann: outDV origin = centroid of dorsal / ventral nodes
    //           inDV origin  = point on sphere with same x, y=0, z=sqrt(R²-x²)
    double d_sx = 0, d_sy = 0, d_sz = 0; int d_n = 0;
    double v_sx = 0, v_sy = 0, v_sz = 0; int v_n = 0;

    for (int i = 0; i < N; ++i) {
        double y = hY[i];
        if (y < -DV_hw) { d_sx += hX[i]; d_sy += hY[i]; d_sz += hZ[i]; d_n++; }
        else if (y > DV_hw) { v_sx += hX[i]; v_sy += hY[i]; v_sz += hZ[i]; v_n++; }
    }

    RegionOrigins org;
    org.DV_half_width = DV_hw;
    org.center_x = cx; org.center_y = cy; org.center_z = cz;

// OLD INCORRECT CODE. SEE BELOW FOR CORRECTED CODE. 
//    if (d_n > 0) { org.OD_x = d_sx/d_n; org.OD_y = d_sy/d_n; org.OD_z = d_sz/d_n; }
//    else         { org.OD_x = cx; org.OD_y = -DV_hw - R*0.25; org.OD_z = cz; }
//
//    if (v_n > 0) { org.OV_x = v_sx/v_n; org.OV_y = v_sy/v_n; org.OV_z = v_sz/v_n; }
//    else         { org.OV_x = cx; org.OV_y = DV_hw + R*0.25; org.OV_z = cz; }
    // --- REPLACE the centroid computation block ---
// OLD (wrong — centroid of all dorsal/ventral nodes):
// if (d_n > 0) { org.OD_x = d_sx/d_n; org.OD_y = d_sy/d_n; org.OD_z = d_sz/d_n; }
// if (v_n > 0) { org.OV_x = v_sx/v_n; org.OV_y = v_sy/v_n; org.OV_z = v_sz/v_n; }

// NEW — OD/OV are points on the sphere at the DV boundary:
// OD sits at Y = -DV_hw on the sphere (dorsal edge of DV strip)
// OV sits at Y = +DV_hw on the sphere (ventral edge of DV strip)
    double z_at_dv = sqrt(R * R - DV_hw * DV_hw);
    org.OD_x = cx;      // centered in X
    org.OD_y = -DV_hw;  // dorsal edge of DV boundary
    org.OD_z = z_at_dv; // on sphere surface
    
    org.OV_x = cx;
    org.OV_y = +DV_hw;
    org.OV_z = z_at_dv;
    
    // --- First pass: compute geodesic path lengths and find max per region ---
    std::vector<double> h_pathlength(N, 0.0);
    std::vector<int> h_nodesDV(N, 0);
    org.max_pathlength_dorsal = 0;
    org.max_pathlength_ventral = 0;
    org.max_pathlength_DV = 0;

    for (int i = 0; i < N; ++i) {
        double x = hX[i], y = hY[i], z = hZ[i];
        RegionType region = classifyNodeRegion(y, DV_hw);

        double ox, oy, oz;
        if (region == REGION_DV) {
            h_nodesDV[i] = 1;
            // Fuhrmann: inDV origin = point on sphere with same x, y=0
            double r_xz = sqrt(x*x + z*z);
            if (r_xz > 1e-10) {
                ox = x;  oy = 0.0;  oz = sqrt(R*R - x*x);
                // Renormalize to sphere surface
                double omag = sqrt(ox*ox + oy*oy + oz*oz);
                if (omag > 1e-10) { ox *= R/omag; oy *= R/omag; oz *= R/omag; }
            } else {
                ox = 0; oy = 0; oz = R;  // pole
            }
        } else if (region == REGION_DORSAL) {
            h_nodesDV[i] = 0;
            ox = org.OD_x; oy = org.OD_y; oz = org.OD_z;
        } else {
            h_nodesDV[i] = 0;
            ox = org.OV_x; oy = org.OV_y; oz = org.OV_z;
        }

        // Geodesic distance: d = R * arccos( (A · O) / R² )
        double dotAO = x*ox + y*oy + z*oz;
        double rA = sqrt(x*x + y*y + z*z);
        double rO = sqrt(ox*ox + oy*oy + oz*oz);
        double cosAngle = dotAO / (rA * rO + 1e-14);
        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));  // clamp for numerical safety
        double geodesic = R * acos(cosAngle);

        h_pathlength[i] = geodesic;

        if (region == REGION_DORSAL)  org.max_pathlength_dorsal  = std::max(org.max_pathlength_dorsal,  geodesic);
        if (region == REGION_VENTRAL) org.max_pathlength_ventral = std::max(org.max_pathlength_ventral, geodesic);
        if (region == REGION_DV)      org.max_pathlength_DV      = std::max(org.max_pathlength_DV,      geodesic);
    }

    // Safety floors
    if (org.max_pathlength_dorsal  < 1e-10) org.max_pathlength_dorsal  = 1.0;
    if (org.max_pathlength_ventral < 1e-10) org.max_pathlength_ventral = 1.0;
    if (org.max_pathlength_DV      < 1e-10) org.max_pathlength_DV      = 1.0;

    // --- Second pass: normalize path lengths and compute basis vectors ---
    std::vector<double> h_eR_x(N), h_eR_y(N), h_eR_z(N);
    std::vector<double> h_ePhi_x(N), h_ePhi_y(N), h_ePhi_z(N);
    std::vector<double> h_eH_x(N), h_eH_y(N), h_eH_z(N);

    int cnt_D = 0, cnt_V = 0, cnt_DV = 0;

    for (int i = 0; i < N; ++i) {
        double x = hX[i], y = hY[i], z = hZ[i];
        RegionType region = classifyNodeRegion(y, DV_hw);

        // Normalize path length to [0,1]
        if (region == REGION_DORSAL) {
            h_pathlength[i] /= org.max_pathlength_dorsal;
            cnt_D++;
        } else if (region == REGION_VENTRAL) {
            h_pathlength[i] /= org.max_pathlength_ventral;
            cnt_V++;
        } else {
            h_pathlength[i] /= org.max_pathlength_DV;
            cnt_DV++;
        }
        h_pathlength[i] = std::max(0.0, std::min(1.0, h_pathlength[i]));

        // ----------------------------------------------------------------
        // Basis vectors — matching Fuhrmann's add_basis_vectors_to_Sph()
        // ----------------------------------------------------------------

        // e_h = outward surface normal = position / |position|
        // (Fuhrmann: e_h = np.matrix(top_balls[["x","y","z"]]); normalize row-wise)
        double r_node = sqrt(x*x + y*y + z*z);
        double eh_x, eh_y, eh_z;
        if (r_node > 1e-10) {
            eh_x = x / r_node; eh_y = y / r_node; eh_z = z / r_node;
        } else {
            eh_x = 0; eh_y = 0; eh_z = 1;
        }

        // Origin for this vertex
        double ox, oy, oz;
        if (region == REGION_DV) {
            double r_xz = sqrt(x*x + z*z);
            if (r_xz > 1e-10) {
                ox = x; oy = 0.0; oz = z * R / r_xz;
                double omag = sqrt(ox*ox + oy*oy + oz*oz);
                if (omag > 1e-10) { ox *= R/omag; oy *= R/omag; oz *= R/omag; }
            } else {
                ox = 0; oy = 0; oz = R;
            }
        } else if (region == REGION_DORSAL) {
            ox = org.OD_x; oy = org.OD_y; oz = org.OD_z;
        } else {
            ox = org.OV_x; oy = org.OV_y; oz = org.OV_z;
        }

        // e_OA = (vertex - origin) / |vertex - origin|
        double eoa_x = x - ox, eoa_y = y - oy, eoa_z = z - oz;
        double eoa_mag = sqrt(eoa_x*eoa_x + eoa_y*eoa_y + eoa_z*eoa_z);
        if (eoa_mag > 1e-10) {
            eoa_x /= eoa_mag; eoa_y /= eoa_mag; eoa_z /= eoa_mag;
        } else {
            // At the origin itself — pick an arbitrary tangent direction
            if (fabs(eh_x) < 0.9) { eoa_x = 1; eoa_y = 0; eoa_z = 0; }
            else                   { eoa_x = 0; eoa_y = 1; eoa_z = 0; }
        }

        // e_R via Gram-Schmidt: project e_OA onto tangent plane (remove e_h component)
        // Fuhrmann: M*e_R = e_OA - (e_h · e_OA) * e_h; then normalize
        double dot_h_OA = eh_x*eoa_x + eh_y*eoa_y + eh_z*eoa_z;
        double eR_x = eoa_x - dot_h_OA * eh_x;
        double eR_y = eoa_y - dot_h_OA * eh_y;
        double eR_z = eoa_z - dot_h_OA * eh_z;
        double eR_mag = sqrt(eR_x*eR_x + eR_y*eR_y + eR_z*eR_z);

        if (eR_mag > 1e-10) {
            eR_x /= eR_mag; eR_y /= eR_mag; eR_z /= eR_mag;
        } else {
            // e_OA was parallel to e_h — pick arbitrary tangent vector
            if (fabs(eh_x) < 0.9) { eR_x = 1; eR_y = 0; eR_z = 0; }
            else                   { eR_x = 0; eR_y = 1; eR_z = 0; }
            double d = eR_x*eh_x + eR_y*eh_y + eR_z*eh_z;
            eR_x -= d*eh_x; eR_y -= d*eh_y; eR_z -= d*eh_z;
            eR_mag = sqrt(eR_x*eR_x + eR_y*eR_y + eR_z*eR_z);
            if (eR_mag > 1e-10) { eR_x /= eR_mag; eR_y /= eR_mag; eR_z /= eR_mag; }
        }

        // e_phi = e_h × e_R
        double ePhi_x = eh_y*eR_z - eh_z*eR_y;
        double ePhi_y = eh_z*eR_x - eh_x*eR_z;
        double ePhi_z = eh_x*eR_y - eh_y*eR_x;
        double ePhi_mag = sqrt(ePhi_x*ePhi_x + ePhi_y*ePhi_y + ePhi_z*ePhi_z);
        if (ePhi_mag > 1e-10) {
            ePhi_x /= ePhi_mag; ePhi_y /= ePhi_mag; ePhi_z /= ePhi_mag;
        }

        h_eH_x[i]   = eh_x;   h_eH_y[i]   = eh_y;   h_eH_z[i]   = eh_z;
        h_eR_x[i]   = eR_x;   h_eR_y[i]   = eR_y;   h_eR_z[i]   = eR_z;
        h_ePhi_x[i] = ePhi_x; h_ePhi_y[i] = ePhi_y; h_ePhi_z[i] = ePhi_z;
    }

    // --- Print diagnostics ---
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  FUHRMANN-STYLE BASIS VECTORS COMPUTED" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "  DV half-width: " << DV_hw << std::endl;
    std::cout << "  OD: (" << org.OD_x << ", " << org.OD_y << ", " << org.OD_z << ")" << std::endl;
    std::cout << "  OV: (" << org.OV_x << ", " << org.OV_y << ", " << org.OV_z << ")" << std::endl;
    std::cout << "  Max geodesic pathlength - Dorsal: " << org.max_pathlength_dorsal
              << ", Ventral: " << org.max_pathlength_ventral
              << ", DV: " << org.max_pathlength_DV << std::endl;
    std::cout << "  Nodes - Dorsal: " << cnt_D << ", Ventral: " << cnt_V << ", DV: " << cnt_DV << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // --- Copy to device ---
    coords.pathlength_scaled.resize(N);
    coords.e_R_x.resize(N);  coords.e_R_y.resize(N);  coords.e_R_z.resize(N);
    coords.e_phi_x.resize(N); coords.e_phi_y.resize(N); coords.e_phi_z.resize(N);
    coords.e_h_x.resize(N);  coords.e_h_y.resize(N);  coords.e_h_z.resize(N);

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

    params.nodes_in_DV.resize(N);
    thrust::copy(h_nodesDV.begin(), h_nodesDV.end(), params.nodes_in_DV.begin());
}


// ============================================================================
// CUDA kernel: build lambda field at each vertex
//
// Matching Fuhrmann:
//   ?_rr = ?_iso(?) × ?_aniso(?)
//   ?_ff = ?_iso(?) / ?_aniso(?)
//   ?_hh = 1.0   (always — no volume conservation in the growth tensor)
//
// Currently uses linear interpolation for ?(?): ? = ?_center + (?_edge - ?_center) × ?
// To match Fuhrmann's polynomial exactly, replace with polynomial evaluation.
// ============================================================================

__global__
void k_buildLambda_Fuhrmann(
    int    N,
    // Per-vertex data (read from device)
    const double* pathlength_scaled,  // ? ? [0,1], pre-computed geodesic
    const int*    nodesDV,            // 1 = inDV, 0 = outDV
    // Basis vectors (pre-computed, stored in LambdaField)
    const CVec3* e_R, const CVec3* e_phi, const CVec3* e_h,
    // Lambda schedule: outDV
    double lam_iso_outDV_center,   double lam_iso_outDV_edge,
    double lam_aniso_outDV_center, double lam_aniso_outDV_edge,
    // Lambda schedule: inDV
    double lam_iso_inDV_center,    double lam_iso_inDV_edge,
    double lam_aniso_inDV_center,  double lam_aniso_inDV_edge,
    // Outputs
    double* lam_rr, double* lam_pp, double* lam_ss,
    double* rho_out,
    Mat_3x3* lam_alpha)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;

    double p = pathlength_scaled[tid];
    rho_out[tid] = p;

    // Choose lambda schedule based on DV membership
    const bool inDV = (nodesDV[tid] != 0);

    // Linear interpolation: ?(?) = ?_center + (?_edge - ?_center) × ?
    // This is equivalent to Fuhrmann's np.poly1d with 2 coefficients:
    //   coeffs = [slope, intercept] = [(?_edge - ?_center), ?_center]
    double lamIso, lamAni;
    if (inDV) {
        lamIso = lam_iso_inDV_center  + (lam_iso_inDV_edge  - lam_iso_inDV_center)  * p;
        lamAni = lam_aniso_inDV_center + (lam_aniso_inDV_edge - lam_aniso_inDV_center) * p;
    } else {
        lamIso = lam_iso_outDV_center  + (lam_iso_outDV_edge  - lam_iso_outDV_center)  * p;
        lamAni = lam_aniso_outDV_center + (lam_aniso_outDV_edge - lam_aniso_outDV_center) * p;
    }

    // Safety: clamp away from zero to avoid division singularity
    const double LAM_FLOOR = 0.0005;
    if (fabs(lamAni) < LAM_FLOOR) lamAni = (lamAni >= 0.0) ? LAM_FLOOR : -LAM_FLOOR;
    if (fabs(lamIso) < LAM_FLOOR) lamIso = (lamIso >= 0.0) ? LAM_FLOOR : -LAM_FLOOR;

    
    if (inDV){
        // Principal stretch ratios
        lam_rr[tid] = 1.0;//lamIso * lamAni;
        lam_pp[tid] = lamIso / lamAni;
        lam_ss[tid] = 1.0;   // <--- FUHRMANN: ?_hh = 1.0, always
    }
    else{
        // Principal stretch ratios
        lam_rr[tid] = lamIso * lamAni;
        lam_pp[tid] = lamIso / lamAni;
        lam_ss[tid] = 1.0;   // <--- FUHRMANN: ?_hh = 1.0, always

    }
    
    // Assemble full tensor: ? = ?_rr (e_R?e_R) + ?_ff (e_f?e_f) + ?_hh (e_h?e_h)
    Mat_3x3 L = Mat_3x3{ CVec3(0,0,0), CVec3(0,0,0), CVec3(0,0,0) };
    axpy(lam_rr[tid], outer(e_R[tid],   e_R[tid]),   L);
    axpy(lam_pp[tid], outer(e_phi[tid], e_phi[tid]), L);
    axpy(lam_ss[tid], outer(e_h[tid],   e_h[tid]),   L);
    lam_alpha[tid] = L;
}


// ============================================================================
// CUDA kernel: project ? onto each edge to get new rest lengths
//
// Matching Fuhrmann (array_wd.py) exactly:
//   ?_avg = 0.5 × (?_a + ?_ß)
//   l* = |?_avg · d_initial|
//
// ALL edges are processed. No vertical skip. Since ?_hh = 1, vertical
// edges (aligned with e_h) naturally produce l* ˜ l_initial.
// ============================================================================

__global__
void k_edgeRestProj_Fuhrmann(
    int    E,
    const int*    e2n1,  const int* e2n2,
    const double* x,     const double* y,     const double* z,
    const Mat_3x3* lam_alpha,
    const double* L_initial,
    double* L_final)
{
    int eid = blockIdx.x * blockDim.x + threadIdx.x;
    if (eid >= E) return;

    int a = e2n1[eid];
    int b = e2n2[eid];

    // Skip invalid edges
    if (a < 0 || b < 0 || a == INT_MAX || b == INT_MAX) {
        L_final[eid] = L_initial[eid];
        return;
    }

    // Edge vector in initial configuration
    CVec3 dX = CVec3(x[a] - x[b], y[a] - y[b], z[a] - z[b]);

    // Average tensor at edge endpoints: ?_avg = 0.5 × (?_a + ?_ß)
    Mat_3x3 La = lam_alpha[a];
    Mat_3x3 Lb = lam_alpha[b];
    Mat_3x3 Lavg;

    thrust::get<0>(Lavg) = CVec3(
        0.5 * (c<0>(thrust::get<0>(La)) + c<0>(thrust::get<0>(Lb))),
        0.5 * (c<1>(thrust::get<0>(La)) + c<1>(thrust::get<0>(Lb))),
        0.5 * (c<2>(thrust::get<0>(La)) + c<2>(thrust::get<0>(Lb))));

    thrust::get<1>(Lavg) = CVec3(
        0.5 * (c<0>(thrust::get<1>(La)) + c<0>(thrust::get<1>(Lb))),
        0.5 * (c<1>(thrust::get<1>(La)) + c<1>(thrust::get<1>(Lb))),
        0.5 * (c<2>(thrust::get<1>(La)) + c<2>(thrust::get<1>(Lb))));

    thrust::get<2>(Lavg) = CVec3(
        0.5 * (c<0>(thrust::get<2>(La)) + c<0>(thrust::get<2>(Lb))),
        0.5 * (c<1>(thrust::get<2>(La)) + c<1>(thrust::get<2>(Lb))),
        0.5 * (c<2>(thrust::get<2>(La)) + c<2>(thrust::get<2>(Lb))));

    // l* = |?_avg · d|
    CVec3 stretched = matVec(Lavg, dX);
    L_final[eid] = norm3(stretched);
}


// ============================================================================
// Public wrappers
// ============================================================================

namespace StrainTensorGPU {

void buildVertexLambda(GeneralParams& gp,
                       CoordInfoVecs& coord,
                       LambdaField&   field,
                       double         tFrac)
{
    const int N = static_cast<int>(coord.nodeLocX.size());

    // Ensure field is sized
    if (field.lam_rr.size() != (size_t)N) field.resize(N);
    if (gp.rho.size() != (size_t)N) gp.rho.resize(N);

    dim3 grid((N + BLOCK_SZ - 1) / BLOCK_SZ);

    // Print lambda values being used (first call only)
    static bool printed = false;
    if (!printed) {
        std::cout << "\n  buildVertexLambda (Fuhrmann-style):" << std::endl;
        std::cout << "    outDV: iso=[" << gp.lambda_iso_center_outDV << " ? " << gp.lambda_iso_edge_outDV
                  << "], aniso=[" << gp.lambda_aniso_center_outDV << " ? " << gp.lambda_aniso_edge_outDV << "]" << std::endl;
        std::cout << "    inDV:  iso=[" << gp.lambda_iso_center_inDV << " ? " << gp.lambda_iso_edge_inDV
                  << "], aniso=[" << gp.lambda_aniso_center_inDV << " ? " << gp.lambda_aniso_edge_inDV << "]" << std::endl;
        std::cout << "    ?_hh = 1.0 (no volume conservation in growth tensor)" << std::endl;
        printed = true;
    }

    // Need pathlength_scaled on device — it was stored in coordInfoVecs by computeBasisVectors
    // The basis vectors are already in field.e_R, field.e_phi, field.e_h
    // (copied from coordInfoVecs in System.cu)

    k_buildLambda_Fuhrmann<<<grid, BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.pathlength_scaled.data()),
        thrust::raw_pointer_cast(gp.nodes_in_DV.data()),
        thrust::raw_pointer_cast(field.e_R.data()),
        thrust::raw_pointer_cast(field.e_phi.data()),
        thrust::raw_pointer_cast(field.e_h.data()),
        // outDV
        gp.lambda_iso_center_outDV,   gp.lambda_iso_edge_outDV,
        gp.lambda_aniso_center_outDV, gp.lambda_aniso_edge_outDV,
        // inDV
        gp.lambda_iso_center_inDV,    gp.lambda_iso_edge_inDV,
        gp.lambda_aniso_center_inDV,  gp.lambda_aniso_edge_inDV,
        // outputs
        thrust::raw_pointer_cast(field.lam_rr.data()),
        thrust::raw_pointer_cast(field.lam_pp.data()),
        thrust::raw_pointer_cast(field.lam_ss.data()),
        thrust::raw_pointer_cast(gp.rho.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()));

    cudaDeviceSynchronize();
}


void updateEdgeRestLengths(CoordInfoVecs&        coord,
                           GeneralParams&        gp,
                           const LambdaField&    field,
                           LinearSpringInfoVecs& lsInfo,
                           int                   /*targetLayer*/)
{
    const int E = static_cast<int>(coord.num_edges);
    dim3 grid((E + BLOCK_SZ - 1) / BLOCK_SZ);

    // No strainMode gating — ALL edges get processed, matching Fuhrmann
    k_edgeRestProj_Fuhrmann<<<grid, BLOCK_SZ>>>(
        E,
        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()),
        thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
        thrust::raw_pointer_cast(coord.nodeLocX.data()),
        thrust::raw_pointer_cast(coord.nodeLocY.data()),
        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()),
        thrust::raw_pointer_cast(lsInfo.edge_initial_length.data()),
        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()));

    cudaDeviceSynchronize();

    // Print diagnostics (first call only)
    static bool diag_printed = false;
    if (!diag_printed) {
        thrust::host_vector<double> h_init = lsInfo.edge_initial_length;
        thrust::host_vector<double> h_final = lsInfo.edge_final_length;

        double min_ratio = 1e10, max_ratio = -1e10, sum_ratio = 0;
        int valid = 0;
        for (int e = 0; e < E; ++e) {
            if (h_init[e] > 1e-10) {
                double ratio = h_final[e] / h_init[e];
                min_ratio = std::min(min_ratio, ratio);
                max_ratio = std::max(max_ratio, ratio);
                sum_ratio += ratio;
                valid++;
            }
        }
        if (valid > 0) {
            std::cout << "\n  Edge rest-length update (Fuhrmann-style, ALL edges):" << std::endl;
            std::cout << "    l*/l0 - min: " << min_ratio << ", max: " << max_ratio
                      << ", avg: " << sum_ratio / valid << std::endl;
            std::cout << "    Total edges: " << E << ", valid: " << valid << std::endl;
        }
        diag_printed = true;
    }
}

// Unused stub
void updatePreferredAngles(BendingTriangleInfoVecs&,
                           const CoordInfoVecs&,
                           const LambdaField&,
                           const GeneralParams&,
                           const LinearSpringInfoVecs&) {}

} // namespace StrainTensorGPU




// Previously working Strain field that has been commented out for now. Nav 03/12/26






//// StrainTensor.cu  -- DV-aware basis construction and rest-length projection
////
//// Spontaneous strain tensor Lambda = lambda_rr (e_R x e_R) + lambda_ff (e_phi x e_phi) + lambda_hh (e_h x e_h)
//// with lambda_rr = lambda_iso * lambda_aniso,  lambda_ff = lambda_iso / lambda_aniso,  lambda_hh computed from volume conservation.
////
//// SIGN CONVENTION FIX (v2):
////   Lambda values from the XML can be either:
////     (a) Absolute stretch ratios  (all positive, typically 0.5-2.0)
////     (b) Strain deltas            (can be negative, represent Deltalambda = lambda_target - 1)
////
////   If any center lambda is negative, we auto-detect mode (b) and convert:
////     lambda_actual = 1 + lambda_xml
////
////   This prevents the division-by-zero singularity that occurs when
////   lambda_aniso interpolates through zero between negative center and positive edge values.
////
//// Inside the DV stripe:
////   - e_R is the *radial* direction away from the line segment joining O_D and O_V
////     (the DV boundary line).
////   - Only along-stripe strain is applied: lambda_pp = lambda_rr (computed), lambda_rr = 1, lambda_ss via vol conservation.
//// Outside the DV stripe:
////   - e_R is constructed as before using OA from OD (dorsal) or OV (ventral).
////
//// Vertical/pillar edges are flagged as -1 and are not altered in the rest-length update.
////
//// Region-aware strain gating (strainMode):
////   0 = apply strain to ALL regions (both inDV and outDV)
////   1 = apply strain to inDV edges only  (outDV edges untouched)
////   2 = apply strain to outDV edges only (inDV edges untouched)
//// The mode is determined automatically from lambda values on the host side.
//
//#include "StrainTensor.h"
//#include "SystemStructures.h"
//#include "System.h"
//
//#include <cuda_runtime.h>
//#include <thrust/tuple.h>
//#include <thrust/device_vector.h>
//#include <thrust/iterator/counting_iterator.h>
//#include <thrust/sequence.h>
//#include <cmath>
//
//#define BLOCK_SZ 256
//
//// ------------ tuple helpers -----------------
//template<int I>  __host__ __device__
//inline double c(const CVec3& v) { return thrust::get<I>(v); }
//
//__host__ __device__
//inline Mat_3x3 outer(const CVec3& a, const CVec3& b) {
//    return Mat_3x3(
//        CVec3( c<0>(a)*c<0>(b), c<0>(a)*c<1>(b), c<0>(a)*c<2>(b) ),
//        CVec3( c<1>(a)*c<0>(b), c<1>(a)*c<1>(b), c<1>(a)*c<2>(b) ),
//        CVec3( c<2>(a)*c<0>(b), c<2>(a)*c<1>(b), c<2>(a)*c<2>(b) )
//    );
//}
//
//__host__ __device__
//inline void axpy(double s, const Mat_3x3& A, Mat_3x3& C) {
//    {
//        CVec3 r = thrust::get<0>(C);
//        CVec3 a = thrust::get<0>(A);
//        thrust::get<0>(r) += s * c<0>(a);
//        thrust::get<1>(r) += s * c<1>(a);
//        thrust::get<2>(r) += s * c<2>(a);
//        thrust::get<0>(C)   = r;
//    }
//    {
//        CVec3 r = thrust::get<1>(C);
//        CVec3 a = thrust::get<1>(A);
//        thrust::get<0>(r) += s * c<0>(a);
//        thrust::get<1>(r) += s * c<1>(a);
//        thrust::get<2>(r) += s * c<2>(a);
//        thrust::get<1>(C)   = r;
//    }
//    {
//        CVec3 r = thrust::get<2>(C);
//        CVec3 a = thrust::get<2>(A);
//        thrust::get<0>(r) += s * c<0>(a);
//        thrust::get<1>(r) += s * c<1>(a);
//        thrust::get<2>(r) += s * c<2>(a);
//        thrust::get<2>(C)   = r;
//    }
//}
//
//__host__ __device__
//inline CVec3 matVec(const Mat_3x3& M, const CVec3& v) {
//    return CVec3(
//        c<0>( thrust::get<0>(M) )*c<0>(v) + c<1>( thrust::get<0>(M) )*c<1>(v) + c<2>( thrust::get<0>(M) )*c<2>(v),
//        c<0>( thrust::get<1>(M) )*c<0>(v) + c<1>( thrust::get<1>(M) )*c<1>(v) + c<2>( thrust::get<1>(M) )*c<2>(v),
//        c<0>( thrust::get<2>(M) )*c<0>(v) + c<1>( thrust::get<2>(M) )*c<1>(v) + c<2>( thrust::get<2>(M) )*c<2>(v)
//    );
//}
//
//__host__ __device__
//inline double norm3(const CVec3& v) {
//    return sqrt(thrust::get<0>(v)*thrust::get<0>(v) +
//                thrust::get<1>(v)*thrust::get<1>(v) +
//                thrust::get<2>(v)*thrust::get<2>(v)) + 1e-14;
//}
//
//__host__ __device__
//inline CVec3 normalize(const CVec3& v) {
//    double n = norm3(v);
//    return CVec3( c<0>(v)/n, c<1>(v)/n, c<2>(v)/n );
//}
//
//__host__ __device__
//inline CVec3 cross(const CVec3& a, const CVec3& b) {
//    return CVec3(
//        c<1>(a)*c<2>(b) - c<2>(a)*c<1>(b),
//        c<2>(a)*c<0>(b) - c<0>(a)*c<2>(b),
//        c<0>(a)*c<1>(b) - c<1>(a)*c<0>(b)
//    );
//}
//
//__host__ __device__
//inline double dot3(const CVec3& a, const CVec3& b) {
//    return c<0>(a)*c<0>(b) + c<1>(a)*c<1>(b) + c<2>(a)*c<2>(b);
//}
//
//
//// ============================================================================
//// Region Classification
//// ============================================================================
//
//enum RegionType {
//    REGION_DORSAL = 0,
//    REGION_DV = 1,
//    REGION_VENTRAL = 2
//};
//
//static RegionType classifyNodeRegion(double y, double DV_half_width) {
//    if (y < -DV_half_width) {
//        return REGION_DORSAL;
//    } else if (y > DV_half_width) {
//        return REGION_VENTRAL;
//    } else {
//        return REGION_DV;
//    }
//}
//
//// ============================================================================
//// Structure to hold region-specific origins
//// ============================================================================
//
//struct RegionOrigins {
//    double OD_x, OD_y, OD_z;
//    double OV_x, OV_y, OV_z;
//    double center_x, center_y, center_z;
//    double DV_half_width;
//    double max_r_dorsal;
//    double max_r_ventral;
//    double max_rho_DV;
//};
//
//// ============================================================================
//// Compute Basis Vectors with Region-Specific Origins
//// Works entirely with device vectors - no HostSetInfoVecs
//// ============================================================================
//
//void StrainTensorGPU::computeBasisVectorsWithDVSeparation(
//    GeneralParams& params,
//    CoordInfoVecs& coords,
//    double theta_DV,
//    double R)
//{
//    int N = params.maxNodeCount;
//    
//    // ========================================================================
//    // STEP 1: Copy device vectors to local std::vectors for processing
//    // ========================================================================
//    
//    std::vector<double> h_nodeLocX(N);
//    std::vector<double> h_nodeLocY(N);
//    std::vector<double> h_nodeLocZ(N);
//    
//    thrust::copy(coords.nodeLocX.begin(), coords.nodeLocX.begin() + N, h_nodeLocX.begin());
//    thrust::copy(coords.nodeLocY.begin(), coords.nodeLocY.begin() + N, h_nodeLocY.begin());
//    thrust::copy(coords.nodeLocZ.begin(), coords.nodeLocZ.begin() + N, h_nodeLocZ.begin());
//    
//    // ========================================================================
//    // STEP 2: Compute region origins
//    // ========================================================================
//    
//    RegionOrigins origins;
//    origins.DV_half_width = R * sin(theta_DV);
//    
//    // Compute mesh center
//    double sum_x = 0, sum_y = 0, sum_z = 0;
//    for (int i = 0; i < N; ++i) {
//        sum_x += h_nodeLocX[i];
//        sum_y += h_nodeLocY[i];
//        sum_z += h_nodeLocZ[i];
//    }
//    origins.center_x = sum_x / N;
//    origins.center_y = sum_y / N;
//    origins.center_z = sum_z / N;
//    
//    // Compute dorsal and ventral centroids
//    double dorsal_sum_x = 0, dorsal_sum_y = 0, dorsal_sum_z = 0;
//    double ventral_sum_x = 0, ventral_sum_y = 0, ventral_sum_z = 0;
//    int dorsal_count = 0, ventral_count = 0;
//    
//    for (int i = 0; i < N; ++i) {
//        double y = h_nodeLocY[i];
//        
//        if (y < -origins.DV_half_width) {
//            dorsal_sum_x += h_nodeLocX[i];
//            dorsal_sum_y += h_nodeLocY[i];
//            dorsal_sum_z += h_nodeLocZ[i];
//            dorsal_count++;
//        } else if (y > origins.DV_half_width) {
//            ventral_sum_x += h_nodeLocX[i];
//            ventral_sum_y += h_nodeLocY[i];
//            ventral_sum_z += h_nodeLocZ[i];
//            ventral_count++;
//        }
//    }
//    
//    if (dorsal_count > 0) {
//        origins.OD_x = dorsal_sum_x / dorsal_count;
//        origins.OD_y = dorsal_sum_y / dorsal_count;
//        origins.OD_z = dorsal_sum_z / dorsal_count;
//    } else {
//        origins.OD_x = origins.center_x;
//        origins.OD_y = -origins.DV_half_width - R * 0.25;
//        origins.OD_z = origins.center_z;
//    }
//    
//    if (ventral_count > 0) {
//        origins.OV_x = ventral_sum_x / ventral_count;
//        origins.OV_y = ventral_sum_y / ventral_count;
//        origins.OV_z = ventral_sum_z / ventral_count;
//    } else {
//        origins.OV_x = origins.center_x;
//        origins.OV_y = origins.DV_half_width + R * 0.25;
//        origins.OV_z = origins.center_z;
//    }
//    
//    // Compute maximum radial distances for normalization
//    origins.max_r_dorsal = 0;
//    origins.max_r_ventral = 0;
//    origins.max_rho_DV = 0;
//    
//    for (int i = 0; i < N; ++i) {
//        double x = h_nodeLocX[i];
//        double y = h_nodeLocY[i];
//        double z = h_nodeLocZ[i];
//        
//        RegionType region = classifyNodeRegion(y, origins.DV_half_width);
//        
//        if (region == REGION_DORSAL) {
//            double dx = x - origins.OD_x;
//            double dy = y - origins.OD_y;
//            double dz = z - origins.OD_z;
//            double r = sqrt(dx*dx + dy*dy + dz*dz);
//            origins.max_r_dorsal = std::max(origins.max_r_dorsal, r);
//        } else if (region == REGION_VENTRAL) {
//            double dx = x - origins.OV_x;
//            double dy = y - origins.OV_y;
//            double dz = z - origins.OV_z;
//            double r = sqrt(dx*dx + dy*dy + dz*dz);
//            origins.max_r_ventral = std::max(origins.max_r_ventral, r);
//        } else {
//            double rho = fabs(y);
//            origins.max_rho_DV = std::max(origins.max_rho_DV, rho);
//        }
//    }
//    
//    if (origins.max_r_dorsal < 1e-10) origins.max_r_dorsal = R;
//    if (origins.max_r_ventral < 1e-10) origins.max_r_ventral = R;
//    if (origins.max_rho_DV < 1e-10) origins.max_rho_DV = origins.DV_half_width;
//    
//    std::cout << "\n" << std::string(60, '=') << std::endl;
//    std::cout << "DV BOUNDARY SEPARATION - Region Origins" << std::endl;
//    std::cout << std::string(60, '=') << std::endl;
//    std::cout << "DV half-width: " << origins.DV_half_width << std::endl;
//    std::cout << "Mesh center: (" << origins.center_x << ", " 
//              << origins.center_y << ", " << origins.center_z << ")" << std::endl;
//    std::cout << "OD (Dorsal origin): (" << origins.OD_x << ", " 
//              << origins.OD_y << ", " << origins.OD_z << ")" << std::endl;
//    std::cout << "OV (Ventral origin): (" << origins.OV_x << ", " 
//              << origins.OV_y << ", " << origins.OV_z << ")" << std::endl;
//    std::cout << "Max radial distances:" << std::endl;
//    std::cout << "  Dorsal: " << origins.max_r_dorsal << std::endl;
//    std::cout << "  Ventral: " << origins.max_r_ventral << std::endl;
//    std::cout << "  DV (rho): " << origins.max_rho_DV << std::endl;
//    std::cout << std::string(60, '=') << std::endl;
//    
//    // ========================================================================
//    // STEP 3: Compute basis vectors and store in local vectors
//    // ========================================================================
//    
//    std::vector<double> h_pathlength(N);
//    std::vector<double> h_eR_x(N), h_eR_y(N), h_eR_z(N);
//    std::vector<double> h_ePhi_x(N), h_ePhi_y(N), h_ePhi_z(N);
//    std::vector<double> h_eH_x(N), h_eH_y(N), h_eH_z(N);
//    std::vector<int> h_nodes_in_DV(N);
//    
//    int count_dorsal = 0, count_ventral = 0, count_DV = 0;
//    
//    for (int i = 0; i < N; ++i) {
//        double x = h_nodeLocX[i];
//        double y = h_nodeLocY[i];
//        double z = h_nodeLocZ[i];
//        
//        RegionType region = classifyNodeRegion(y, origins.DV_half_width);
//        
//        double radial_x, radial_y, radial_z;
//        double pathlength;
//        
//        if (region == REGION_DV) {
//            h_nodes_in_DV[i] = 1;
//            count_DV++;
//            
//            radial_x = 0;
//            radial_y = (y >= 0) ? 1.0 : -1.0;
//            radial_z = 0;
//            
//            double rho = fabs(y);
//            pathlength = rho / origins.max_rho_DV;
//            
//        } else if (region == REGION_DORSAL) {
//            h_nodes_in_DV[i] = 0;
//            count_dorsal++;
//            
//            double dx = x - origins.OD_x;
//            double dy = y - origins.OD_y;
//            double dz = z - origins.OD_z;
//            double r = sqrt(dx*dx + dy*dy + dz*dz);
//            
//            if (r > 1e-10) {
//                radial_x = dx / r;
//                radial_y = dy / r;
//                radial_z = dz / r;
//            } else {
//                radial_x = 0;
//                radial_y = -1;
//                radial_z = 0;
//            }
//            
//            pathlength = r / origins.max_r_dorsal;
//            
//        } else {
//            h_nodes_in_DV[i] = 0;
//            count_ventral++;
//            
//            double dx = x - origins.OV_x;
//            double dy = y - origins.OV_y;
//            double dz = z - origins.OV_z;
//            double r = sqrt(dx*dx + dy*dy + dz*dz);
//            
//            if (r > 1e-10) {
//                radial_x = dx / r;
//                radial_y = dy / r;
//                radial_z = dz / r;
//            } else {
//                radial_x = 0;
//                radial_y = 1;
//                radial_z = 0;
//            }
//            
//            pathlength = r / origins.max_r_ventral;
//        }
//        
//        h_eR_x[i] = radial_x;
//        h_eR_y[i] = radial_y;
//        h_eR_z[i] = radial_z;
//        h_pathlength[i] = pathlength;
//        
//        // Compute raw surface normal
//        double nx = x - origins.center_x;
//        double ny = y - origins.center_y;
//        double nz = z - origins.center_z;
//        double n_mag = sqrt(nx*nx + ny*ny + nz*nz);
//        
//        if (n_mag > 1e-10) {
//            nx /= n_mag; ny /= n_mag; nz /= n_mag;
//        } else {
//            nx = 0; ny = 0; nz = 1;
//        }
//        
//        // Gram-Schmidt: orthogonalize e_h against e_R
//        double dot_h_R = nx*radial_x + ny*radial_y + nz*radial_z;
//        
//        double hx = nx - dot_h_R * radial_x;
//        double hy = ny - dot_h_R * radial_y;
//        double hz = nz - dot_h_R * radial_z;
//        
//        double h_mag = sqrt(hx*hx + hy*hy + hz*hz);
//        
//        if (h_mag > 1e-10) {
//            hx /= h_mag; hy /= h_mag; hz /= h_mag;
//        } else {
//            if (fabs(radial_x) < 0.9) {
//                hx = 1; hy = 0; hz = 0;
//            } else {
//                hx = 0; hy = 1; hz = 0;
//            }
//            double dot = hx*radial_x + hy*radial_y + hz*radial_z;
//            hx -= dot * radial_x;
//            hy -= dot * radial_y;
//            hz -= dot * radial_z;
//            h_mag = sqrt(hx*hx + hy*hy + hz*hz);
//            if (h_mag > 1e-10) {
//                hx /= h_mag; hy /= h_mag; hz /= h_mag;
//            }
//        }
//        
//        h_eH_x[i] = hx;
//        h_eH_y[i] = hy;
//        h_eH_z[i] = hz;
//        
//        // e_phi = e_h x e_R
//        double phi_x = hy * radial_z - hz * radial_y;
//        double phi_y = hz * radial_x - hx * radial_z;
//        double phi_z = hx * radial_y - hy * radial_x;
//        
//        double phi_mag = sqrt(phi_x*phi_x + phi_y*phi_y + phi_z*phi_z);
//        
//        if (phi_mag > 1e-10) {
//            phi_x /= phi_mag; phi_y /= phi_mag; phi_z /= phi_mag;
//        } else {
//            phi_x = 0; phi_y = 0; phi_z = 1;
//        }
//        
//        h_ePhi_x[i] = phi_x;
//        h_ePhi_y[i] = phi_y;
//        h_ePhi_z[i] = phi_z;
//    }
//    
//    std::cout << "\nNode Region Distribution:" << std::endl;
//    std::cout << "  Dorsal nodes: " << count_dorsal << std::endl;
//    std::cout << "  DV boundary nodes: " << count_DV << std::endl;
//    std::cout << "  Ventral nodes: " << count_ventral << std::endl;
//    std::cout << "  Total: " << (count_dorsal + count_DV + count_ventral) << std::endl;
//    
//    // ========================================================================
//    // STEP 4: Copy computed values to device vectors
//    // ========================================================================
//    
//    coords.pathlength_scaled.resize(N);
//    coords.e_R_x.resize(N);  coords.e_R_y.resize(N);  coords.e_R_z.resize(N);
//    coords.e_phi_x.resize(N); coords.e_phi_y.resize(N); coords.e_phi_z.resize(N);
//    coords.e_h_x.resize(N);  coords.e_h_y.resize(N);  coords.e_h_z.resize(N);
//    
//    thrust::copy(h_pathlength.begin(), h_pathlength.end(), coords.pathlength_scaled.begin());
//    thrust::copy(h_eR_x.begin(), h_eR_x.end(), coords.e_R_x.begin());
//    thrust::copy(h_eR_y.begin(), h_eR_y.end(), coords.e_R_y.begin());
//    thrust::copy(h_eR_z.begin(), h_eR_z.end(), coords.e_R_z.begin());
//    thrust::copy(h_ePhi_x.begin(), h_ePhi_x.end(), coords.e_phi_x.begin());
//    thrust::copy(h_ePhi_y.begin(), h_ePhi_y.end(), coords.e_phi_y.begin());
//    thrust::copy(h_ePhi_z.begin(), h_ePhi_z.end(), coords.e_phi_z.begin());
//    thrust::copy(h_eH_x.begin(), h_eH_x.end(), coords.e_h_x.begin());
//    thrust::copy(h_eH_y.begin(), h_eH_y.end(), coords.e_h_y.begin());
//    thrust::copy(h_eH_z.begin(), h_eH_z.end(), coords.e_h_z.begin());
//    
//    params.nodes_in_DV.resize(N);
//    thrust::copy(h_nodes_in_DV.begin(), h_nodes_in_DV.end(), params.nodes_in_DV.begin());
//    
//    std::cout << "Basis vectors computed and copied to device. N = " << N << std::endl;
//}
//
//
//// ============================================================================
//// Build lambda at vertices
//// rho := radial distance from disc centre in the plane, normalized by disc radius.
////
//// KEY FIX: Lambda values are now interpreted correctly regardless of whether
//// the XML provides absolute stretch ratios or strain deltas.
//// Volume conservation: lambda_ss = 1 / (lambda_rr * lambda_pp) for full 3D incompressibility.
//// ============================================================================
//
//__global__
//void k_buildLambda(int    N,
//                   const double *x, const double *y, const double *z,
//                   double cx, double cy, double cz,
//                   // outDV schedule (ALREADY CORRECTED on host side)
//                   double lam_iso_outDV_center,   double lam_iso_outDV_edge,
//                   double lam_aniso_outDV_center, double lam_aniso_outDV_edge,
//                   // inDV schedule  (ALREADY CORRECTED on host side)
//                   double lam_iso_inDV_center,    double lam_iso_inDV_edge,
//                   double lam_aniso_inDV_center,  double lam_aniso_inDV_edge,
//                   // geometry
//                   double disc_radius,
//                   // outputs / work
//                   double *rho,
//                   double *lam_rr, double *lam_pp, double *lam_ss,
//                   const CVec3* e_R, const CVec3* e_phi, const CVec3* e_h,
//                   Mat_3x3* lam_alpha,
//                   const int * /*layerFlag*/,
//                   const int * DVflag)
//{
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (tid >= N) return;
//
//    // radial coordinate (flat projection from disc centre)
//    double dx = x[tid] - cx;
//    double dy = y[tid] - cy;
//    (void)z; (void)cz;
//    double r = sqrt(dx*dx + dy*dy) + 1e-14;
//    double p = (disc_radius > 0.0) ? (r / disc_radius) : 0.0;
//    p = (p < 0.0) ? 0.0 : ((p > 1.0) ? 1.0 : p);
//    rho[tid] = p;
//
//    // choose schedule based on DV membership
//    const bool inDV = (DVflag[tid] != 0);
//
//    double lamIso = inDV ?
//        (lam_iso_inDV_center    + (lam_iso_inDV_edge    - lam_iso_inDV_center)   * (p*p)) :
//        (lam_iso_outDV_center   + (lam_iso_outDV_edge   - lam_iso_outDV_center)  * (p*p));
//
//    double lamAni = inDV ?
//        (lam_aniso_inDV_center  + (lam_aniso_inDV_edge  - lam_aniso_inDV_center) * (p*p)) :
//        (lam_aniso_outDV_center + (lam_aniso_outDV_edge - lam_aniso_outDV_center)* (p*p));
//
//    // ========================================================================
//    // SINGULARITY PROTECTION: Clamp lambda_aniso away from zero
//    // Division by zero occurs if lamAni passes through zero during interpolation.
//    // After the host-side delta correction, this should not happen, but we
//    // keep this as a safety net.
//    // ========================================================================
//    const double LAM_FLOOR = 0.05;  // minimum allowed stretch ratio
//    if (fabs(lamAni) < LAM_FLOOR) {
//        lamAni = (lamAni >= 0.0) ? LAM_FLOOR : -LAM_FLOOR;
//    }
//    if (fabs(lamIso) < LAM_FLOOR) {
//        lamIso = (lamIso >= 0.0) ? LAM_FLOOR : -LAM_FLOOR;
//    }
//
//    // Compute principal stretch ratios
//    lam_rr[tid] = lamIso * lamAni;
//    lam_pp[tid] = lamIso / lamAni;
//
//    // ========================================================================
//    // VOLUME CONSERVATION: lambda_ss = 1 / (lambda_rr * lambda_pp)
//    //
//    // For an incompressible tissue: lambda_rr * lambda_pp * lambda_ss = 1
//    // This drives tissue thinning when in-plane stretches expand the tissue,
//    // which is the key missing physics for the wL3 -> 4hAPF transition.
//    //
//    // Fuhrmann et al. (Science Advances, 2024) use a thin-shell model where
//    // height changes arise from in-plane strain incompatibility. In our 3D
//    // model, explicit volume conservation gives the equivalent effect.
//    // ========================================================================
//    double in_plane_product = fabs(lam_rr[tid] * lam_pp[tid]);
//    if (in_plane_product < 0.01) in_plane_product = 0.01;  // safety floor
//    lam_ss[tid] = 1.0 ;/// in_plane_product;
//
//    // Inside DV stripe: redirect strain to along-stripe direction
//    if (inDV) {
//        double along_stripe_strain = lam_rr[tid];
//        lam_rr[tid] = along_stripe_strain;//1.0;                    // kill across-stripe (e_R = y^)
//        lam_pp[tid] = 1.0;     // apply along-stripe (e_phi ~ x^)
//        
//        // Recompute volume conservation from ACTUAL applied stretches
//       // in_plane_product = fabs(lam_rr[tid] * lam_pp[tid]);
//        //if (in_plane_product < 0.01) in_plane_product = 0.01;
//        //lam_ss[tid] = 1.0 / in_plane_product;
//    }
//    
//    // Safety clamp: tissue can't change height by >10x
//    //if (lam_ss[tid] > 10.0) lam_ss[tid] = 10.0;
//    //if (lam_ss[tid] < 0.1)  lam_ss[tid] = 0.1;
//
//    // assemble tensor Lambda = lambda_rr (e_Rxe_R) + lambda_pp (e_phixe_phi) + lambda_ss (e_hxe_h)
//    Mat_3x3 L = Mat_3x3{ CVec3(0,0,0), CVec3(0,0,0), CVec3(0,0,0) };
//    axpy(lam_rr[tid], outer(e_R [tid], e_R [tid]), L);
//    axpy(lam_pp[tid], outer(e_phi[tid], e_phi[tid]), L);
//    axpy(lam_ss[tid], outer(e_h  [tid], e_h  [tid]), L);
//    lam_alpha[tid] = L;
//}
//
//// ============================================================================
//// Rest-length update by full projection of Lambda onto the edge direction.
//// Vertical edges are flagged as -1 and are skipped.
////
//// strainMode controls which region's edges get processed:
////   0 = ALL edges (both inDV and outDV get strain)
////   1 = inDV ONLY  (outDV edges are left untouched)
////   2 = outDV ONLY (inDV edges are left untouched)
//// ============================================================================
//
//__global__
//void k_edgeRestProj(int    E,
//                    const int    *e2n1,  const int *e2n2,
//                    const double *x,     const double *y,   const double *z,
//                    const Mat_3x3 *lam_alpha,
//                    double *L0, double *Lstar,
//                    const int *edgeLayerFlags,
//                    const int *nodesDVflag,
//                    int    strainMode)
//{
//    int eid = blockIdx.x*blockDim.x + threadIdx.x;
//    if (eid >= E) return;
//
//    // Do not alter vertical/pillar edges
//    if (edgeLayerFlags[eid] == -1) return;
//
//    int a = e2n1[eid];
//    int b = e2n2[eid];
//
//    const bool a_inDV = (nodesDVflag[a] != 0);
//    const bool b_inDV = (nodesDVflag[b] != 0);
//
//    if (strainMode == 1) {
//        if (!a_inDV && !b_inDV) return;
//    } else if (strainMode == 2) {
//        if (a_inDV && b_inDV) return;
//    }
//
//    CVec3 dX = CVec3(x[a]-x[b], y[a]-y[b], z[a]-z[b]);
//    L0[eid] = norm3(dX);
//
//    // average tensor at edge endpoints
//    Mat_3x3 La = lam_alpha[a], Lb = lam_alpha[b], Lp;
//
//    thrust::get<0>(Lp) = CVec3(
//        0.5*( c<0>(thrust::get<0>(La)) + c<0>(thrust::get<0>(Lb)) ),
//        0.5*( c<1>(thrust::get<0>(La)) + c<1>(thrust::get<0>(Lb)) ),
//        0.5*( c<2>(thrust::get<0>(La)) + c<2>(thrust::get<0>(Lb)) ) );
//
//    thrust::get<1>(Lp) = CVec3(
//        0.5*( c<0>(thrust::get<1>(La)) + c<0>(thrust::get<1>(Lb)) ),
//        0.5*( c<1>(thrust::get<1>(La)) + c<1>(thrust::get<1>(Lb)) ),
//        0.5*( c<2>(thrust::get<1>(La)) + c<2>(thrust::get<1>(Lb)) ) );
//
//    thrust::get<2>(Lp) = CVec3(
//        0.5*( c<0>(thrust::get<2>(La)) + c<0>(thrust::get<2>(Lb)) ),
//        0.5*( c<1>(thrust::get<2>(La)) + c<1>(thrust::get<2>(Lb)) ),
//        0.5*( c<2>(thrust::get<2>(La)) + c<2>(thrust::get<2>(Lb)) ) );
//
//    CVec3 dX_stretch = matVec(Lp, dX);
//    Lstar[eid] = norm3(dX_stretch);
//}
//
//
//// ============================================================================
//// Public wrappers
//// ============================================================================
//
//namespace StrainTensorGPU {
//
//// ============================================================================
//// Helper: Auto-detect and correct lambda sign convention
////
//// If any center lambda is negative, we interpret ALL lambdas as deltas:
////   lambda_actual = 1.0 + lambda_xml
////
//// This converts:
////   -0.124  -> 0.876  (12.4% contraction)
////   -0.062  -> 0.938  (6.2% contraction)
////    1.208  -> 2.208  (120.8% expansion)  -- but this seems too large
////
//// ALTERNATIVE: If center values are negative but edge values are > 1,
//// it's more likely that center values are deltas and edge values are
//// absolute. We handle this by checking both.
////
//// Returns true if correction was applied.
//// ============================================================================
//static bool correctLambdaSignConvention(GeneralParams& gp)
//{
//    bool any_negative =
//        (gp.lambda_iso_center_outDV < 0.0) ||
//        (gp.lambda_aniso_center_outDV < 0.0) ||
//        (gp.lambda_iso_center_inDV < 0.0) ||
//        (gp.lambda_aniso_center_inDV < 0.0);
//
//    if (!any_negative) {
//        // All positive -- check if they're reasonable stretch ratios
//        bool all_reasonable =
//            (gp.lambda_iso_center_outDV > 0.05 && gp.lambda_iso_center_outDV < 20.0) &&
//            (gp.lambda_aniso_center_outDV > 0.05 && gp.lambda_aniso_center_outDV < 20.0) &&
//            (gp.lambda_iso_edge_outDV > 0.05 && gp.lambda_iso_edge_outDV < 20.0) &&
//            (gp.lambda_aniso_edge_outDV > 0.05 && gp.lambda_aniso_edge_outDV < 20.0);
//
//        if (all_reasonable) {
//            std::cout << "  Lambda sign convention: ABSOLUTE stretch ratios (all positive, no correction needed)" << std::endl;
//            return false;
//        }
//    }
//
//    // ========================================================================
//    // CORRECTION: Interpret as deltas -> convert to absolute stretch ratios
//    // ========================================================================
//    std::cout << "\n  *** LAMBDA SIGN CONVENTION CORRECTION ***" << std::endl;
//    std::cout << "  Detected negative center lambda(s) -- interpreting as strain deltas." << std::endl;
//    std::cout << "  Converting: lambda_actual = 1.0 + lambda_xml" << std::endl;
//    std::cout << std::endl;
//
//    auto fix = [](double& val, const char* name) {
//        double old = val;
//        val = 1.0 + val;
//        printf("    %s: %+.6f -> %.6f\n", name, old, val);
//    };
//
//    fix(gp.lambda_iso_center_outDV,   "lambda_iso_center_outDV  ");
//    fix(gp.lambda_iso_edge_outDV,     "lambda_iso_edge_outDV    ");
//    fix(gp.lambda_aniso_center_outDV, "lambda_aniso_center_outDV");
//    fix(gp.lambda_aniso_edge_outDV,   "lambda_aniso_edge_outDV  ");
//    fix(gp.lambda_iso_center_inDV,    "lambda_iso_center_inDV   ");
//    fix(gp.lambda_iso_edge_inDV,      "lambda_iso_edge_inDV     ");
//    fix(gp.lambda_aniso_center_inDV,  "lambda_aniso_center_inDV ");
//    fix(gp.lambda_aniso_edge_inDV,    "lambda_aniso_edge_inDV   ");
//
//    std::cout << std::endl;
//
//    // Validate: all corrected values should now be positive
//    bool valid = true;
//    if (gp.lambda_iso_center_outDV <= 0.0 || gp.lambda_aniso_center_outDV <= 0.0) {
//        std::cout << "  *** WARNING: Corrected outDV center lambdas still non-positive! ***" << std::endl;
//        std::cout << "  *** The XML values may not be strain deltas. Check your data source. ***" << std::endl;
//        valid = false;
//    }
//    if (gp.lambda_iso_edge_outDV <= 0.0 || gp.lambda_aniso_edge_outDV <= 0.0) {
//        std::cout << "  *** WARNING: Corrected outDV edge lambdas are non-positive! ***" << std::endl;
//        valid = false;
//    }
//
//    // Print the resulting strain profile at center and edge
//    double rr_center = gp.lambda_iso_center_outDV * gp.lambda_aniso_center_outDV;
//    double pp_center = gp.lambda_iso_center_outDV / gp.lambda_aniso_center_outDV;
//    double ss_center = 1.0 / (rr_center * pp_center);
//    double rr_edge = gp.lambda_iso_edge_outDV * gp.lambda_aniso_edge_outDV;
//    double pp_edge = gp.lambda_iso_edge_outDV / gp.lambda_aniso_edge_outDV;
//    double ss_edge = 1.0 / (rr_edge * pp_edge);
//
//    std::cout << "  Resulting outDV strain field:" << std::endl;
//    std::cout << "    Center: lambda_rr=" << rr_center << ", lambda_pp=" << pp_center 
//              << ", lambda_ss=" << ss_center << std::endl;
//    std::cout << "    Edge:   lambda_rr=" << rr_edge << ", lambda_pp=" << pp_edge 
//              << ", lambda_ss=" << ss_edge << std::endl;
//    std::cout << "    (lambda_ss = 1/(lambda_rr*lambda_pp) for volume conservation)" << std::endl;
//
//    return valid;
//}
//
//
//void buildVertexLambda(GeneralParams& gp,
//                       CoordInfoVecs& coord,
//                       LambdaField&   field,
//                       double         tFrac /*unused but kept for API parity*/)
//{
//    // size guards
//    if (field.lam_rr.size() != coord.nodeLocX.size())
//        field.resize(coord.nodeLocX.size());
//    if (gp.rho.size() != coord.nodeLocX.size())
//        gp.rho.resize(coord.nodeLocX.size());
//
//    const int N = static_cast<int>(coord.nodeLocX.size());
//    dim3 grid((N + BLOCK_SZ - 1) / BLOCK_SZ);
//
//    // ========================================================================
//    // AUTO-CORRECT lambda sign convention BEFORE building the field.
//    // This only modifies GeneralParams if negative values are detected.
//    // It is safe to call multiple times -- once corrected, values stay positive
//    // and the function becomes a no-op.
//    // ========================================================================
//    static bool lambda_corrected = false;
//    if (!lambda_corrected) {
//        lambda_corrected = correctLambdaSignConvention(gp);
//        // Even if correction returns false (already valid), mark as done
//        lambda_corrected = true;
//    }
//
//    // --- build lambda field (uses DV mask & the basis we just built) ---
//    k_buildLambda<<<grid,BLOCK_SZ>>>(
//        N,
//        thrust::raw_pointer_cast(coord.nodeLocX.data()),
//        thrust::raw_pointer_cast(coord.nodeLocY.data()),
//        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
//        gp.centerX, gp.centerY, gp.centerZ,
//        // outDV
//        gp.lambda_iso_center_outDV,   gp.lambda_iso_edge_outDV,
//        gp.lambda_aniso_center_outDV, gp.lambda_aniso_edge_outDV,
//        // inDV
//        gp.lambda_iso_center_inDV,    gp.lambda_iso_edge_inDV,
//        gp.lambda_aniso_center_inDV,  gp.lambda_aniso_edge_inDV,
//        // geometry
//        gp.disc_radius,
//        // outputs
//        thrust::raw_pointer_cast(gp.rho.data()),
//        thrust::raw_pointer_cast(field.lam_rr.data()),
//        thrust::raw_pointer_cast(field.lam_pp.data()),
//        thrust::raw_pointer_cast(field.lam_ss.data()),
//        thrust::raw_pointer_cast(field.e_R.data()),
//        thrust::raw_pointer_cast(field.e_phi.data()),
//        thrust::raw_pointer_cast(field.e_h.data()),
//        thrust::raw_pointer_cast(field.lam_alpha.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_upperhem.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_DV.data()) );
//    cudaDeviceSynchronize();
//}
//
//// ============================================================================
//// Helper: Determine strain mode from lambda values
//// ============================================================================
//static int determineStrainMode(const GeneralParams& gp)
//{
//    const double tol = 1e-12;
//
//    bool outDV_is_identity =
//        (fabs(gp.lambda_iso_center_outDV   - 1.0) < tol) &&
//        (fabs(gp.lambda_iso_edge_outDV     - 1.0) < tol) &&
//        (fabs(gp.lambda_aniso_center_outDV - 1.0) < tol) &&
//        (fabs(gp.lambda_aniso_edge_outDV   - 1.0) < tol);
//
//    bool inDV_is_identity =
//        (fabs(gp.lambda_iso_center_inDV    - 1.0) < tol) &&
//        (fabs(gp.lambda_iso_edge_inDV      - 1.0) < tol) &&
//        (fabs(gp.lambda_aniso_center_inDV  - 1.0) < tol) &&
//        (fabs(gp.lambda_aniso_edge_inDV    - 1.0) < tol);
//
//    int mode;
//    if (!inDV_is_identity && !outDV_is_identity) {
//        mode = 0;
//    } else if (!inDV_is_identity && outDV_is_identity) {
//        mode = 1;
//    } else if (inDV_is_identity && !outDV_is_identity) {
//        mode = 2;
//    } else {
//        mode = -1; 
//    }
//
//    const char* modeNames[] = {"BOTH", "inDV ONLY", "outDV ONLY", "NONE (identity)"};
//    int modeIdx = (mode == -1) ? 3 : mode;
//    std::cout << "  Strain mode: " << modeNames[modeIdx] 
//              << " (outDV_identity=" << outDV_is_identity 
//              << ", inDV_identity=" << inDV_is_identity << ")" << std::endl;
//
//    return mode;
//}
//
//void updateEdgeRestLengths(CoordInfoVecs&  coord,
//                           GeneralParams&  gp,
//                           const LambdaField& field,
//                           LinearSpringInfoVecs& lsInfo,
//                           int /*targetLayer, unused now*/)
//{
//    const int E = static_cast<int>(coord.num_edges);
//    dim3 grid((E + BLOCK_SZ - 1) / BLOCK_SZ);
//
//    int strainMode = determineStrainMode(gp);
//
//    if (strainMode == -1) {
//        std::cout << "  Skipping edge rest-length update (all lambdas are identity)." << std::endl;
//        thrust::copy(lsInfo.edge_initial_length.begin(),
//                     lsInfo.edge_initial_length.begin() + E,
//                     lsInfo.edge_final_length.begin());
//        return;
//    }
//
//    k_edgeRestProj<<<grid,BLOCK_SZ>>>(
//        E,
//        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()),
//        thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
//        thrust::raw_pointer_cast(coord.nodeLocX.data()),
//        thrust::raw_pointer_cast(coord.nodeLocY.data()),
//        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
//        thrust::raw_pointer_cast(field.lam_alpha.data()),
//        thrust::raw_pointer_cast(lsInfo.edge_initial_length.data()),
//        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
//        thrust::raw_pointer_cast(gp.edges_in_upperhem.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_DV.data()),
//        strainMode
//    );
//    cudaDeviceSynchronize();
//}
//
//} // namespace StrainTensorGPU