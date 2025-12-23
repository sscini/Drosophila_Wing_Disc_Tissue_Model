// StrainTensor.cu  — DV-aware basis construction and rest-length projection
//
// Spontaneous strain tensor ? = ?_rr (e_R ? e_R) + ?_ff (e_f ? e_f) + ?_hh (e_h ? e_h)
// with ?_rr = ?_iso * ?_aniso,  ?_ff = ?_iso / ?_aniso,  ?_hh = 1.
//
// Inside the DV stripe:
//   - e_R is locked to the DV axis (the blue line from one DV edge to the other).
// Outside the DV stripe:
//   - e_R is constructed exactly as before BUT using OA from OD (dorsal) or OV (ventral).
//   - The side (D vs V) is decided by the sign of projection onto a fixed in-plane vector
//     perpendicular to the DV axis.
//
// Vertical/pillar edges are flagged as -1 and are not altered here.

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
// Mark DV stripe (independent of layer).
// Stripe is |x - centerX| = R * sin(theta_DV/2).
// ============================================================================

__global__
void k_markDVstripe(int N,
                    const double* x,
                    double centerX, double R, double thetaDV,
                    const int* /*isUpper, unused for gating*/,
                    int* DVflag)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i>=N) return;

    double halfw = R * sin(0.5*thetaDV);
    DVflag[i] = (fabs(x[i] - centerX) <= halfw) ? 1 : 0;
}

// ============================================================================
// Build local basis vectors (e_h, e_R, e_phi) at every vertex
// DV-aware variant:
//   - inside DV stripe: e_R is the (global) unit DV axis
//   - outside stripe:   e_R from OD/OV edge-centers on the appropriate side
// ============================================================================

__global__
void k_buildBasis_DVaware(
                  int N,
                  const double* x, const double *y, const double *z,
                  // sphere center for e_h (normal)
                  double s_cx, double s_cy, double s_cz,
                  // legacy disc center (kept for ABI; not used here)
                  double /*cx*/, double /*cy*/, double /*cz*/,
                  // DV axis endpoints (edge-to-edge through the stripe)
                  CVec3 DV_A, CVec3 DV_B,
                  // Dorsal/Ventral edge-center origins
                  CVec3 O_D, CVec3 O_V,
                  // outputs
                  const int* DVflag,
                  CVec3* e_h,
                  CVec3* e_R,
                  CVec3* e_phi)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i>=N) return;

    // -------- per-vertex position
    CVec3 P(x[i], y[i], z[i]);

    // -------- surface normal e_h (use true sphere center)
    CVec3 PS = CVec3(x[i]-s_cx, y[i]-s_cy, z[i]-s_cz);
    CVec3 eh = normalize(PS);
    e_h[i] = eh;

    // -------- constant (global) DV axis and its in-plane perpendicular
    CVec3 u = normalize( CVec3( c<0>(DV_B)-c<0>(DV_A),
                                 c<1>(DV_B)-c<1>(DV_A),
                                 c<2>(DV_B)-c<2>(DV_A) ) );   // axis along stripe
    CVec3 Omid = CVec3( 0.5*(c<0>(DV_A)+c<0>(DV_B)),
                        0.5*(c<1>(DV_A)+c<1>(DV_B)),
                        0.5*(c<2>(DV_A)+c<2>(DV_B)) );

    // disc normal (approx) from sphere center to DV mid
    CVec3 n_disc = normalize( CVec3( c<0>(Omid)-s_cx, c<1>(Omid)-s_cy, c<2>(Omid)-s_cz ) );

    // in-plane perpendicular to u, pointing from DV axis toward (say) Ventral
    CVec3 v_perp = normalize( cross(n_disc, u) );  // n×u -> in-plane, ? to u

    // -------- choose construction by region
    CVec3 eR;
    if (DVflag[i]) {
        // Inside the stripe: e_R is locked to the DV axis
        eR = u;
    } else {
        // Outside stripe: decide side by signed distance along v_perp
        double sgn = dot3( CVec3( c<0>(P)-c<0>(Omid),
                                  c<1>(P)-c<1>(Omid),
                                  c<2>(P)-c<2>(Omid) ), v_perp );

        CVec3 O = (sgn >= 0.0) ? O_V : O_D;  // +side -> Ventral, -side -> Dorsal

        // OA from region origin, projected to tangent plane, normalized
        CVec3 OA = CVec3( c<0>(P)-c<0>(O), c<1>(P)-c<1>(O), c<2>(P)-c<2>(O) );
        double hdot = dot3(eh, OA);
        CVec3 OA_tan = CVec3( c<0>(OA) - hdot*c<0>(eh),
                              c<1>(OA) - hdot*c<1>(eh),
                              c<2>(OA) - hdot*c<2>(eh) );
        eR = normalize(OA_tan);
    }

    e_R[i] = eR;
    // right-handed complement
    CVec3 ephi = normalize( cross(eh, eR) );
    e_phi[i] = ephi;
}

// ============================================================================
// Build ? at vertices (unchanged except it consumes the e_R/e_phi/e_h from above)
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

    // Do not alter vertical/pillar edges
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

    // allocate (once) the DV mask
    if (gp.nodes_in_DV.size() != N) gp.nodes_in_DV.resize(N);

    // --- mark DV stripe: all layers participate in DV ---
    k_markDVstripe<<<grid,BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.nodeLocX.data()),
        gp.centerX,
        gp.disc_radius,      // radius used for stripe half-width
        gp.theta_DV,         // DV angular width (radians)
        thrust::raw_pointer_cast(gp.nodes_in_upperhem.data()), // kept for ABI
        thrust::raw_pointer_cast(gp.nodes_in_DV.data()) );
    cudaDeviceSynchronize();

    // Precompute constant tuples for the DV-aware basis kernel
    CVec3 DV_A(gp.DV_ax, gp.DV_ay, gp.DV_az);
    CVec3 DV_B(gp.DV_bx, gp.DV_by, gp.DV_bz);
    CVec3 O_D(gp.ODx,    gp.ODy,    gp.ODz);
    CVec3 O_V(gp.OVx,    gp.OVy,    gp.OVz);

    // --- build local bases on every vertex (DV-aware) ---
    k_buildBasis_DVaware<<<grid,BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.nodeLocX.data()),
        thrust::raw_pointer_cast(coord.nodeLocY.data()),
        thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        gp.c_dx, gp.c_dy, gp.c_dz,       // sphere center FOR e_h
        gp.centerX, gp.centerY, gp.centerZ, // (kept for ABI, not used here)
        DV_A, DV_B, O_D, O_V,
        thrust::raw_pointer_cast(gp.nodes_in_DV.data()),
        thrust::raw_pointer_cast(field.e_h.data()),
        thrust::raw_pointer_cast(field.e_R.data()),
        thrust::raw_pointer_cast(field.e_phi.data()) );
    cudaDeviceSynchronize();

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
