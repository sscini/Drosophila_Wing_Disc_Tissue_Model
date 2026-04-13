#pragma once
/******************************************************************************************************
 *  Strain tensor engine — rewritten to match Fuhrmann et al. (Science Advances, 2024) exactly.
 *
 *  What Fuhrmann does (Python, methods.py + array_wd.py):
 *
 *    1. Basis vectors (e_h, e_R, e_phi) are computed on the INITIAL spherical geometry.
 *       - e_h  = outward surface normal = position / |position|
 *       - e_R  = Gram-Schmidt projection of (vertex - origin) onto the tangent plane
 *       - e_phi = e_h × e_R
 *       Origins differ for inDV (DV midline) vs outDV (OD or OV centroid).
 *
 *    2. Normalized path length ? ? [0,1] is the geodesic distance on the sphere from
 *       the region-specific origin, divided by the max path length in that region.
 *
 *    3. ?_iso(?) and ?_aniso(?) are polynomial functions of ? (np.poly1d).
 *       For now we keep linear interpolation (2 coefficients = slope + intercept)
 *       to match your XML input format, but the infrastructure supports polynomials.
 *
 *    4. The full 3×3 tensor is:
 *         ? = ?_rr (e_R?e_R) + ?_ff (e_f?e_f) + ?_hh (e_h?e_h)
 *       with ?_rr = ?_iso×?_aniso,  ?_ff = ?_iso/?_aniso,  ?_hh = 1.0
 *
 *    5. For each edge (a,ß):
 *         ?_avg = 0.5 × (?_a + ?_ß)
 *         l*    = |?_avg · d|      where d = x_a - x_ß (initial edge vector)
 *       ALL edges go through this — vertical edges get ?_hh=1 so they don't change.
 *
 *  Author: Navaira Sherwani, 2025 (rewrite by Claude, March 2026)
 ******************************************************************************************************/

#include <thrust/device_vector.h>
#include "SystemStructures.h"
#include "System.h"

/* ================================================================================================ */
/*  Per-vertex growth field container                                                                */
/* ================================================================================================ */
struct LambdaField {

    // Orthonormal basis at each vertex
    thrust::device_vector<CVec3> e_h;     // surface normal
    thrust::device_vector<CVec3> e_R;     // in-surface radial direction
    thrust::device_vector<CVec3> e_phi;   // in-surface azimuthal direction

    // Diagonal components of ? in that basis
    thrust::device_vector<double> lam_rr; // radial      (?_iso × ?_aniso)
    thrust::device_vector<double> lam_pp; // circumferen. (?_iso / ?_aniso)
    thrust::device_vector<double> lam_ss; // thickness   (= 1.0, matching Fuhrmann)

    // Full 3×3 tensor in world coordinates
    thrust::device_vector<Mat_3x3> lam_alpha;

    // Normalised path length ? ? [0,1]
    thrust::device_vector<double> rho;

    inline void resize(std::size_t N) {
        rho.resize(N);
        lam_rr.resize(N);
        lam_pp.resize(N);
        lam_ss.resize(N);
        e_h.resize(N);
        e_R.resize(N);
        e_phi.resize(N);
        lam_alpha.resize(N);
    }
};

/* ================================================================================================ */
/*  Public GPU interface                                                                            */
/* ================================================================================================ */
namespace StrainTensorGPU {

    /**
     * Compute basis vectors with DV separation — called ONCE on initial geometry.
     * Stores results in coordInfoVecs (e_R_x/y/z, e_phi_x/y/z, e_h_x/y/z, pathlength_scaled)
     * and also sets params.nodes_in_DV.
     */
    void computeBasisVectorsWithDVSeparation(
        GeneralParams& params,
        CoordInfoVecs& coords,
        double theta_DV = 0.1931,
        double R = 1.0);

    /**
     * Build per-vertex ? field from current lambda parameters.
     * Uses the basis vectors already stored in the LambdaField (copied from coordInfoVecs).
     * ?_hh = 1.0 (no volume conservation in the growth tensor, matching Fuhrmann).
     */
    void buildVertexLambda(GeneralParams& gp,
                           CoordInfoVecs& coord,
                           LambdaField&   field,
                           double         tFrac);

    /**
     * Update ALL edge rest lengths by projecting the full ? tensor on each edge.
     *
     * Matching Fuhrmann exactly:
     *   ?_avg = 0.5 * (?_a + ?_ß)
     *   l*    = |?_avg · d_initial|
     *
     * ALL edges are processed (no vertical skip). Since ?_hh = 1, vertical
     * edges naturally get l* ˜ l_initial.
     *
     * Results go into lsInfo.edge_final_length.
     */
    void updateEdgeRestLengths(CoordInfoVecs&        coord,
                               GeneralParams&        gp,
                               const LambdaField&    field,
                               LinearSpringInfoVecs& lsInfo,
                               int                   targetLayer);

    // Unused but kept for ABI compatibility
    void updatePreferredAngles(BendingTriangleInfoVecs&  btiv,
                               const CoordInfoVecs&      coord,
                               const LambdaField&        field,
                               const GeneralParams&      gp,
                               const LinearSpringInfoVecs& lsInfo);

} // namespace StrainTensorGPU


// Previously working strain field commented out. 03/12/26






//#pragma once
///******************************************************************************************************
// *  Spontaneous-strain engine for planar / weak-curvature spring meshes
// *
// *  Implements the axi-symmetric growth tensor ?(x) used in the wing-eversion model.
// *
// *      1) StrainTensorGPU::buildVertexLambda   ? per-vertex bases (e_h,e_R,e_f) and ? tensor
// *      2) StrainTensorGPU::updateEdgeRestLengths ? edge rest lengths from full ? : (e ? e)
// *         (vertical/pillar edges are skipped when edgeLayerFlags == -1)
// *
// *  DV stripe membership is computed independently of layer: all stacks participate in DV.
// *
// *  Author: Navaira Sherwani, 2025
// ******************************************************************************************************/
//
//#include <thrust/device_vector.h>
//#include "SystemStructures.h"   // defines CVec3, Mat_3x3, etc.
//#include "System.h"             // GeneralParams, CoordInfoVecs, LinearSpringInfoVecs, ...
//
///* ================================================================================================ */
///*  Per-vertex growth field container                                                                */
///* ================================================================================================ */
//struct LambdaField {
//
//    // Orthonormal basis at each vertex
//    thrust::device_vector<CVec3> e_h;     // surface normal
//    thrust::device_vector<CVec3> e_R;     // in-surface radial direction
//    thrust::device_vector<CVec3> e_phi;   // in-surface azimuthal direction
//
//    // Diagonal components of ? in that basis
//    thrust::device_vector<double> lam_rr; // radial      (?_rr = ?_iso * ?_aniso)
//    thrust::device_vector<double> lam_pp; // circumferen. (?_ff = ?_iso / ?_aniso)
//    thrust::device_vector<double> lam_ss; // thickness   (?_hh = 1)
//
//    // Full 3×3 tensor in world coordinates (assembled from basis above)
//    thrust::device_vector<Mat_3x3> lam_alpha;
//
//    // Normalised planar radius ? = r / disc_radius, stored for diagnostics
//    thrust::device_vector<double> rho;
//
//    inline void resize(std::size_t N) {
//        rho.resize(N);
//        lam_rr.resize(N);
//        lam_pp.resize(N);
//        lam_ss.resize(N);
//        e_h.resize(N);
//        e_R.resize(N);
//        e_phi.resize(N);
//        lam_alpha.resize(N);
//    }
//};
//
///* ================================================================================================ */
///*  Public GPU interface (wrappers are defined in StrainTensor.cu)                                  */
///* ================================================================================================ */
//namespace StrainTensorGPU {
//    
//    
//    void computeBasisVectorsWithDVSeparation(
//        GeneralParams& params,
//        CoordInfoVecs& coords,
//        double theta_DV = 0.1931,  // ~11 degrees
//        double R = 1.0);
//    /**
//     * Build local bases and ? field at all vertices.
//     * Internally also marks the DV stripe (all layers participate in DV).
//     *
//     * Inputs:
//     *   gp     : General parameters (centers, disc_radius, theta_DV,
//     *            lambda_*_{inDV,outDV}, nodes_in_upperhem (used as edge flags), etc.)
//     *   coord  : Node/edge connectivity and positions
//     *   tFrac  : Pseudo-time fraction (kept for API parity; not used in current kernels)
//     *
//     * Outputs (written on device):
//     *   field.e_* , field.lam_* , field.lam_alpha , gp.rho , gp.nodes_in_DV
//     */
//    void buildVertexLambda(GeneralParams& gp,
//                           CoordInfoVecs& coord,
//                           LambdaField&   field,
//                           double         tFrac);
//
//    /**
//     * Update linear-spring rest lengths by projecting the full ? tensor on each edge.
//     * Vertical/pillar edges are not modified (they are identified by edgeLayerFlags == -1
//     * which you provide in gp.edges_in_upperhem for this call).
//     */
//    void updateEdgeRestLengths(CoordInfoVecs&        coord,
//                               GeneralParams&        gp,
//                               const LambdaField&    field,
//                               LinearSpringInfoVecs& lsInfo,
//                               int                   targetLayer /*kept for compatibility*/);
//
//    // Optional bending term (not used in your current runs). Declaration left for ABI compatibility.
//    void updatePreferredAngles(BendingTriangleInfoVecs&  btiv,
//                               const CoordInfoVecs&      coord,
//                               const LambdaField&        field,
//                               const GeneralParams&      gp,
//                               const LinearSpringInfoVecs& lsInfo);
//
//} // namespace StrainTensorGPU
