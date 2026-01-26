#pragma once
/******************************************************************************************************
 *  Spontaneous-strain engine for planar / weak-curvature spring meshes
 *
 *  Implements the axi-symmetric growth tensor ?(x) used in the wing-eversion model.
 *
 *      1) StrainTensorGPU::buildVertexLambda   ? per-vertex bases (e_h,e_R,e_f) and ? tensor
 *      2) StrainTensorGPU::updateEdgeRestLengths ? edge rest lengths from full ? : (e ? e)
 *         (vertical/pillar edges are skipped when edgeLayerFlags == -1)
 *
 *  DV stripe membership is computed independently of layer: all stacks participate in DV.
 *
 *  Author: Navaira Sherwani, 2025
 ******************************************************************************************************/

#include <thrust/device_vector.h>
#include "SystemStructures.h"   // defines CVec3, Mat_3x3, etc.
#include "System.h"             // GeneralParams, CoordInfoVecs, LinearSpringInfoVecs, ...

/* ================================================================================================ */
/*  Per-vertex growth field container                                                                */
/* ================================================================================================ */
struct LambdaField {

    // Orthonormal basis at each vertex
    thrust::device_vector<CVec3> e_h;     // surface normal
    thrust::device_vector<CVec3> e_R;     // in-surface radial direction
    thrust::device_vector<CVec3> e_phi;   // in-surface azimuthal direction

    // Diagonal components of ? in that basis
    thrust::device_vector<double> lam_rr; // radial      (?_rr = ?_iso * ?_aniso)
    thrust::device_vector<double> lam_pp; // circumferen. (?_ff = ?_iso / ?_aniso)
    thrust::device_vector<double> lam_ss; // thickness   (?_hh = 1)

    // Full 3×3 tensor in world coordinates (assembled from basis above)
    thrust::device_vector<Mat_3x3> lam_alpha;

    // Normalised planar radius ? = r / disc_radius, stored for diagnostics
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
/*  Public GPU interface (wrappers are defined in StrainTensor.cu)                                  */
/* ================================================================================================ */
namespace StrainTensorGPU {

    /**
     * Build local bases and ? field at all vertices.
     * Internally also marks the DV stripe (all layers participate in DV).
     *
     * Inputs:
     *   gp     : General parameters (centers, disc_radius, theta_DV,
     *            lambda_*_{inDV,outDV}, nodes_in_upperhem (used as edge flags), etc.)
     *   coord  : Node/edge connectivity and positions
     *   tFrac  : Pseudo-time fraction (kept for API parity; not used in current kernels)
     *
     * Outputs (written on device):
     *   field.e_* , field.lam_* , field.lam_alpha , gp.rho , gp.nodes_in_DV
     */
    void buildVertexLambda(GeneralParams& gp,
                           CoordInfoVecs& coord,
                           LambdaField&   field,
                           double         tFrac);

    /**
     * Update linear-spring rest lengths by projecting the full ? tensor on each edge.
     * Vertical/pillar edges are not modified (they are identified by edgeLayerFlags == -1
     * which you provide in gp.edges_in_upperhem for this call).
     */
    void updateEdgeRestLengths(CoordInfoVecs&        coord,
                               GeneralParams&        gp,
                               const LambdaField&    field,
                               LinearSpringInfoVecs& lsInfo,
                               int                   targetLayer /*kept for compatibility*/);

    // Optional bending term (not used in your current runs). Declaration left for ABI compatibility.
    void updatePreferredAngles(BendingTriangleInfoVecs&  btiv,
                               const CoordInfoVecs&      coord,
                               const LambdaField&        field,
                               const GeneralParams&      gp,
                               const LinearSpringInfoVecs& lsInfo);

} // namespace StrainTensorGPU
