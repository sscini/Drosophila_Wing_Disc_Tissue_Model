#pragma once
/******************************************************************************************************
 *  Spontaneous-strain engine for planar / weak-curvature spring meshes                               *
 *                                                                                                    *
 *  Implements the axi-symmetric growth tensor lambda(x,t) used in the Science Advances               *
 *  2024 wing-eversion paper (Eqs 22-29).                                                             *
 *                                                                                                    *
 *      1. buildVertexLambda()     lambda  at every vertex (radial, circumf., through-thick)          *
 *      2. updateEdgeRestLengths()     lambda-projected edge rest lengths (Eq 26)                     *
 *      3. updatePreferredAngles()     optional anisotropic spontaneous curvature  (not used here)    *
 *                                                                                                    *
 *  Author: Navaira Sherwani, 2025                                                                    *
 ******************************************************************************************************/

#include <thrust/device_vector.h>
#include "SystemStructures.h"
#include "System.h"


//__host__ __device__
//inline Mat_3x3 makeIdentity3x3()
//{
//    return thrust::make_tuple<CVec3>(
//    /* row 1 */ thrust::make_tuple<double>(1.0, 0.0, 0.0),
//    /* row 2 */ thrust::make_tuple<double>(0.0, 1.0, 0.0),
//    /* row 3 */ thrust::make_tuple<double>(0.0, 0.0, 1.0)
//    );
//}
//
//__host__ __device__
//inline CVec3 makeOnes3()
//{
//    return thrust::make_tuple<double>(1.0, 1.0, 1.0);
//}


/* ============================================================================= */
/*  GPU container that stores three diagonal entries of ? at every vertex        */
/* ============================================================================= */
struct LambdaField {
    
    // make lambda field 
    
    // thrust::device_vector<double> lam_rr, lam_pp, lam_ss; // <- should be a matrix with diagonal entries
    
    // make basis vectors. per-vertex orthonormal basis, now as a 3x1 column vector. 
    
    thrust::device_vector<CVec3> e_h, e_R, e_phi; // 3x1 basis at v_i     

    thrust::device_vector<double> lam_rr;   ///< radial component   lambda_rr
    thrust::device_vector<double> lam_pp;   ///< circumf. component lambda_ff
    thrust::device_vector<double> lam_ss;   ///< thickness          lambda_ss
    
    thrust::device_vector<Mat_3x3> lam_alpha; // full lambda tensor 3x3
    
    thrust::device_vector<double> rho;      ///< normalized radius at vertex i 
    // basis vectors:
   // thrust::device_vector<double> e_h_x, e_h_y, e_h_z;
   // thrust::device_vector<double> e_R_x, e_R_y, e_R_z;
   // thrust::device_vector<double> e_phi_x, e_phi_y, e_phi_z;
    void resize(std::size_t N)
    {    
        rho.resize(N);
        lam_rr.resize(N);
        lam_pp.resize(N);
        lam_ss.resize(N);
        e_h.resize(N);// e_h_y.resize(N); e_h_z.resize(N);
        e_R.resize(N);// e_R_y.resize(N); e_R_z.resize(N);
        e_phi.resize(N);// e_phi_y.resize(N); e_phi_z.resize(N);
        lam_alpha.resize(N);
    }
};

/* ============================================================================= */
/*  Public interface (all runs on GPU; wrapper functions are in this namespace)  */
/* ============================================================================= */
namespace StrainTensorGPU { 

    /** Build the basis vectors and lambda field at the current pseudo-time fraction tFrac (T/Tf) at each vertex. */
    void buildVertexLambda(GeneralParams& gp,
                           CoordInfoVecs& coord,
                           LambdaField&         field,
                           double               tFrac);

    /** Update linear-spring rest lengths with edge-wise ?-projection
        (plain diagonal average). */
    void updateEdgeRestLengths(CoordInfoVecs&      coord,
                               GeneralParams& gp,
                               const LambdaField&        field,
                               LinearSpringInfoVecs&     lsInfo, 
                               int targetLayer);


    //--------------------------------------------------------------
    //  Preferred dihedral angles  ?0  for every bending triangle
    //--------------------------------------------------------------
    void updatePreferredAngles(                     //  <<<  HOST  >>>
            BendingTriangleInfoVecs&  btiv,   // t2e*, theta0
            const CoordInfoVecs& coord,     // e2n*, edgeLen
            const LambdaField& field,
            const GeneralParams& gp,
            const LinearSpringInfoVecs& lsInfo);     // gp.thickness

} // namespace StrainTensorGPU



//#pragma once
///******************************************************************************************************
// *  Spontaneous-strain engine for planar / weak-curvature spring meshes                                *
// *                                                                                                     *
// *  Implements the axi-symmetric growth tensor ?(x,t) used in the Science Advances                     *
// *  2024 wing-eversion paper (Eqs 22-29).                                                              *
// *                                                                                                     *
// *      1. buildVertexLambda()      ?  at every vertex (radial, circumf., through-thick)               *
// *      2. updateEdgeRestLengths()  ?-projected edge rest lengths (Eq 26)                              *
// *      3. updatePreferredAngles()  optional anisotropic spontaneous curvature (not used here)         *
// *                                                                                                     *
// *  Author: Navaira Sherwani, 2025                                                                     *
// ******************************************************************************************************/
//
//#include <thrust/device_vector.h>
//#include "SystemStructures.h"
//#include "System.h"
//
///* ------------------------------------------------------------------------- */
///*  CUDA kernel that tags vertices inside the DV boundary stripe              */
///* ------------------------------------------------------------------------- */
//namespace StrainTensorGPU {
//
///**
// * Mark vertices that lie within ±(frac·R) of the DV axis (y-direction) on the
// * apical layer.  Results are written to gp.nodes_in_DV (0 = outside, 1 = inside).
// *
// *  N             : total number of vertices
// *  x, y          : raw pointers to node coordinates
// *  cx, cy        : centre of apical layer
// *  R             : apical radius measured along +x axis
// *  frac          : half-width of stripe expressed as fraction of R (e.g. 0.05)
// *  upperFlag[i]  : 1 if vertex i belongs to the apical (upper) hemi-layer
// *  DVflag[i]     : output mask (0/1)
// */
////__global__
////void k_markDVstripe(int N,
////                    const double* __restrict__ x,
////                    const double* __restrict__ y,
////                    double cx, double cy,
////                    double R,   double frac,
////                    const int*  __restrict__ upperFlag,
////                    int*              DVflag);
////
////} // namespace StrainTensorGPU
//
//
///* ============================================================================= */
///*  GPU container that stores three diagonal entries of ? at every vertex        */
///* ============================================================================= */
//struct LambdaField {
//    /* basis vectors (per vertex) */
//    thrust::device_vector<CVec3> e_h, e_R, e_phi;        // 3-component each
//
//    /* diagonal components of lambda */
//    thrust::device_vector<double> lam_rr;   // radial
//    thrust::device_vector<double> lam_pp;   // circumferential
//    thrust::device_vector<double> lam_ss;   // thickness
//
//    /* full tensor (optional) and polar radius */
//    thrust::device_vector<Mat_3x3> lam_alpha;
//    thrust::device_vector<double>  rho;     // normalised radius
//
//    void resize(std::size_t N)
//    {
//        rho.resize(N);
//        lam_rr.resize(N);  lam_pp.resize(N);  lam_ss.resize(N);
//        e_h.resize(N);     e_R.resize(N);     e_phi.resize(N);
//        lam_alpha.resize(N);
//    }
//};
//
///* ============================================================================= */
///*  Public interface (all runs on GPU; wrapper functions are in this namespace)  */
///* ============================================================================= */
//namespace StrainTensorGPU {
//
///** Build the basis vectors and lambda field at the current pseudo-time
//    fraction tFrac (T/Tf) for every vertex.  Internally launches k_markDVstripe. */
//void buildVertexLambda(GeneralParams&     gp,
//                       CoordInfoVecs&     coord,
//                       LambdaField&       field,
//                       double             tFrac);
//
///** Update linear-spring rest lengths with edge-wise ? projection (Eq 26). */
//void updateEdgeRestLengths(CoordInfoVecs&         coord,
//                           GeneralParams&         gp,
//                           const LambdaField&     field,
//                           LinearSpringInfoVecs&  lsInfo,
//                           int                    targetLayer);
//
///** Preferred dihedral angles ?0 for every bending triangle. */
//void updatePreferredAngles( BendingTriangleInfoVecs&  btiv,
//                            const CoordInfoVecs&      coord,
//                            const LambdaField&        field,
//                            const GeneralParams&      gp,
//                            const LinearSpringInfoVecs& lsInfo );
//
//} // namespace StrainTensorGPU
