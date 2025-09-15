/*

      SPONTANEOUS STRAIN TENSOR
      
      ? = | ?11 ?12 ?13 |
          | ?21 ?22 ?23 |
          | ?31 ?32 ?33 |
          
      ? = ?11 (e1 x e1) + ?22 (e2 x e2) + ?33 (e3 x e3)   -- e1 and e2 are surface tangents. 
      
      ?11 = ?_iso * ?_aniso
      ?22 = ?_iso/?_aniso
      
      ? = ?_iso*?_aniso (e1 x e1) + ?_iso/?_aniso (e2 x e2) + (e3 x e3) <- THIS IS OUR LAMBDA FIELD. 
      
      1. So start by taking user values for ?_iso and ?_aniso
      
      2. Create the strain field. 
      
      3. apply the strain field to the vertices of the edge in question so:
                
                
                ?_a = 0.5 * [ ?*(X_a) + ?*(X_b)] (a and b are vertices of the spring a)
                
      4. Calculate the length of the spring and store in L0. 
             
                da = || X_a - X_b || = del_Xa -> spring length
      
      5. Apply the decomposed ?_a to (not the scalar) spring length.  
                
                da_F = || ?_a . del_Xa ||
                
                This should be then stored in edge_current_length. 
                
                The reason this is not edge_final_length is because the strain tensor is divided
                into segments so that the model can be relaxed incrementally.                                                                                                 
                
                       

*/

/*

For Gradient Descent - currently the LinearSpring files compute the forces on all the nodes. To implement gradient descent with overdamped dynamics add those to the node positions with tau_over_gamma as a variable that helps with damping. 

*/





#include "StrainTensor.h"
#include "SystemStructures.h"
#include "System.h"
#include <cuda_runtime.h>
#include <cmath>

#define BLOCK_SZ 256


// tuple helpers 

template<int I>  __host__ __device__
inline double c(const CVec3& v) { return thrust::get<I>(v); }

__host__ __device__
inline Mat_3x3 outer(const CVec3& a, const CVec3& b)
{
    return Mat_3x3(
        CVec3( c<0>(a)*c<0>(b), c<0>(a)*c<1>(b), c<0>(a)*c<2>(b) ),
        CVec3( c<1>(a)*c<0>(b), c<1>(a)*c<1>(b), c<1>(a)*c<2>(b) ),
        CVec3( c<2>(a)*c<0>(b), c<2>(a)*c<1>(b), c<2>(a)*c<2>(b) )
    );
}


// scaled add: C+= s * A (element-wise)

__host__ __device__
inline void axpy(double s, const Mat_3x3& A, Mat_3x3& C)
{
    /* row 0 */
    CVec3 r0 = thrust::get<0>(C);
    CVec3 a0 = thrust::get<0>(A);
    thrust::get<0>(r0) += s * c<0>(a0);
    thrust::get<1>(r0) += s * c<1>(a0);
    thrust::get<2>(r0) += s * c<2>(a0);
    thrust::get<0>(C)   = r0;

    /* row 1 */
    CVec3 r1 = thrust::get<1>(C);
    CVec3 a1 = thrust::get<1>(A);
    thrust::get<0>(r1) += s * c<0>(a1);
    thrust::get<1>(r1) += s * c<1>(a1);
    thrust::get<2>(r1) += s * c<2>(a1);
    thrust::get<1>(C)   = r1;

    /* row 2 */
    CVec3 r2 = thrust::get<2>(C);
    CVec3 a2 = thrust::get<2>(A);
    thrust::get<0>(r2) += s * c<0>(a2);
    thrust::get<1>(r2) += s * c<1>(a2);
    thrust::get<2>(r2) += s * c<2>(a2);
    thrust::get<2>(C)   = r2;
}



/* 3×3 · 3×1  ------------------------------------------------- */
__host__ __device__
inline CVec3 matVec(const Mat_3x3& M, const CVec3& v)
{
    return CVec3(
        c<0>( thrust::get<0>(M) )*c<0>(v) + c<1>( thrust::get<0>(M) )*c<1>(v) + c<2>( thrust::get<0>(M) )*c<2>(v),
        c<0>( thrust::get<1>(M) )*c<0>(v) + c<1>( thrust::get<1>(M) )*c<1>(v) + c<2>( thrust::get<1>(M) )*c<2>(v),
        c<0>( thrust::get<2>(M) )*c<0>(v) + c<1>( thrust::get<2>(M) )*c<1>(v) + c<2>( thrust::get<2>(M) )*c<2>(v)
    );
}


/* ?v?  ------------------------------------------------------- */
__host__ __device__
double norm3(const CVec3& v)
{
    return sqrt(thrust::get<0>(v)*thrust::get<0>(v) +
                thrust::get<1>(v)*thrust::get<1>(v) +
                thrust::get<2>(v)*thrust::get<2>(v)) + 1e-14;
}



/*----------------------------------------------------------------------------------
1. Function to build basis vectors
-----------------------------------------------------------------------------------*/

__global__
void k_buildBasis(int N, const double* x, const double *y, const double *z,
                  double c_dx, double c_dy, double c_dz, // projected center of full sphere 
                  double cx, double cy, double cz,// center of layer being computed
                  CVec3* e_h,// double* e_h_y, double* e_h_z, 
                  CVec3* e_R,// double* e_R_y, double* e_R_z, 
                  CVec3* e_phi)//, double* e_phi_y, double* e_phi_z)
{
                  
    int i = blockIdx.x*blockDim.x + threadIdx.x;
                      
    if (i>=N) return;
    
    // position and normal - eh vector
    double px = x[i], py = y[i], pz = z[i]; // c_d is the center of the completed sphere. 
    double nrm = sqrt(px*px + py*py + pz*pz) + 1e-14;
    double hx = px/nrm, hy = py/nrm, hz = pz/nrm;
    e_h[i] = CVec3(hx, hy, hz);  //hx; e_h_y[i] = hy; e_h_z[i] = hz;  
    
    // vector OA - where O is center of disc (not center of sphere) 
    double ox = x[i] - cx, oy = y[i] - cy, oz = z[i] - cz;
    double on = sqrt(ox*ox + oy*oy + oz*oz) +1e-14;
    ox/=on; oy/=on; oz/=on;
    
    // R = OA - (h.OA) h
    double dot_h_OA = hx*ox + hy*oy + hz*oz;
    double rx = ox - dot_h_OA*hx;
    double ry = oy - dot_h_OA*hy;
    double rz = oz - dot_h_OA*hz;
    double rn = sqrt(rx*rx + ry*ry + rz*rz) +1e-14;
    rx/=rn; ry/=rn; rz/=rn;
    
    e_R[i] = CVec3(rx, ry, rz);// rx; e_R_y[i] = ry; e_R_z[i] = rz;
    
    // e_phi = h x R
    double phix = hy*rz - hz*ry; 
    double phiy = hz*rx - hx*rz;
    double phiz = hx*ry - hy*rx;
    double phi_n = sqrt(phix*phix + phiy*phiy + phiz*phiz) +1e-14;
    e_phi[i] = CVec3(phix/phi_n, phiy/phi_n, phiz/phi_n);
                  
                  
}
                  

/* ----------------------------------------------------------------------------
   1.  build lambda at vertices  (one thread per vertex)                           */
__global__
void k_buildLambda(int    N,
                   const double *x,
                   const double *y,
                   const double *z,
                   double cx, double cy, double cz,          // mesh centre (GP)
                   double lam_iso_outDV_center, double lam_iso_outDV_edge, double lam_aniso_outDV_center, double lam_aniso_outDV_edge, double disc_radius, double *rho,   // Strain tensor field parameters from System.h generalParams +disc Radius 
                   double tFrac,
                   double *lam_rr,
                   double *lam_pp,
                   double *lam_ss, 
                   const CVec3* e_R, const CVec3* e_phi, const CVec3* e_h, 
                   Mat_3x3* lam_alpha, const int *layerFlag)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;
    
    //if (layerFlag == 0 || layerflag == 1) return;
    
    /*----  axi-symmetric basis vectors (flat sheet approximation)  ----*/
    double dx = x[tid] - cx;
    double dy = y[tid] - cy;
    double dz = z[tid] - cz;
    
    double r = sqrt(dx*dx + dy*dy) + 1e-14;
     
    rho[tid] = r/disc_radius;
                
    (void)r; // not used here but left for consistency / extensions

    /* linear schedule  ?(t) = I + e t                                          */
    double lamIso = (lam_iso_outDV_center + (lam_iso_outDV_edge - lam_iso_outDV_center)*(rho[tid]*rho[tid]));// * tFrac;   // radial   (_rr)
    double lamAni = (lam_aniso_outDV_center + (lam_aniso_outDV_edge - lam_aniso_outDV_center)*(rho[tid]*rho[tid])); //* tFrac;   // circumf. (_ff)
    
    //bool isBasal = (layerFlag[tid] < 0);
    //if (isBasal) lamAni = 1.0/lamAni;
    
    lam_rr[tid] = (lamIso*lamAni);
    lam_pp[tid] = (lamIso/lamAni);
    lam_ss[tid] = 1.0;                     // no through-thickness growth
    
    
    // Tensor at v_i
    
    Mat_3x3 L = Mat_3x3{
        CVec3(0,0,0), CVec3(0,0,0), CVec3(0,0,0),
    };
    
    axpy(lam_rr[tid], outer(e_R [tid], e_R [tid]), L);
    axpy(lam_pp[tid], outer(e_phi[tid], e_phi[tid]), L);
    axpy(lam_ss[tid], outer(e_h  [tid], e_h  [tid]), L);

    lam_alpha[tid] = L;   // store (row-major) 3×3 tensor
    
}



/* ===== Full projection lambda:e_e ============================ */
__global__
void k_edgeRestProj(int    E,
                    const int    *e2n1,  const int *e2n2,
                    const double *x,     const double *y,   const double *z,
                    const Mat_3x3 *lam_alpha,          // NEW  ? one ? per vertex
                    double *L0,                  // ? original, still here
                    double *Lstar,
                    const int *edgeLayerFlags, int targetLayer)
{
    int eid = blockIdx.x * blockDim.x + threadIdx.x;
    if (eid >= E) return;
    
    // skip edges if not in the desired layer
    if(edgeLayerFlags[eid]==2 ) return;

    int a = e2n1[eid]; 
    int b = e2n2[eid];
    
    CVec3 dX = CVec3(x[a] - x[b], y[a] - y[b], z[a] - z[b]);
    
    // store initial stretch in edge_initial_length
    L0[eid] = norm3(dX);
    
    Mat_3x3 La = lam_alpha[a];
    Mat_3x3 Lb = lam_alpha[b];
   Mat_3x3 Lp;
thrust::get<0>(Lp) = CVec3(
    0.5*( c<0>( thrust::get<0>(La) ) + c<0>( thrust::get<0>(Lb) )),
    0.5*( c<1>( thrust::get<0>(La) ) + c<1>( thrust::get<0>(Lb) )),
    0.5*( c<2>( thrust::get<0>(La) ) + c<2>( thrust::get<0>(Lb) )) );

thrust::get<1>(Lp) = CVec3(
    0.5*( c<0>( thrust::get<1>(La) ) + c<0>( thrust::get<1>(Lb) )),
    0.5*( c<1>( thrust::get<1>(La) ) + c<1>( thrust::get<1>(Lb) )),
    0.5*( c<2>( thrust::get<1>(La) ) + c<2>( thrust::get<1>(Lb) )) );

thrust::get<2>(Lp) = CVec3(
    0.5*( c<0>( thrust::get<2>(La) ) + c<0>( thrust::get<2>(Lb) )),
    0.5*( c<1>( thrust::get<2>(La) ) + c<1>( thrust::get<2>(Lb) )),
    0.5*( c<2>( thrust::get<2>(La) ) + c<2>( thrust::get<2>(Lb) )) );

     /* ---- stretched edge vector ----------------------------- */
    CVec3 dX_stretch = matVec(Lp, dX);

    /* ---- new rest length ----------------------------------- */
    Lstar[eid] = norm3(dX_stretch);

//    /* edge direction */
//    double ex = x[b]-x[a], ey = y[b]-y[a], ez = z[b]-z[a];
//    double len = sqrt(ex*ex + ey*ey + ez*ez) + 1e-14;
//    ex/=len;  ey/=len;  ez/=len;
//
//    /* average ? diagonal */
//    double lamR = 0.5*(lam_rr[a]+lam_rr[b]);
//    double lamP = 0.5*(lam_pp[a]+lam_pp[b]);
//    double lamS = 0.5*(lam_ss[a]+lam_ss[b]);
//
//
//    double erx = e_R_x[a], ery = e_R_y[a], erz = e_R_z[a];
//    double epx = e_phi_x[a], epy = e_phi_y[a], epz = e_phi_z[a];
//    double ehx = e_h_x[a], ehy = e_h_y[a], ehz = e_h_z[a];
//    
//
////    /* local polar basis at vertex *a* (flat sheet) */
////    double rlen = sqrt(x[a]*x[a] + y[a]*y[a]) + 1e-14;
////    double erx = x[a]/rlen, ery = y[a]/rlen;          // e_r
////    double epx = -ery,    epy = erx;                  // e_f
//
//    double dR = ex*erx + ey*ery + ez*erz;   // e*e_r
//    double dP = ex*epx + ey*epy + ez*epz;   // e*e_f
//    double dH = ex*ehx + ey*ehy + ez*ehz;   // e*e_h
//
//    double lamProj = lamR*dR*dR + lamP*dP*dP + lamS*dH*dH;
//    Lstar[eid] = lamProj * L0[eid];
}
//
///* ----------------------------------------------------------------------------
//   3.  preferred bending angle (anisotropic growth, Eq?25)                    */
//__global__
//void k_prefAngles(int T,
//                  const int    *t2e1, const int *t2e2, const int *t2e3,
//                  const int    *e2n1, const int *e2n2,
//                  const double *lam_rr, const double *lam_pp,
//                  const double *edgeLen,
//                  double thickness,
//                  double* theta0)
//{
//    //--------------------
//    // 1. collect vertices
//    //--------------------
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (tid >= T) return;
//    
//    int e1 = t2e1[tid], e2 = t2e2[tid], e3 = t2e3[tid];
//
//    int v1 = e2n1[e1];
//    int v2 = e2n2[e1];
//    /* pick the third vertex: it must be the one that is NOT v1 or v2 */
//    int v3 = (e2n1[e2]!=v1 && e2n1[e2]!=v2) ? e2n1[e2] :
//             (e2n2[e2]!=v1 && e2n2[e2]!=v2) ? e2n2[e2] :
//             (e2n1[e3]!=v1 && e2n1[e3]!=v2) ? e2n1[e3] : e2n2[e3];
//
//    //-----------------------------
//    // 2. average ? over the triple
//    //-----------------------------
//    double lamR = (lam_rr[v1] + lam_rr[v2] + lam_rr[v3]) / 3.0;
//    double lamP = (lam_pp[v1] + lam_pp[v2] + lam_pp[v3]) / 3.0;
//
//    //--------------------------------------------
//    // 3. spontaneous curvature   C0  (half-plate)
//    //--------------------------------------------
//    double C0 = 0.5 * (lamP - lamR) / thickness;   // ? ½ factor
//
//    /* mean surrounding edge length (more robust than triangle edge-specific) */
//    double Lm = (edgeLen[e1] + edgeLen[e2] + edgeLen[e3]) / 3.0;
//
//    theta0[tid] = C0 * Lm;                     // dihedral offset
//}

/* ----------------------------------------------------------------------------
   ====   PUBLIC WRAPPERS   ================================================== */
namespace StrainTensorGPU {

void buildVertexLambda(GeneralParams& gp,
                       CoordInfoVecs& coord,
                       LambdaField&         field,
                       double               tFrac)
{
    if (field.lam_rr.size() != coord.nodeLocX.size())
        field.resize(coord.nodeLocX.size());
    if (gp.rho.size() != coord.nodeLocX.size())
        gp.rho.resize(coord.nodeLocX.size());

    int N = static_cast<int>(coord.nodeLocX.size());
    dim3 grid((N + BLOCK_SZ - 1) / BLOCK_SZ);

    
    // build basis vectors
    k_buildBasis<<<grid,BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.nodeLocX.data()), thrust::raw_pointer_cast(coord.nodeLocY.data()), thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        gp.c_dx, gp.c_dy, gp.c_dz,
        gp.centerX, gp.centerY, gp.centerZ,
        thrust::raw_pointer_cast(field.e_h.data()), thrust::raw_pointer_cast(field.e_R.data()), thrust::raw_pointer_cast(field.e_phi.data()));
        
        cudaDeviceSynchronize();
    
    // build lambda field 
    k_buildLambda<<<grid,BLOCK_SZ>>>(
        N,
        thrust::raw_pointer_cast(coord.nodeLocX.data()), thrust::raw_pointer_cast(coord.nodeLocY.data()), thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        gp.centerX, gp.centerY, gp.centerZ,
        gp.lambda_iso_center_outDV, gp.lambda_iso_edge_outDV,
        gp.lambda_aniso_center_outDV, gp.lambda_aniso_edge_outDV,
        gp.disc_radius,
        thrust::raw_pointer_cast(gp.rho.data()),
        tFrac,
        thrust::raw_pointer_cast(field.lam_rr.data()), thrust::raw_pointer_cast(field.lam_pp.data()), thrust::raw_pointer_cast(field.lam_ss.data()),
        thrust::raw_pointer_cast(field.e_R.data()), thrust::raw_pointer_cast(field.e_phi.data()), thrust::raw_pointer_cast(field.e_h.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()),
        thrust::raw_pointer_cast(gp.edges_in_upperhem.data()));
    cudaDeviceSynchronize();
}

/* ------------------------------------------------------------------------- */
void updateEdgeRestLengths(CoordInfoVecs&  coord,
                           GeneralParams& gp,
                           const LambdaField&    field,
                           LinearSpringInfoVecs& lsInfo, int targetLayer)
{
    int E = static_cast<int>(coord.num_edges);
    dim3 grid((E + BLOCK_SZ - 1) / BLOCK_SZ); 


  /* full lambda:e_e projection */
    k_edgeRestProj<<<grid,BLOCK_SZ>>>(
        E,
        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()), thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
        thrust::raw_pointer_cast(coord.nodeLocX.data()),      thrust::raw_pointer_cast(coord.nodeLocY.data()),      thrust::raw_pointer_cast(coord.nodeLocZ.data()),
        thrust::raw_pointer_cast(field.lam_alpha.data()),        // NEW — whole tensor array
        thrust::raw_pointer_cast(lsInfo.edge_initial_length.data()),
        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
        thrust::raw_pointer_cast(gp.edges_in_upperhem.data()),
        targetLayer );

    cudaDeviceSynchronize();
}
//
////==============================================================
////  StrainTensorGPU::updatePreferredAngles (in .cu file)
////==============================================================
//void updatePreferredAngles(BendingTriangleInfoVecs& btiv,
//                           const CoordInfoVecs& coord,
//                           const LambdaField& field,
//                           const GeneralParams& gp,
//                           const LinearSpringInfoVecs& lsInfo)
//{
//    int  T      = coord.triangles2Edges_1.size();
//    int  TPB    = 128;
//    int  nBlk   = (T + TPB - 1) / TPB;
//    btiv.initial_angle.resize(T);
//
//    k_prefAngles<<<nBlk, TPB>>>(
//        T,
//        thrust::raw_pointer_cast(coord.triangles2Edges_1.data()),
//        thrust::raw_pointer_cast(coord.triangles2Edges_2.data()),
//        thrust::raw_pointer_cast(coord.triangles2Edges_3.data()),
//        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()),
//        thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
//        thrust::raw_pointer_cast(field.lam_rr.data()),
//        thrust::raw_pointer_cast(field.lam_pp.data()),
//        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
//        gp.thickness,
//        thrust::raw_pointer_cast(btiv.initial_angle.data()) );
//
//   // CUDA_CHECK_ERROR();
//}
//

} // namespace StrainTensorGPU





///*
//
//
//      SPONTANEOUS STRAIN TENSOR
//      
//      ? = | ?11 ?12 ?13 |
//          | ?21 ?22 ?23 |
//          | ?31 ?32 ?33 |
//          
//      ? = ?11 (e1 x e1) + ?22 (e2 x e2) + ?33 (e3 x e3)   -- e1 and e2 are surface tangents. 
//      
//      ?11 = ?_iso * ?_aniso
//      ?22 = ?_iso/?_aniso
//      
//      ? = ?_iso*?_aniso (e1 x e1) + ?_iso/?_aniso (e2 x e2) + (e3 x e3) <- THIS IS OUR LAMBDA FIELD. 
//      
//      1. So start by taking user values for ?_iso and ?_aniso
//      
//      2. Create the strain field. 
//      
//      3. apply the strain field to the vertices of the edge in question so:
//                
//                
//                ?_a = 0.5 * [ ?*(X_a) + ?*(X_b)] (a and b are vertices of the spring a)
//                
//      4. Calculate the length of the spring and store in L0. 
//             
//                da = || X_a - X_b || = del_Xa -> spring length
//      
//      5. Apply the decomposed ?_a to (not the scalar) spring length.  
//                
//                da_F = || ?_a . del_Xa ||
//                
//                This should be then stored in edge_current_length. 
//                
//                The reason this is not edge_final_length is because the strain tensor is divided
//                into segments so that the model can be relaxed incrementally.                                                                                                 
//                
//                       
//
//*/
//
///*
//
//For Gradient Descent - currently the LinearSpring files compute the forces on all the nodes. To implement gradient descent with overdamped dynamics add those to the node positions with tau_over_gamma as a variable that helps with damping. 
//
//*/
//
//
//
//
//
//#include "StrainTensor.h"
//#include "SystemStructures.h"
//#include "System.h"
//#include <cuda_runtime.h>
//#include <cmath>
//
//#define BLOCK_SZ 256
//
//
//// tuple helpers 
//
//template<int I>  __host__ __device__
//inline double c(const CVec3& v) { return thrust::get<I>(v); }
//
//__host__ __device__
//inline Mat_3x3 outer(const CVec3& a, const CVec3& b)
//{
//    return Mat_3x3(
//        CVec3( c<0>(a)*c<0>(b), c<0>(a)*c<1>(b), c<0>(a)*c<2>(b) ),
//        CVec3( c<1>(a)*c<0>(b), c<1>(a)*c<1>(b), c<1>(a)*c<2>(b) ),
//        CVec3( c<2>(a)*c<0>(b), c<2>(a)*c<1>(b), c<2>(a)*c<2>(b) )
//    );
//}
//
//
//// scaled add: C+= s * A (element-wise)
//
//__host__ __device__
//inline void axpy(double s, const Mat_3x3& A, Mat_3x3& C)
//{
//    /* row 0 */
//    CVec3 r0 = thrust::get<0>(C);
//    CVec3 a0 = thrust::get<0>(A);
//    thrust::get<0>(r0) += s * c<0>(a0);
//    thrust::get<1>(r0) += s * c<1>(a0);
//    thrust::get<2>(r0) += s * c<2>(a0);
//    thrust::get<0>(C)   = r0;
//
//    /* row 1 */
//    CVec3 r1 = thrust::get<1>(C);
//    CVec3 a1 = thrust::get<1>(A);
//    thrust::get<0>(r1) += s * c<0>(a1);
//    thrust::get<1>(r1) += s * c<1>(a1);
//    thrust::get<2>(r1) += s * c<2>(a1);
//    thrust::get<1>(C)   = r1;
//
//    /* row 2 */
//    CVec3 r2 = thrust::get<2>(C);
//    CVec3 a2 = thrust::get<2>(A);
//    thrust::get<0>(r2) += s * c<0>(a2);
//    thrust::get<1>(r2) += s * c<1>(a2);
//    thrust::get<2>(r2) += s * c<2>(a2);
//    thrust::get<2>(C)   = r2;
//}
//
//
//
///* 3×3 · 3×1  ------------------------------------------------- */
//__host__ __device__
//inline CVec3 matVec(const Mat_3x3& M, const CVec3& v)
//{
//    return CVec3(
//        c<0>( thrust::get<0>(M) )*c<0>(v) + c<1>( thrust::get<0>(M) )*c<1>(v) + c<2>( thrust::get<0>(M) )*c<2>(v),
//        c<0>( thrust::get<1>(M) )*c<0>(v) + c<1>( thrust::get<1>(M) )*c<1>(v) + c<2>( thrust::get<1>(M) )*c<2>(v),
//        c<0>( thrust::get<2>(M) )*c<0>(v) + c<1>( thrust::get<2>(M) )*c<1>(v) + c<2>( thrust::get<2>(M) )*c<2>(v)
//    );
//}
//
//
///* ?v?  ------------------------------------------------------- */
//__host__ __device__
//double norm3(const CVec3& v)
//{
//    return sqrt(thrust::get<0>(v)*thrust::get<0>(v) +
//                thrust::get<1>(v)*thrust::get<1>(v) +
//                thrust::get<2>(v)*thrust::get<2>(v)) + 1e-14;
//}
//
//
//                  
//__global__
//void k_markDVstripe(int N,
//                    const double* x, const double* y,
//                    double cx, double cy,
//                    double R, double frac,
//                    const int* upperFlag,
//                    int* DVflag)
//{
//    int i = blockIdx.x*blockDim.x + threadIdx.x;
//    if(i>=N) return;
//
//    // only apical vertices participate
//    if(!upperFlag[i]){ DVflag[i]=0; return; }
//
//    // signed distance to x–axis passing through the centre
//    double distAxis = fabs(y[i]-cy);
//
//    DVflag[i] = (distAxis <= frac*R) ? 1 : 0;
//}
//
//
///*----------------------------------------------------------------------------------
//1. Function to build basis vectors
//-----------------------------------------------------------------------------------*/
//
//__global__
//void k_buildBasis(int N, const double* x, const double *y, const double *z,
//                  double c_dx, double c_dy, double c_dz, // projected center of full sphere 
//                  double cx, double cy, double cz,// center of layer being computed
//                  CVec3* e_h,// double* e_h_y, double* e_h_z, 
//                  CVec3* e_R,// double* e_R_y, double* e_R_z, 
//                  CVec3* e_phi)//, double* e_phi_y, double* e_phi_z)
//{
//                  
//    int i = blockIdx.x*blockDim.x + threadIdx.x;
//                      
//    if (i>=N) return;
//    
//    // position and normal - eh vector
//    double px = x[i], py = y[i], pz = z[i]; // c_d is the center of the completed sphere. 
//    double nrm = sqrt(px*px + py*py + pz*pz) + 1e-14;
//    double hx = px/nrm, hy = py/nrm, hz = pz/nrm;
//    e_h[i] = CVec3(hx, hy, hz);  //hx; e_h_y[i] = hy; e_h_z[i] = hz;  
//    
//    // vector OA - where O is center of disc (not center of sphere) 
//    double ox = x[i] - cx, oy = y[i] - cy, oz = z[i] - cz;
//    double on = sqrt(ox*ox + oy*oy + oz*oz) +1e-14;
//    ox/=on; oy/=on; oz/=on;
//    
//    // R = OA - (h.OA) h
//    double dot_h_OA = hx*ox + hy*oy + hz*oz;
//    double rx = ox - dot_h_OA*hx;
//    double ry = oy - dot_h_OA*hy;
//    double rz = oz - dot_h_OA*hz;
//    double rn = sqrt(rx*rx + ry*ry + rz*rz) +1e-14;
//    rx/=rn; ry/=rn; rz/=rn;
//    
//    e_R[i] = CVec3(rx, ry, rz);// rx; e_R_y[i] = ry; e_R_z[i] = rz;
//    
//    // e_phi = h x R
//    double phix = hy*rz - hz*ry; 
//    double phiy = hz*rx - hx*rz;
//    double phiz = hx*ry - hy*rx;
//    double phi_n = sqrt(phix*phix + phiy*phiy + phiz*phiz) +1e-14;
//    e_phi[i] = CVec3(phix/phi_n, phiy/phi_n, phiz/phi_n);
//                  
//                  
//}
//
//
///* ----------------------------------------------------------------------------
//   1.  build lambda at vertices  (one thread per vertex)                           */
//__global__
//void k_buildLambda(int    N,
//                   const double *x,
//                   const double *y,
//                   const double *z,
//                   double cx, double cy, double cz,          // mesh centre (GP)
//                   double lam_iso_outDV_center, double lam_iso_outDV_edge, double lam_aniso_outDV_center, double lam_aniso_outDV_edge,
//                   double lambda_iso_center_inDV, double lambda_iso_edge_inDV, double lambda_aniso_center_inDV, double lambda_aniso_edge_inDV,
//                   double disc_radius, double *rho,   // Strain tensor field parameters from System.h generalParams +disc Radius 
//                   double tFrac,
//                   double *lam_rr,
//                   double *lam_pp,
//                   double *lam_ss, 
//                   const CVec3* e_R, const CVec3* e_phi, const CVec3* e_h, 
//                   Mat_3x3* lam_alpha, const int *layerFlag, const int* DVflag)
//{
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (tid >= N) return;
//    
//    //if (layerFlag == 0 || layerflag == 1) return;
//    
//    /*----  axi-symmetric basis vectors (flat sheet approximation)  ----*/
//    double dx = x[tid] - cx;
//    double dy = y[tid] - cy;
//    double dz = z[tid] - cz;
//    
//    double r = sqrt(dx*dx + dy*dy) + 1e-14;
//     
//    rho[tid] = r/disc_radius;
//                
//    (void)r; // not used here but left for consistency / extensions
//    
//    bool inDV = DVflag[tid];
//    double lamIso, lamAni;
//    
//    if(inDV){
//        lamIso = (lambda_iso_center_inDV   + (lambda_iso_edge_inDV   - lambda_iso_center_inDV)*(rho[tid]*rho[tid]));
//        lamAni = (lambda_aniso_center_inDV + (lambda_aniso_edge_inDV - lambda_aniso_center_inDV)*(rho[tid]*rho[tid]));
//    }else{
//        lamIso = (lam_iso_outDV_center  + (lam_iso_outDV_edge  - lam_iso_outDV_center)*(rho[tid]*rho[tid]));
//        lamAni = (lam_aniso_outDV_center+ (lam_aniso_outDV_edge- lam_aniso_outDV_center)*(rho[tid]*rho[tid]));
//    }
//
//  
//  //    /* linear schedule  ?(t) = I + e t                                          */
//  //    double lamIso = (lam_iso_outDV_center + (lam_iso_outDV_edge - lam_iso_outDV_center)*(rho[tid]*rho[tid]));// * tFrac;   // radial   (_rr)
//  //    double lamAni = (lam_aniso_outDV_center + (lam_aniso_outDV_edge - lam_aniso_outDV_center)*(rho[tid]*rho[tid])); //* tFrac;   // circumf. (_ff)
//    
//    //bool isBasal = (layerFlag[tid] < 0);
//    //if (isBasal) lamAni = 1.0/lamAni;
//    
//    lam_rr[tid] = (lamIso*lamAni);
//    lam_pp[tid] = (lamIso/lamAni);
//    lam_ss[tid] = 1.0;                     // no through-thickness growth
//    
//    
//    // Tensor at v_i
//    
//    Mat_3x3 L = Mat_3x3{
//        CVec3(0,0,0), CVec3(0,0,0), CVec3(0,0,0),
//    };
//    
//    axpy(lam_rr[tid], outer(e_R [tid], e_R [tid]), L);
//    axpy(lam_pp[tid], outer(e_phi[tid], e_phi[tid]), L);
//    axpy(lam_ss[tid], outer(e_h  [tid], e_h  [tid]), L);
//
//    lam_alpha[tid] = L;   // store (row-major) 3×3 tensor
//    
//}
//
//
//
///* ===== Full projection lambda:e_e ============================ */
//__global__
//void k_edgeRestProj(int    E,
//                    const int    *e2n1,  const int *e2n2,
//                    const double *x,     const double *y,   const double *z,
//                    const Mat_3x3 *lam_alpha,          // NEW  ? one ? per vertex
//                    double *L0,                  // ? original, still here
//                    double *Lstar,
//                    const int *edgeLayerFlags, int targetLayer)
//{
//    int eid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (eid >= E) return;
//    
//    // skip edges if not in the desired layer
//    if(edgeLayerFlags[eid]==0 ) return;
//
//    int a = e2n1[eid]; 
//    int b = e2n2[eid];
//    
//    CVec3 dX = CVec3(x[a] - x[b], y[a] - y[b], z[a] - z[b]);
//    
//    // store initial stretch in edge_initial_length
//    L0[eid] = norm3(dX);
//    
//    Mat_3x3 La = lam_alpha[a];
//    Mat_3x3 Lb = lam_alpha[b];
//   Mat_3x3 Lp;
//thrust::get<0>(Lp) = CVec3(
//    0.5*( c<0>( thrust::get<0>(La) ) + c<0>( thrust::get<0>(Lb) )),
//    0.5*( c<1>( thrust::get<0>(La) ) + c<1>( thrust::get<0>(Lb) )),
//    0.5*( c<2>( thrust::get<0>(La) ) + c<2>( thrust::get<0>(Lb) )) );
//
//thrust::get<1>(Lp) = CVec3(
//    0.5*( c<0>( thrust::get<1>(La) ) + c<0>( thrust::get<1>(Lb) )),
//    0.5*( c<1>( thrust::get<1>(La) ) + c<1>( thrust::get<1>(Lb) )),
//    0.5*( c<2>( thrust::get<1>(La) ) + c<2>( thrust::get<1>(Lb) )) );
//
//thrust::get<2>(Lp) = CVec3(
//    0.5*( c<0>( thrust::get<2>(La) ) + c<0>( thrust::get<2>(Lb) )),
//    0.5*( c<1>( thrust::get<2>(La) ) + c<1>( thrust::get<2>(Lb) )),
//    0.5*( c<2>( thrust::get<2>(La) ) + c<2>( thrust::get<2>(Lb) )) );
//
//     /* ---- stretched edge vector ----------------------------- */
//    CVec3 dX_stretch = matVec(Lp, dX);
//
//    /* ---- new rest length ----------------------------------- */
//    Lstar[eid] = norm3(dX_stretch);
//
////    /* edge direction */
////    double ex = x[b]-x[a], ey = y[b]-y[a], ez = z[b]-z[a];
////    double len = sqrt(ex*ex + ey*ey + ez*ez) + 1e-14;
////    ex/=len;  ey/=len;  ez/=len;
////
////    /* average ? diagonal */
////    double lamR = 0.5*(lam_rr[a]+lam_rr[b]);
////    double lamP = 0.5*(lam_pp[a]+lam_pp[b]);
////    double lamS = 0.5*(lam_ss[a]+lam_ss[b]);
////
////
////    double erx = e_R_x[a], ery = e_R_y[a], erz = e_R_z[a];
////    double epx = e_phi_x[a], epy = e_phi_y[a], epz = e_phi_z[a];
////    double ehx = e_h_x[a], ehy = e_h_y[a], ehz = e_h_z[a];
////    
////
//////    /* local polar basis at vertex *a* (flat sheet) */
//////    double rlen = sqrt(x[a]*x[a] + y[a]*y[a]) + 1e-14;
//////    double erx = x[a]/rlen, ery = y[a]/rlen;          // e_r
//////    double epx = -ery,    epy = erx;                  // e_f
////
////    double dR = ex*erx + ey*ery + ez*erz;   // e*e_r
////    double dP = ex*epx + ey*epy + ez*epz;   // e*e_f
////    double dH = ex*ehx + ey*ehy + ez*ehz;   // e*e_h
////
////    double lamProj = lamR*dR*dR + lamP*dP*dP + lamS*dH*dH;
////    Lstar[eid] = lamProj * L0[eid];
//}
//
//
//
////
/////* ----------------------------------------------------------------------------
////   3.  preferred bending angle (anisotropic growth, Eq?25)                    */
////__global__
////void k_prefAngles(int T,
////                  const int    *t2e1, const int *t2e2, const int *t2e3,
////                  const int    *e2n1, const int *e2n2,
////                  const double *lam_rr, const double *lam_pp,
////                  const double *edgeLen,
////                  double thickness,
////                  double* theta0)
////{
////    //--------------------
////    // 1. collect vertices
////    //--------------------
////    int tid = blockIdx.x * blockDim.x + threadIdx.x;
////    if (tid >= T) return;
////    
////    int e1 = t2e1[tid], e2 = t2e2[tid], e3 = t2e3[tid];
////
////    int v1 = e2n1[e1];
////    int v2 = e2n2[e1];
////    /* pick the third vertex: it must be the one that is NOT v1 or v2 */
////    int v3 = (e2n1[e2]!=v1 && e2n1[e2]!=v2) ? e2n1[e2] :
////             (e2n2[e2]!=v1 && e2n2[e2]!=v2) ? e2n2[e2] :
////             (e2n1[e3]!=v1 && e2n1[e3]!=v2) ? e2n1[e3] : e2n2[e3];
////
////    //-----------------------------
////    // 2. average ? over the triple
////    //-----------------------------
////    double lamR = (lam_rr[v1] + lam_rr[v2] + lam_rr[v3]) / 3.0;
////    double lamP = (lam_pp[v1] + lam_pp[v2] + lam_pp[v3]) / 3.0;
////
////    //--------------------------------------------
////    // 3. spontaneous curvature   C0  (half-plate)
////    //--------------------------------------------
////    double C0 = 0.5 * (lamP - lamR) / thickness;   // ? ½ factor
////
////    /* mean surrounding edge length (more robust than triangle edge-specific) */
////    double Lm = (edgeLen[e1] + edgeLen[e2] + edgeLen[e3]) / 3.0;
////
////    theta0[tid] = C0 * Lm;                     // dihedral offset
////}
//
///* ----------------------------------------------------------------------------
//   ====   PUBLIC WRAPPERS   ================================================== */
//namespace StrainTensorGPU {
//
//void buildVertexLambda(GeneralParams& gp,
//                       CoordInfoVecs& coord,
//                       LambdaField&         field,
//                       double               tFrac)
//{
//    if (field.lam_rr.size() != coord.nodeLocX.size())
//        field.resize(coord.nodeLocX.size());
//    if (gp.rho.size() != coord.nodeLocX.size())
//        gp.rho.resize(coord.nodeLocX.size());
//
//    int N = static_cast<int>(coord.nodeLocX.size());
//    dim3 grid((N + BLOCK_SZ - 1) / BLOCK_SZ);
//
//    // allocate once
//    if(gp.nodes_in_DV.size()!=N) gp.nodes_in_DV.resize(N);
//    
//    // --- mark DV -----------------------------------------------
//    double frac = 0.05;                     // 5 % of apical radius
//    k_markDVstripe<<<grid,BLOCK_SZ>>>( N,
//        thrust::raw_pointer_cast(coord.nodeLocX.data()),
//        thrust::raw_pointer_cast(coord.nodeLocY.data()),
//        gp.centerX, gp.centerY,
//        gp.apical_Rx,
//        frac,
//        thrust::raw_pointer_cast(gp.nodes_in_upperhem.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_DV.data()));
//        
//        cudaDeviceSynchronize();
//
//    // build basis vectors
//    k_buildBasis<<<grid,BLOCK_SZ>>>(
//        N,
//        thrust::raw_pointer_cast(coord.nodeLocX.data()), thrust::raw_pointer_cast(coord.nodeLocY.data()), thrust::raw_pointer_cast(coord.nodeLocZ.data()),
//        gp.c_dx, gp.c_dy, gp.c_dz,
//        gp.centerX, gp.centerY, gp.centerZ,
//        thrust::raw_pointer_cast(field.e_h.data()), thrust::raw_pointer_cast(field.e_R.data()), thrust::raw_pointer_cast(field.e_phi.data()));
//        
//        cudaDeviceSynchronize();
//    
//    // build lambda field 
//    k_buildLambda<<<grid,BLOCK_SZ>>>(
//        N,
//        thrust::raw_pointer_cast(coord.nodeLocX.data()), thrust::raw_pointer_cast(coord.nodeLocY.data()), thrust::raw_pointer_cast(coord.nodeLocZ.data()),
//        gp.centerX, gp.centerY, gp.centerZ,
//        gp.lambda_iso_center_outDV, gp.lambda_iso_edge_outDV,
//        gp.lambda_aniso_center_outDV, gp.lambda_aniso_edge_outDV,
//        gp.lambda_iso_center_inDV, gp.lambda_iso_edge_inDV,
//        gp.lambda_aniso_center_inDV, gp.lambda_aniso_edge_inDV,
//        gp.disc_radius,
//        thrust::raw_pointer_cast(gp.rho.data()),
//        tFrac,
//        thrust::raw_pointer_cast(field.lam_rr.data()), thrust::raw_pointer_cast(field.lam_pp.data()), thrust::raw_pointer_cast(field.lam_ss.data()),
//        thrust::raw_pointer_cast(field.e_R.data()), thrust::raw_pointer_cast(field.e_phi.data()), thrust::raw_pointer_cast(field.e_h.data()),
//        thrust::raw_pointer_cast(field.lam_alpha.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_upperhem.data()),
//        thrust::raw_pointer_cast(gp.nodes_in_DV.data()));
//    cudaDeviceSynchronize();
//}
//
///* ------------------------------------------------------------------------- */
//void updateEdgeRestLengths(CoordInfoVecs&  coord,
//                           GeneralParams& gp,
//                           const LambdaField&    field,
//                           LinearSpringInfoVecs& lsInfo, int targetLayer)
//{
//    int E = static_cast<int>(coord.num_edges);
//    dim3 grid((E + BLOCK_SZ - 1) / BLOCK_SZ); 
//
//
//  /* full lambda:e_e projection */
//    k_edgeRestProj<<<grid,BLOCK_SZ>>>(
//        E,
//        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()), thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
//        thrust::raw_pointer_cast(coord.nodeLocX.data()),      thrust::raw_pointer_cast(coord.nodeLocY.data()),      thrust::raw_pointer_cast(coord.nodeLocZ.data()),
//        thrust::raw_pointer_cast(field.lam_alpha.data()),        // NEW — whole tensor array
//        thrust::raw_pointer_cast(lsInfo.edge_initial_length.data()),
//        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
//        thrust::raw_pointer_cast(gp.edges_in_upperhem.data()),
//        targetLayer );
//
//    cudaDeviceSynchronize();
//}
////
//////==============================================================
//////  StrainTensorGPU::updatePreferredAngles (in .cu file)
//////==============================================================
////void updatePreferredAngles(BendingTriangleInfoVecs& btiv,
////                           const CoordInfoVecs& coord,
////                           const LambdaField& field,
////                           const GeneralParams& gp,
////                           const LinearSpringInfoVecs& lsInfo)
////{
////    int  T      = coord.triangles2Edges_1.size();
////    int  TPB    = 128;
////    int  nBlk   = (T + TPB - 1) / TPB;
////    btiv.initial_angle.resize(T);
////
////    k_prefAngles<<<nBlk, TPB>>>(
////        T,
////        thrust::raw_pointer_cast(coord.triangles2Edges_1.data()),
////        thrust::raw_pointer_cast(coord.triangles2Edges_2.data()),
////        thrust::raw_pointer_cast(coord.triangles2Edges_3.data()),
////        thrust::raw_pointer_cast(coord.edges2Nodes_1.data()),
////        thrust::raw_pointer_cast(coord.edges2Nodes_2.data()),
////        thrust::raw_pointer_cast(field.lam_rr.data()),
////        thrust::raw_pointer_cast(field.lam_pp.data()),
////        thrust::raw_pointer_cast(lsInfo.edge_final_length.data()),
////        gp.thickness,
////        thrust::raw_pointer_cast(btiv.initial_angle.data()) );
////
////   // CUDA_CHECK_ERROR();
////}
////
//
//} // namespace StrainTensorGPU
//
//
//
