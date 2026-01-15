#ifndef STRAINTENSOR_H_
#define STRAINTENSOR_H_

#include "SystemStructures.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

/**
 * LambdaField: Stores per-vertex strain tensor information
 * 
 * This structure holds the decomposed lambda values (?_RR, ?_ff, ?_hh) 
 * and the local basis vectors (e_R, e_f, e_h) for each vertex.
 * 
 * The strain tensor at each vertex is:
 *   ? = ?_RR * (e_R ? e_R) + ?_ff * (e_f ? e_f) + ?_hh * (e_h ? e_h)
 * 
 * Following the Fuhrmann et al. paper:
 *   ?_RR = ?_iso * ?_aniso
 *   ?_ff = ?_iso / ?_aniso  (NOTE: division, not multiplication!)
 *   ?_hh = 1.0 (typically, no strain in thickness direction)
 */
struct LambdaField {
    // Per-vertex lambda components
    thrust::host_vector<double> lambda_RR;      // Radial stretch
    thrust::host_vector<double> lambda_phiphi;  // Circumferential stretch
    thrust::host_vector<double> lambda_hh;      // Height/thickness stretch
    
    // Per-vertex pathlength (normalized to [0,1])
    thrust::host_vector<double> pathlength_scaled;
    
    // Per-vertex basis vectors: e_R (radial direction in surface)
    thrust::host_vector<double> e_R_x;
    thrust::host_vector<double> e_R_y;
    thrust::host_vector<double> e_R_z;
    
    // Per-vertex basis vectors: e_phi (circumferential direction in surface)
    thrust::host_vector<double> e_phi_x;
    thrust::host_vector<double> e_phi_y;
    thrust::host_vector<double> e_phi_z;
    
    // Per-vertex basis vectors: e_h (surface normal / height direction)
    thrust::host_vector<double> e_h_x;
    thrust::host_vector<double> e_h_y;
    thrust::host_vector<double> e_h_z;
    
    // Resize all vectors to N vertices
    void resize(int N) {
        lambda_RR.resize(N);
        lambda_phiphi.resize(N);
        lambda_hh.resize(N);
        pathlength_scaled.resize(N);
        e_R_x.resize(N); e_R_y.resize(N); e_R_z.resize(N);
        e_phi_x.resize(N); e_phi_y.resize(N); e_phi_z.resize(N);
        e_h_x.resize(N); e_h_y.resize(N); e_h_z.resize(N);
    }
};

/**
 * StrainTensorGPU namespace: Functions for computing and applying strain
 */
namespace StrainTensorGPU {

    /**
     * computeBasisVectorsAndPathlength
     * 
     * Computes the local basis vectors (e_R, e_phi, e_h) and pathlength_scaled
     * for each vertex, following the Fuhrmann et al. paper methodology.
     * 
     * Key features:
     * - Different center points for DV vs outDV regions
     * - Pathlength normalized separately within each region
     * - Basis vectors computed from surface geometry
     * 
     * @param generalParams  Contains nodes_in_DV and geometry parameters
     * @param coordInfoVecs  Contains node positions
     * @param hostSetInfoVecs Contains nodes_in_DV array
     * @param lambda         Output: LambdaField with basis vectors and pathlength
     * @param theta_DV       Angular width of the DV region (radians)
     * @param R              Radius of the spherical cap
     */
    void computeBasisVectorsAndPathlength(
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LambdaField& lambda,
        double theta_DV,
        double R = 1.0);

    /**
     * buildVertexLambda
     * 
     * Computes the lambda values (?_RR, ?_ff, ?_hh) for each vertex
     * based on the pathlength and region (DV vs outDV).
     * 
     * Uses linear interpolation: ?(p) = center + (edge - center) * p
     * where p = pathlength_scaled ? [0, 1]
     * 
     * @param generalParams  Contains lambda center/edge values for each region
     * @param coordInfoVecs  Contains node positions
     * @param lambda         LambdaField with pathlength (input) and lambda values (output)
     * @param hostSetInfoVecs Contains nodes_in_DV array
     * @param frac           Fraction of strain to apply (for quasi-static loading)
     */
    void buildVertexLambda(
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LambdaField& lambda,
        double frac = 1.0);

    /**
     * updateEdgeRestLengths
     * 
     * Computes the target rest lengths for each edge based on the strain tensor.
     * 
     * For each edge:
     *   1. Get the spring vector v = (x2-x1, y2-y1, z2-z1)
     *   2. Average the lambda tensors at the two endpoints
     *   3. Transform the spring vector: v' = ?_avg · v
     *   4. New rest length = |v'|
     * 
     * @param coordInfoVecs  Contains node positions and edge connectivity
     * @param generalParams  Contains layer flags
     * @param lambda         LambdaField with lambda values and basis vectors
     * @param linearSpringInfoVecs  Contains edge_initial_length, edge_final_length
     * @param hostSetInfoVecs Contains edge layer information
     * @param layerflag      Which layer to apply strain to (-1 = all, 0 = basal, etc.)
     */
    void updateEdgeRestLengths(
        CoordInfoVecs& coordInfoVecs,
        GeneralParams& generalParams,
        LambdaField& lambda,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        int layerflag = -1);

} // namespace StrainTensorGPU

#endif /* STRAINTENSOR_H_ */

