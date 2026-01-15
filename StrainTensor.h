#ifndef STRAINTENSOR_H_
#define STRAINTENSOR_H_

#include <vector>
#include "SystemStructures.h"
#include "System.h"

/**
 * LambdaField: Stores per-vertex strain tensor information on HOST
 * 
 * This structure holds the decomposed lambda values and basis vectors.
 * All vectors are std::vector for CPU-side processing.
 */
struct LambdaField {
    // Per-vertex lambda components
    std::vector<double> lambda_RR;      // Radial stretch
    std::vector<double> lambda_phiphi;  // Circumferential stretch  
    std::vector<double> lambda_hh;      // Height/thickness stretch
    
    // Per-vertex pathlength (normalized to [0,1])
    std::vector<double> pathlength_scaled;
    
    // Per-vertex DV region flag
    std::vector<int> nodes_in_DV;
    
    // Per-vertex basis vectors: e_R (radial direction in surface)
    std::vector<double> e_R_x;
    std::vector<double> e_R_y;
    std::vector<double> e_R_z;
    
    // Per-vertex basis vectors: e_phi (circumferential direction in surface)
    std::vector<double> e_phi_x;
    std::vector<double> e_phi_y;
    std::vector<double> e_phi_z;
    
    // Per-vertex basis vectors: e_h (surface normal / height direction)
    std::vector<double> e_h_x;
    std::vector<double> e_h_y;
    std::vector<double> e_h_z;
    
    // Resize all vectors to N vertices
    void resize(int N) {
        lambda_RR.resize(N, 1.0);
        lambda_phiphi.resize(N, 1.0);
        lambda_hh.resize(N, 1.0);
        pathlength_scaled.resize(N, 0.0);
        nodes_in_DV.resize(N, 0);
        e_R_x.resize(N, 0.0); e_R_y.resize(N, 0.0); e_R_z.resize(N, 0.0);
        e_phi_x.resize(N, 0.0); e_phi_y.resize(N, 0.0); e_phi_z.resize(N, 0.0);
        e_h_x.resize(N, 0.0); e_h_y.resize(N, 0.0); e_h_z.resize(N, 0.0);
    }
};

// Forward declarations - actual structs defined in SystemStructures.h
struct GeneralParams;
struct CoordInfoVecs;
struct LinearSpringInfoVecs;

/**
 * StrainTensorGPU namespace: Functions for computing and applying strain
 */
namespace StrainTensorGPU {

    /**
     * computeBasisVectorsAndPathlength
     * 
     * Computes basis vectors and pathlength for each vertex.
     * Copies data from device to host internally.
     */
    void computeBasisVectorsAndPathlength(
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LambdaField& lambda,
        double theta_DV = 0.1931,
        double R = 1.0);

    /**
     * buildVertexLambda
     * 
     * Computes lambda values for each vertex based on pathlength and region.
     */
    void buildVertexLambda(
        GeneralParams& generalParams,
        LambdaField& lambda,
        double frac = 1.0);

    /**
     * updateEdgeRestLengths
     * 
     * Computes target rest lengths for each edge based on strain tensor.
     */
    void updateEdgeRestLengths(
        CoordInfoVecs& coordInfoVecs,
        GeneralParams& generalParams,
        LambdaField& lambda,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        int layerflag = -1);

    /**
     * Helper: Get converted lambda coefficients for a stage
     */
    void getLambdaCoeffsForStage(
        int stage,
        double& iso_center_outDV, double& iso_edge_outDV,
        double& aniso_center_outDV, double& aniso_edge_outDV,
        double& iso_center_inDV, double& iso_edge_inDV,
        double& aniso_center_inDV, double& aniso_edge_inDV);

} // namespace StrainTensorGPU

#endif /* STRAINTENSOR_H_ */