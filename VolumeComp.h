#ifndef VOLUMECOMP_H_
#define VOLUMECOMP_H_

#include "SystemStructures.h"
#include "System.h"
#include <cmath>   // fabs

// Compute global enclosed volume by summing over an explicit list of prisms.
// Each prism is given by 6 node indices (P1..P6) stored in prismInfoVecs.
//
// Uses the SAME 3-tetrahedron decomposition per prism as in VolumeComp.cu:
//
//   s1 = 6V_tet(P2,P3,P4,P1)
//   s2 = 6V_tet(P2,P4,P6,P1)
//   s3 = 6V_tet(P4,P5,P6,P1)
//
// Stores:
//   generalParams.current_total_volume      = signed volume (can be +/-)
//   generalParams.true_current_total_volume = sum of absolute tet volumes (physical magnitude)
//
// Notes:
// - This function expects prismInfoVecs.num_prisms and prismInfoVecs.P1..P6 to be valid.
// - Prism construction is done elsewhere (e.g., BuildPrismsFromLayerFlags).
void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs,
    PrismInfoVecs& prismInfoVecs);

// Optional device-side functor for a prism list (not used by your host ComputeVolume).
// Returns 6 * |V_prism| for prism p, using the same 3-tet decomposition.
struct PrismTetFunctor {
    const double* X; const double* Y; const double* Z;
    const int* P1; const int* P2; const int* P3;
    const int* P4; const int* P5; const int* P6;

    __host__ __device__
    PrismTetFunctor(const double* _X, const double* _Y, const double* _Z,
                    const int* _P1, const int* _P2, const int* _P3,
                    const int* _P4, const int* _P5, const int* _P6)
        : X(_X), Y(_Y), Z(_Z),
          P1(_P1), P2(_P2), P3(_P3),
          P4(_P4), P5(_P5), P6(_P6) {}

    __device__ inline double sixV_tet(int i, int j, int k, int l) const {
        const double ux = X[i] - X[l];
        const double uy = Y[i] - Y[l];
        const double uz = Z[i] - Z[l];

        const double vx = X[j] - X[l];
        const double vy = Y[j] - Y[l];
        const double vz = Z[j] - Z[l];

        const double wx = X[k] - X[l];
        const double wy = Y[k] - Y[l];
        const double wz = Z[k] - Z[l];

        const double cx = vy * wz - vz * wy;
        const double cy = vz * wx - vx * wz;
        const double cz = vx * wy - vy * wx;

        return ux * cx + uy * cy + uz * cz; // 6 * V_signed
    }

    __device__ double operator()(int p) const {
        const int p1 = P1[p];
        const int p2 = P2[p];
        const int p3 = P3[p];
        const int p4 = P4[p];
        const int p5 = P5[p];
        const int p6 = P6[p];

        const double s1 = fabs(sixV_tet(p2, p3, p4, p1));
        const double s2 = fabs(sixV_tet(p2, p4, p6, p1));
        const double s3 = fabs(sixV_tet(p4, p5, p6, p1)); 

        return s1 + s2 + s3; // 6 * |V_prism| (sum of unsigned tet contributions)
    }
};

#endif // VOLUMECOMP_H_
