#ifndef VOLUMECOMP_H_
#define VOLUMECOMP_H_

#include "SystemStructures.h"
#include "System.h"

// Compute global enclosed volume for a single 2-layer hexagonal unit
// using the 6-prism, 3-tet decomposition from your notes.
//
// On return:
//   generalParams.current_total_volume      = signed volume (can be +/-)
//   generalParams.true_current_total_volume = |volume|
void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs,
    PrismInfoVecs& prismInfoVecs);

// Optional device-side prism functor (not currently used by the host version).
// This one just returns 6 * |V_prism|, i.e. unsigned contribution.
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

        return ux * cx + uy * cy + uz * cz;
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

        return s1 + s2 + s3; // 6 * |V_prism|
    }
};

#endif // VOLUMECOMP_H_
