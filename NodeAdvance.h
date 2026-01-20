#ifndef NODEADVANCE_H_
#define NODEADVANCE_H_

#include "SystemStructures.h"

void AdvancePositions(
	CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
	DomainParams& domainParams);


struct SaxpyFunctorPrimary {

    double dt;
    double mass;
    int maxNodeCount;
    double domainLengthX;
    double domainLengthY;
    double domainLengthZ;
    
    __host__ __device__
    SaxpyFunctorPrimary(
    
        double& _dt,
        double& _mass, 
        int& _maxNodeCount,
        double& _domainLengthX,
        double& _domainLengthY,
        double& _domainLengthZ):
        
        dt(_dt),
        mass(_mass),
        maxNodeCount(_maxNodeCount),
        domainLengthX(_domainLengthX),
        domainLengthY(_domainLengthY),
        domainLengthZ(_domainLengthZ) {}
        
    __device__
    CVec3 operator()(const BoolCVec3 &p3, const CVec3 &f3) {
        bool isFixed = thrust::get<0>(p3); // true if fixed, false if movable

        // Get current position
        double xLoc = thrust::get<1>(p3);
        double yLoc = thrust::get<2>(p3);
        double zLoc = thrust::get<3>(p3);
        
        // Get forces
        double fx = thrust::get<0>(f3);
        double fy = thrust::get<1>(f3);
        double fz = thrust::get<2>(f3);
        
        // ==================== SAFEGUARD: Check for NaN/Inf forces ====================
        if (isnan(fx) || isnan(fy) || isnan(fz) ||
            isinf(fx) || isinf(fy) || isinf(fz)) {
            // Return current position unchanged if force is invalid
            return thrust::make_tuple(xLoc, yLoc, zLoc);
        }

        if (isFixed == false) {
            // Compute displacement
            double d_x = dt * fx;
            double d_y = dt * fy;
            double d_z = dt * fz;
            
            // ==================== SAFEGUARD: Clamp maximum displacement ====================
            // This prevents nodes from moving too far in a single timestep,
            // which can cause mesh inversion and numerical instability.
            // Adjust max_disp based on your typical edge length (should be ~0.1 to 0.5 of min edge length)
            const double max_disp = 0.5;  // Maximum displacement per timestep
            
            double disp_mag = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
            if (disp_mag > max_disp) {
                double scale = max_disp / disp_mag;
                d_x *= scale;
                d_y *= scale;
                d_z *= scale;
            }
            
            double xLocNew = xLoc + d_x;
            double yLocNew = yLoc + d_y;
            double zLocNew = zLoc + d_z;
            
            return thrust::make_tuple(xLocNew, yLocNew, zLocNew);
        }
        else {
            // Fixed nodes: very small movement (essentially fixed but with tiny compliance)
            // Using multiplication instead of addition (your original had a bug with +)
            double scale_fixed = 0.00001;
            double xLocNew = xLoc + scale_fixed * fx;
            double yLocNew = yLoc + scale_fixed * fy;
            double zLocNew = zLoc + scale_fixed * fz;
            
            return thrust::make_tuple(xLocNew, yLocNew, zLocNew);
        }
    }
};

#endif /*NODEADVANCE_H_*/
