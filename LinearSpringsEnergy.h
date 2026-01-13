#ifndef LINEARSPRINGSENERGY_H_
#define LINEARSPRINGSENERGY_H_ 

#include "SystemStructures.h"
#include <cstdio>

double ComputeLinearSpringsEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct LinearSpringEnergyFunctor {
    
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    
	__host__ __device__ LinearSpringEnergyFunctor(
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        int* _idKey
        ) :
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuud &u3d) {
        		
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u3d);
    		int place = 2 * counter;//represents location in write to vector.

        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);
        double length_zero = thrust::get<3>(u3d);

        // compute forces.
        double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
        double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
        double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
   

        double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                                    (yLoc_LR) * (yLoc_LR)  + 
                                    (zLoc_LR) * (zLoc_LR) );

        double energy;

        if (length_current >= length_zero){		
			
				idKey[place] = edgeL;
                energy = (spring_constant/2.0) * (length_current - length_zero) * (length_current - length_zero);
		}
		else{	
			
				idKey[place] = edgeL;
                energy = 0.0;

		}
       
    //std::cout << "Current Length (linear springs energy) = "<<length_current<<". Zero length (same file) = " << length_zero<< std::endl;
        printf("current_len = %g, zero_len = %g\n", length_current, length_zero);
        //double energy = (spring_constant/2.0) * (length_current - length_zero) * (length_current - length_zero);
        return energy;


    }
};

#endif