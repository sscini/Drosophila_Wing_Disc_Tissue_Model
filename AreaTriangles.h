#ifndef AREATRIANGLES_H_
#define AREATRIANGLES_H_ 

#include "SystemStructures.h"

void ComputeAreaTriangleSprings(  
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);


struct AreaSpringFunctor {
    // Parameters
    int SCALE_TYPE;
    bool nonuniform_wall_weakening_area;
    double maxSpringScaler_area;
    double scaling_pow;
    double gausssigma;
    double hilleqnconst;
    double hilleqnpow;
    double* scaling_per_edge;
	double area_0;
    double spring_constant;
    double spring_constant_weak;
    int* triangles_in_upperhem;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    // Output arrays
    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ AreaSpringFunctor(
        int& _SCALE_TYPE,
        bool& _nonuniform_wall_weakening_area,
        double& _maxSpringScaler_area,
        double& _scaling_pow,
        double& _gausssigma,
        double& _hilleqnconst,
        double& _hilleqnpow,
        double* _scaling_per_edge,
        double& _area_0,
        double& _spring_constant,
        double& _spring_constant_weak,
        int* _triangles_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr):
        SCALE_TYPE(_SCALE_TYPE),
        nonuniform_wall_weakening_area(_nonuniform_wall_weakening_area),
        maxSpringScaler_area(_maxSpringScaler_area),
        scaling_pow(_scaling_pow),
        gausssigma(_gausssigma),
        hilleqnconst(_hilleqnconst),
        hilleqnpow(_hilleqnpow),
        scaling_per_edge(_scaling_per_edge),
        area_0(_area_0),
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        triangles_in_upperhem(_triangles_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	      // Operator for computing area triangle springs
        //hand in counting iterator and id of triangle
	__device__ double operator()(const Tuuuuuuu &u7) {
        //test placing the ids of the nodes and then get positions. 
        //double scaling_pow = 4.0;
		int counter = thrust::get<0>(u7);
		int place = 3 * counter; // Represents location in write to vector.

        // Check if the triangle ids are valid
        int id_i = thrust::get<1>(u7);
        int id_j = thrust::get<2>(u7);
        int id_k = thrust::get<3>(u7);
        int e_id_i = thrust::get<4>(u7);
        int e_id_j = thrust::get<5>(u7);
        int e_id_k = thrust::get<6>(u7);
        double target_area;
        
        if ((id_i < (INT_MAX-100) && id_i >= 0) && (id_j < (INT_MAX-100) && id_j >= 0) && (id_k < (INT_MAX-100) && id_k >= 0)){   
            double what_spring_constant;
            if (SCALE_TYPE == 0){
             what_spring_constant = spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*gausssigma))*exp(-pow(scaling_per_edge[e_id_i],2.0)/gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*gausssigma))*exp(-pow(scaling_per_edge[e_id_j],2.0)/gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*gausssigma))*exp(-pow(scaling_per_edge[e_id_k],2.0)/gausssigma))))/3.0;
                                            if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant_weak;}
            }
            else if (SCALE_TYPE == 1){
                 what_spring_constant = ((spring_constant_weak*2.0)*pow(scaling_per_edge[e_id_i],scaling_pow) + spring_constant_weak*(1-pow(scaling_per_edge[e_id_i],scaling_pow)) +
                                            (spring_constant_weak*2.0)*pow(scaling_per_edge[e_id_j],scaling_pow) + spring_constant_weak*(1-pow(scaling_per_edge[e_id_j],scaling_pow)) +
                                            (spring_constant_weak*2.0)*pow(scaling_per_edge[e_id_k],scaling_pow) + spring_constant_weak*(1-pow(scaling_per_edge[e_id_k],scaling_pow)))/3.0;
                            }
            else if (SCALE_TYPE == 2){
                what_spring_constant = (spring_constant - (spring_constant - spring_constant_weak)*scaling_per_edge[e_id_i] +
                                            spring_constant - (spring_constant - spring_constant_weak)*scaling_per_edge[e_id_j] +
                                            spring_constant - (spring_constant - spring_constant_weak)*scaling_per_edge[e_id_k])/3.0;
            }
            else if (SCALE_TYPE == 3){        
                if (triangles_in_upperhem[counter] == 1){
                    what_spring_constant = spring_constant_weak;
                    target_area = area_0*1.0;
                }
                else if (triangles_in_upperhem[counter] == 0){
                    what_spring_constant = (spring_constant_weak + spring_constant)/2.0;
                    target_area = area_0*1.0;
                }
                else{
                    what_spring_constant = spring_constant;
                    target_area = area_0*1.0;
                }
            }
            else if (SCALE_TYPE == 4){
                if (nonuniform_wall_weakening_area == true){
                    //double scaling = 0.0;//spring_constant_weak/spring_constant;
                    double spectrum = maxSpringScaler_area*spring_constant - spring_constant_weak;
                    what_spring_constant = (spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[e_id_i], hilleqnpow)))*spectrum) +
                                            spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[e_id_j], hilleqnpow)))*spectrum) +
                                            spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[e_id_k], hilleqnpow)))*spectrum))/3.0;
                    if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant_weak;}
                    target_area = area_0*1.0;
                }
                else{
                    if (triangles_in_upperhem[counter] == 1){// && triangles_in_tip[counter]==1){
                        what_spring_constant = spring_constant_weak;
                        target_area = area_0*1.0;
                    }
                    else if (triangles_in_upperhem[counter] == 0){// && triangles_in_tip[counter]!=1){
                        //what_spring_constant = (spring_constant_weak*2.0);
                        what_spring_constant = (spring_constant_weak + spring_constant)/2.0;
                        target_area = area_0*1.0;
                    }
                    else{
                        what_spring_constant = spring_constant;
                        target_area = area_0*1.0;
                    }
                }
		    }

        
            // Calculate positions of the nodes
            CVec3 ri = thrust::make_tuple(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
            CVec3 rj = thrust::make_tuple(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
            CVec3 rk = thrust::make_tuple(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);

            CVec3 rkj = CVec3_subtract(rk, rj);
            CVec3 rij = CVec3_subtract(ri, rj);

            double area_current = sqrt( CVec3_dot( CVec3_cross(rkj, rij), CVec3_cross(rkj, rij) ) )/2.0;

            // Compute derivative wrt to area
            CVec3 A = CVec3_cross(rkj, rij);//rkj must come first
            double A1 = thrust::get<0>(A);
            double A2 = thrust::get<1>(A);
            double A3 = thrust::get<2>(A);
            
            CVec3 A1Rj = thrust::make_tuple(
                0.0, -thrust::get<2>(rij) + thrust::get<2>(rkj), -thrust::get<1>(rkj) + thrust::get<1>(rij));
                //[0, -Rij(3)+Rkj(3), -Rkj(2)+Rij(2)];
            CVec3 A2Rj = thrust::make_tuple(    
                thrust::get<2>(rij) - thrust::get<2>(rkj), 0.0, thrust::get<0>(rkj) - thrust::get<0>(rij));
                //[Rij(3)-Rkj(3), 0, Rkj(1)-Rij(1)];
            CVec3 A3Rj = thrust::make_tuple(
                -thrust::get<1>(rij) + thrust::get<1>(rkj), -thrust::get<0>(rkj) + thrust::get<0>(rij), 0.0);
                //[-Rij(2)+Rkj(2), -Rkj(1)+Rij(1), 0];

            CVec3 A1Rk = thrust::make_tuple(
                0.0, thrust::get<2>(rij), -thrust::get<1>(rij));
                //[0, Rij(3), -Rij(2)];
            //Derivative of A1 with respect to [Rkx, Rky, Rkz].
            CVec3 A2Rk = thrust::make_tuple(
                -thrust::get<2>(rij), 0.0, thrust::get<0>(rij));
                //[-Rij(3), 0, Rij(1)];
            CVec3 A3Rk = thrust::make_tuple(
                thrust::get<1>(rij), -thrust::get<0>(rij), 0.0);
                //[Rij(2), -Rij(1), 0];
            CVec3 A1Ri = thrust::make_tuple(
                0.0, -thrust::get<2>(rkj), thrust::get<1>(rkj));
                //[0, -Rkj(3), Rkj(2)];
            CVec3 A2Ri = thrust::make_tuple(
                thrust::get<2>(rkj), 0.0, -thrust::get<0>(rkj));
                //[Rkj(3), 0, -Rkj(1)];
            CVec3 A3Ri = thrust::make_tuple(
                -thrust::get<1>(rkj), thrust::get<0>(rkj), 0.0);
                //[-Rkj(2), Rkj(1), 0];
        
            double magnitude = -((what_spring_constant) * (area_current - area_0)/area_0) / (2.0 * area_current) ;
            CVec3 rj_force = CVec3_scalermult(magnitude, CVec3_plus( 
                CVec3_scalermult(A1, A1Rj), CVec3_scalermult(A2, A2Rj), CVec3_scalermult(A3, A3Rj)));
                // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(1) + A2*A2Rj(1) + A3*A3Rj(1))/(2*AREA(i));
                // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(2) + A2*A2Rj(2) + A3*A3Rj(2))/(2*AREA(i));
                // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(3) + A2*A2Rj(3) + A3*A3Rj(3))/(2*AREA(i));
                
            CVec3 rk_force = CVec3_scalermult(magnitude, CVec3_plus( 
                CVec3_scalermult(A1, A1Rk), CVec3_scalermult(A2, A2Rk), CVec3_scalermult(A3, A3Rk)));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(1) + A2*A2Rk(1) + A3*A3Rk(1))/(2*AREA(i));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(2) + A2*A2Rk(2) + A3*A3Rk(2))/(2*AREA(i));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(3) + A2*A2Rk(3) + A3*A3Rk(3))/(2*AREA(i));

            CVec3 ri_force = CVec3_scalermult(magnitude, CVec3_plus( 
                CVec3_scalermult(A1, A1Ri), CVec3_scalermult(A2, A2Ri), CVec3_scalermult(A3, A3Ri)));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(1) + A2*A2Ri(1) + A3*A3Ri(1))/(2*AREA(i));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(2) + A2*A2Ri(2) + A3*A3Ri(2))/(2*AREA(i));
                //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(3) + A2*A2Ri(3) + A3*A3Ri(3))/(2*AREA(i));
            
            // Store the forces and ids in the output arrays
            idKey[place] = id_i;
            forceXAddr[place] = thrust::get<0>( ri_force );
            forceYAddr[place] = thrust::get<1>( ri_force );
            forceZAddr[place] = thrust::get<2>( ri_force );

            idKey[place+1] = id_j;
            forceXAddr[place+1] = thrust::get<0>( rj_force );
            forceYAddr[place+1] = thrust::get<1>( rj_force );
            forceZAddr[place+1] = thrust::get<2>( rj_force );
            
            idKey[place+2] = id_k;
            forceXAddr[place+2] = thrust::get<0>( rk_force );
            forceYAddr[place+2] = thrust::get<1>( rk_force );
            forceZAddr[place+2] = thrust::get<2>( rk_force );

            // Calculate energy
            double energy =  (what_spring_constant/(2.0)) * (area_current - target_area) * (area_current - target_area) / target_area;
            return energy;
        }
        else{
            // If the triangle ids are invalid, return 0 energy
            double energy = 0.0;
            return energy;
        }

    };
};
#endif //AREATRIANGLES_H_
