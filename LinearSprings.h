#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"
#include <stdio.h>
// Declare the function ComputeLinearSprings, which is responsible for calculating linear spring forces. It takes references to various data structures related to the simulation as arguments.
void ComputeLinearSprings(
  
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);

// Define a functor named LinearSpringFunctor, which will be used to calculate linear spring forces for each edge (linear spring) between two connected nodes.    
struct LinearSpringFunctor {
    // Member variables representing parameters required by the functor.
    int SCALE_TYPE;
    bool nonuniform_wall_weakening_linear;
    double maxSpringScaler_linear;
    double scaling_pow;
    double gausssigma;
    double hilleqnconst;
    double hilleqnpow;
    double* scaling_per_edge;
    double spring_constant;
    double spring_constant_weak;
    double spring_constant_vertical;
   // double length_zero;
    //double length_zero_growth;
    double* rest_length;
    int* edges_in_upperhem;
    //int* boundaries_in_upperhem;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
// Constructor for the LinearSpringFunctor, which initializes member variables with the provided arguments.    
	__host__ __device__ LinearSpringFunctor(
        int& _SCALE_TYPE,
        bool& _nonuniform_wall_weakening_linear,
        double& _maxSpringScaler_linear,
        double& _scaling_pow,
        double& _gausssigma,
        double& _hilleqnconst,
        double& _hilleqnpow,
        double* _scaling_per_edge,
        double& _spring_constant,
        double& _spring_constant_weak,
        double& _spring_constant_vertical,
       // double& _length_zero,
       // double& _length_zero_growth,
        double* _rest_length,
        int* _edges_in_upperhem,
       // int* _boundaries_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :
        SCALE_TYPE(_SCALE_TYPE),
        nonuniform_wall_weakening_linear(_nonuniform_wall_weakening_linear),
        maxSpringScaler_linear(_maxSpringScaler_linear),
        scaling_pow(_scaling_pow),
        gausssigma(_gausssigma),
        hilleqnconst(_hilleqnconst),
        hilleqnpow(_hilleqnpow),
        scaling_per_edge(_scaling_per_edge),
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        spring_constant_vertical(_spring_constant_vertical),
        //length_zero(_length_zero),
       // length_zero_growth(_length_zero_growth),
        rest_length(_rest_length),
        edges_in_upperhem(_edges_in_upperhem),
       // boundaries_in_upperhem(_boundaries_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
// Functor operator for the LinearSpringFunctor, which calculates the linear spring forces for each edge.
	__device__ double operator()(const Tuuu &u3d) {
        		
             //   double scaling_pow = 4.0;
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u3d);
		int place = 2 * counter;//represents location in write to vector.

        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

//        if (edgeL != INT_MAX && edgeL >= 0 && edgeR != INT_MAX && edgeR >= 0){
//            if (edgeL != INT_MAX && edgeL >= 0 edgeR != INT_MAX && edgeR>= 0){
//                double length_zero = rest_length[counter];    
//            }
      if (edgeL != INT_MAX && edgeL >= 0 && edgeR != INT_MAX && edgeR >= 0){
          double length_zero = rest_length[counter];      
          double what_spring_constant;
          if (SCALE_TYPE == 0){
              what_spring_constant = spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*gausssigma))*exp(-(scaling_per_edge[counter]*scaling_per_edge[counter])/gausssigma)));
              if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant;}
          }
          else if (SCALE_TYPE == 1){
              what_spring_constant = spring_constant_weak*pow(scaling_per_edge[counter],scaling_pow) +
               spring_constant_weak*(1-pow(scaling_per_edge[counter], scaling_pow));
          }
          else if (SCALE_TYPE == 2){
              what_spring_constant = spring_constant - (spring_constant - spring_constant_weak)*scaling_per_edge[counter];
          }
          else if (SCALE_TYPE == 3){
              if (edges_in_upperhem[counter] == 1){
                  what_spring_constant = spring_constant_weak;
                  //length_zero = length_zero_growth;
              }
              else if (edges_in_upperhem[counter] == 0){
                  what_spring_constant = (spring_constant_weak + spring_constant)/2.0;
              }
              else{
                  what_spring_constant = spring_constant;
              }
          }
          else if (SCALE_TYPE == 4){
              if (nonuniform_wall_weakening_linear == true){
                  //double scaling = 0.0;//spring_constant_weak/spring_constant;
                  double spectrum = maxSpringScaler_linear*spring_constant - spring_constant_weak;
                  //what_spring_constant = spring_constant*((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[counter], hilleqnpow)))*(1-scaling) + scaling);
                  what_spring_constant = spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[counter], hilleqnpow)))*spectrum);
                  if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant_weak;}
              }
              else{
                  if (edges_in_upperhem[counter] == 1){
                      what_spring_constant = spring_constant_weak; // currently we're making the top and bottom layers weak while keeping the spring constant for the vertical springs very high. 
                  //length_zero = length_zero_growth;
                  }
                  else if (edges_in_upperhem[counter] == -1){
                      what_spring_constant = spring_constant_vertical;//(spring_constant_weak + spring_constant)/2.0;
                  }
                  else{
                      what_spring_constant = spring_constant;
                  }
              }
          }  
  
      
          //double length_zero = thrust::get<3>(u3d);
  
          // compute forces.
          double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
          double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
          double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
  
  
          double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                                      (yLoc_LR) * (yLoc_LR)  + 
                                      (zLoc_LR) * (zLoc_LR) );
  
          double energy = 0.0;
          if (length_current != length_zero){
              //double magnitude = -(what_spring_constant/(length_zero*length_zero)) * (length_current - length_zero);
              double magnitude = -(what_spring_constant) * (length_current - length_zero);
      
              idKey[place] = edgeL;
              printf("what_spring_constant = %f\n", what_spring_constant);
              printf("magnitude = %f\n", magnitude);
              printf("Linear Springs length_current = %f\n", length_current);
              printf("xLoc_LR = %f\n", xLoc_LR);
              printf("yLoc_LR = %f\n", yLoc_LR);
              printf("zLoc_LR = %f\n", zLoc_LR);
              
              //issue here writing to vectors with force????
              forceXAddr[place] = magnitude * (xLoc_LR/length_current);
              forceYAddr[place] = magnitude * (yLoc_LR/length_current);
              forceZAddr[place] = magnitude * (zLoc_LR/length_current);
  
              idKey[place + 1] = edgeR;
              forceXAddr[place + 1] = -magnitude * (xLoc_LR/length_current);
              forceYAddr[place + 1] = -magnitude * (yLoc_LR/length_current);
              forceZAddr[place + 1] = -magnitude * (zLoc_LR/length_current);
              
             // std::cout<< "Force at node "<< place << " = " << forceXAddr[place + 1] << ", "<<forceYAddr[place + 1]<< ", "<< forceZAddr[place + 1] <<std::endl;
      
              //energy = (what_spring_constant/(2.0*length_zero*length_zero)) * (length_current - length_zero) * (length_current - length_zero);
               energy = (what_spring_constant/(2.0)) * (length_current - length_zero) * (length_current - length_zero);
          }
          return energy;
      }
  
      else{
          double energy = 0.0;
          return energy;
      }
  }

};

#endif

