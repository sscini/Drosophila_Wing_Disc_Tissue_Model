#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_ 

#include "SystemStructures.h"
#include <math.h>

// Function declaration for computing volume springs.
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
	  LinearSpringInfoVecs& linearSpringInfoVecs, 
	  CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);


struct VolumeSpringFunctor {
      // Member variables representing parameters required by the functor.
      double current_total_volume;
      double true_current_total_volume;
      double eq_total_volume;
      double spring_constant;
      int num_of_triangles;
      double Rmin;

      int* triangles2Nodes_1;
      int* triangles2Nodes_2;
      int* triangles2Nodes_3;
      double* membraneNodeXAddr;
      double* membraneNodeYAddr;
      double* membraneNodeZAddr;
  
      int* id_value_expanded;
      int* keyBegin;
      int* keyEnd;
      
      // Constructor for the VolumeSpringsFunctor, which initiates member variables with the provided arguments. 
      __host__ __device__
      VolumeSpringFunctor(
          double& _current_total_volume,
          double& _true_current_total_volume,
          double& _eq_total_volume,
          double& _spring_constant,
          int& _num_of_triangles,
          double& _Rmin,
  
          int* _triangles2Nodes_1,
          int* _triangles2Nodes_2,
          int* _triangles2Nodes_3,
          double* _membraneNodeXAddr,
          double* _membraneNodeYAddr,
          double* _membraneNodeZAddr,
  
          int* _id_value_expanded,
          int* _keyBegin,
          int* _keyEnd):
          
          //initialize member variables with constructor arguments.
          
          current_total_volume(_current_total_volume),
          true_current_total_volume(_true_current_total_volume),
          eq_total_volume(_eq_total_volume),
          spring_constant(_spring_constant),
          num_of_triangles(_num_of_triangles),
          Rmin(_Rmin),
  
          triangles2Nodes_1(_triangles2Nodes_1),
          triangles2Nodes_2(_triangles2Nodes_2),
          triangles2Nodes_3(_triangles2Nodes_3),
          membraneNodeXAddr(_membraneNodeXAddr),
          membraneNodeYAddr(_membraneNodeYAddr),
          membraneNodeZAddr(_membraneNodeZAddr),
  
          id_value_expanded(_id_value_expanded),
          keyBegin(_keyBegin),
          keyEnd(_keyEnd) {}
          
          // The functor's operator() is called for each node and calculates the volume spring forces.
    __device__ 
    CVec3 operator()(const U2CVec3& u2d3) {
    
        int node_id = thrust::get<0>(u2d3);
        
        int bucketId = thrust::get<1>(u2d3); // bucket containing nodeId

        double r1_dx = 0.0;
        double r1_dy = 0.0;
        double r1_dz = 0.0;
        double dV_r1_dx = 0.0;
        double dV_r1_dy = 0.0;
        double dV_r1_dz = 0.0;
        
        // Initialize variables to accumulate force components for the current node.
        double forceX = thrust::get<2>(u2d3);
        double forceY = thrust::get<3>(u2d3);
        double forceZ = thrust::get<4>(u2d3);

        // Loop through all triangles, calculating forces only for the triangles that contain the current node.
        for (int i = 0; i < num_of_triangles; i++){ 
            // Determine the triangle's vertices (r1, r2, and r3) based on whether the current node belongs to the triangle.
            int r1, r2, r3;
        
            // Check if the current node is the first vertex of the triangle.
            if (triangles2Nodes_1[i] == node_id){
                r1 = triangles2Nodes_1[i];
                r2 = triangles2Nodes_2[i];
                r3 = triangles2Nodes_3[i];
                }
            // Check if the current node is the second vertex of the triangle.
            else if (triangles2Nodes_2[i] == node_id){
                r1 = triangles2Nodes_2[i];
                r2 = triangles2Nodes_3[i];
                r3 = triangles2Nodes_1[i];
                }
            // Check if the current node is the third vertex of the triangle.
            else if (triangles2Nodes_3[i] == node_id){
                r1 = triangles2Nodes_3[i];
                r2 = triangles2Nodes_1[i];
                r3 = triangles2Nodes_2[i];
                }   
            else{ // If the current node does not belong to the current triangle, skip the computation.
                forceX += 0.0;
                forceY += 0.0;
                forceZ += 0.0;
                continue;
                }

            // Check if all three vertices (r1, r2, and r3) of the triangle are valid (not INT_MAX).
            if (r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
  
                // Extract coordinates of the vertices of the triangle.
                double r1x = membraneNodeXAddr[r1];
                double r1y = membraneNodeYAddr[r1];
                double r1z = membraneNodeZAddr[r1];
                double r2x = membraneNodeXAddr[r2];
                double r2y = membraneNodeYAddr[r2];
                double r2z = membraneNodeZAddr[r2];
                double r3x = membraneNodeXAddr[r3];
                double r3y = membraneNodeYAddr[r3];
                double r3z = membraneNodeZAddr[r3];
            // Compute the normal vector components of the triangle (un-normalized).
                double n1 = (r2y - r1y)*(r3z - r1z) - (r3y - r1y)*(r2z - r1z);
                double n2 = -(r2x - r1x)*(r3z - r1z) + (r3x - r1x)*(r2z - r1z);
                double n3 = (r2x - r1x)*(r3y - r1y) - (r3x - r1x)*(r2y - r1y);
                double norm_n = sqrt(n1*n1 + n2*n2 + n3*n3);
      
                // Normalize the components of the normal vector.
                double N1 = n1/norm_n;
                double N2 = n2/norm_n;
                double N3 = n3/norm_n;

                // Compute derivatives of the components of the normal vector with respect to vertex r1 (r1x, r1y, r1z).
                double n1_x = 0.0;
                double n1_y = -(r3z - r1z) + (r2z - r1z);
                double n1_z = -(r2y - r1y) + (r3y - r1y);        
                double n2_x = (r3z - r1z) - (r2z - r1z);
                double n2_y = 0.0;
                double n2_z = (r2x - r1x) - (r3x - r1x);       
                double n3_x = -(r3y - r1y) + (r2y - r1y);
                double n3_y = -(r2x - r1x) + (r3x - r1x);
                double n3_z = 0.0;

                // Compute derivatives of normalized normal vector components (N1, N2, N3) with respect to vertex r1.
                double N1_x = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_x - n1*(1.0/norm_n)*(n1*n1_x + n2*n2_x + n3*n3_x));
                double N1_y = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_y - n1*(1.0/norm_n)*(n1*n1_y + n2*n2_y + n3*n3_y));
                double N1_z = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_z - n1*(1.0/norm_n)*(n1*n1_z + n2*n2_z));
                // ... The rest of the function was not provided, but you should continue it here.
                
                // You need to return a CVec3 at the end of the function.
                return CVec3(/*...*/);
            }
        }
        return CVec3(/* default return value if loop doesn't execute*/);
    }
};
//      
//};
//// Definition of the functor used to calculate volume spring forces.
//struct VolumeSpringFunctor : public thrust::unary_function<U2CVec3,CVec3> {
//      
//      // Member variables representing parameters required by the functor.
//      double current_total_volume;
//      double true_current_total_volume;
//      double eq_total_volume;
//      double spring_constant;
//      int num_of_triangles;
//      double Rmin;
//
//      int* triangles2Nodes_1;
//      int* triangles2Nodes_2;
//      int* triangles2Nodes_3;
//      double* membraneNodeXAddr;
//      double* membraneNodeYAddr;
//      double* membraneNodeZAddr;
//  
//      int* id_value_expanded;
//      int* keyBegin;
//      int* keyEnd;
//    
//	    // Constructor for the VolumeSpringFunctor, which initializes member variables with the provided arguments.
//      __host__ __device__ 
//      VolumeSpringFunctor(
//          double& _current_total_volume,
//          double& _true_current_total_volume,
//          double& _eq_total_volume,
//          double& _spring_constant,
//          int& _num_of_triangles,
//          double& _Rmin,
//
//          int* _triangles2Nodes_1,
//          int* _triangles2Nodes_2,
//          int* _triangles2Nodes_3,
//          double* _membraneNodeXAddr,
//          double* _membraneNodeYAddr,
//          double* _membraneNodeZAddr,
//
//          int* _id_value_expanded,
//          int* _keyBegin,
//          int* _keyEnd):
//
//          // Initialize member variables with constructor arguments.
//          current_total_volume(_current_total_volume),
//          true_current_total_volume(_true_current_total_volume),
//          eq_total_volume(_eq_total_volume),
//          spring_constant(_spring_constant),
//          num_of_triangles(_num_of_triangles),
//          Rmin(_Rmin),
//
//          triangles2Nodes_1(_triangles2Nodes_1),
//          triangles2Nodes_2(_triangles2Nodes_2),
//          triangles2Nodes_3(_triangles2Nodes_3),
//          membraneNodeXAddr(_membraneNodeXAddr),
//          membraneNodeYAddr(_membraneNodeYAddr),
//          membraneNodeZAddr(_membraneNodeZAddr),
//
//          id_value_expanded(_id_value_expanded),
//          keyBegin(_keyBegin),
//          keyEnd(_keyEnd) {}
//
//      //hand in counting iterator and id of two edges and preferred length
//      // Functor operator for the VolumeSpringFunctor, which calculates the volume spring forces for each node.
//      __device__ 
//      CVec3 operator()(const U2CVec3& u2d3) {
//      // The functor's operator() is called for each node and calculates the volume spring forces.
//
//          int node_id = thrust::get<0>(u2d3);
//          
//          int bucketId = thrust::get<1>(u2d3);//bucket containing nodeId
//
//          double r1_dx = 0.0;
//          double r1_dy = 0.0;
//          double r1_dz = 0.0;
//          double dV_r1_dx = 0.0;
//          double dV_r1_dy = 0.0;
//          double dV_r1_dz = 0.0;
//          //double forceX, forceY, forceZ;
//          
//          // Initialize variables to accumulate force components for the current node.
//          double forceX = thrust::get<2>(u2d3);
//          double forceY = thrust::get<3>(u2d3);
//          double forceZ = thrust::get<4>(u2d3);
//
//          // Loop through all triangles, calculating forces only for the triangles that contain the current node.
//          for (int i = 0; i < num_of_triangles; i++){ 
//              // Determine the triangle's vertices (r1, r2, and r3) based on whether the current node belongs to the triangle.
//              // Here, we always set the "counter" vertex to be node 1 (or r1 in this case) for derivative calculations.
//              //Loop through all triangles, but only do computation if the triangle contains the "counter" vertex.
//              //Here we always set the "counter" vertex to be node 1 (or r1 in this case).
//              //We will take the derivative with respect to r1.
//              int r1, r2, r3;
//          
//              // Check if the current node is the first vertex of the triangle.
//              if (triangles2Nodes_1[i] == node_id){
//                  r1 = triangles2Nodes_1[i];
//                  r2 = triangles2Nodes_2[i];
//                  r3 = triangles2Nodes_3[i];
//                  }
//              // Check if the current node is the second vertex of the triangle.
//              else if (triangles2Nodes_2[i] == node_id){
//                  r1 = triangles2Nodes_2[i];
//                  r2 = triangles2Nodes_3[i];
//                  r3 = triangles2Nodes_1[i];
//                  }
//              // Check if the current node is the third vertex of the triangle.
//              else if (triangles2Nodes_3[i] == node_id){
//                  r1 = triangles2Nodes_3[i];
//                  r2 = triangles2Nodes_1[i];
//                  r3 = triangles2Nodes_2[i];
//                  }   
//              else{// If the current node does not belong to the current triangle, skip the computation for this triangle.
//                
//                  forceX += 0.0;
//                  forceY += 0.0;
//                  forceZ += 0.0;
//                  continue;
//                  
//                  }//if "counter" vertex does not belong to the current triangle, skip the computation.
//
//              // Check if all three vertices (r1, r2, and r3) of the triangle are valid (not INT_MAX).
//              if (r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
//    
//                  // Extract coordinates of the vertices of the triangle.
//                  double r1x = membraneNodeXAddr[r1];
//                  double r1y = membraneNodeYAddr[r1];
//                  double r1z = membraneNodeZAddr[r1];
//                  double r2x = membraneNodeXAddr[r2];
//                  double r2y = membraneNodeYAddr[r2];
//                  double r2z = membraneNodeZAddr[r2];
//                  double r3x = membraneNodeXAddr[r3];
//                  double r3y = membraneNodeYAddr[r3];
//                  double r3z = membraneNodeZAddr[r3];
//
//                  // Compute the normal vector components of the triangle (un-normalized).
//                  double n1 = (r2y - r1y)*(r3z - r1z) - (r3y - r1y)*(r2z - r1z);//un-normalized components of the normal vector
//                  double n2 = -(r2x - r1x)*(r3z - r1z) + (r3x - r1x)*(r2z - r1z);
//                  double n3 = (r2x - r1x)*(r3y - r1y) - (r3x - r1x)*(r2y - r1y);
//                  double norm_n = sqrt(n1*n1 + n2*n2 + n3*n3);
//        
//                  // Normalize the components of the normal vector.
//                  double N1 = n1/norm_n;
//                  double N2 = n2/norm_n;
//                  double N3 = n3/norm_n;
//
//                  // Compute derivatives of the components of the normal vector with respect to vertex r1 (r1x, r1y, r1z).
//                  double n1_x = 0.0;
//                  double n1_y = -(r3z - r1z) + (r2z - r1z);
//                  double n1_z = -(r2y - r1y) + (r3y - r1y);        
//                  double n2_x = (r3z - r1z) - (r2z - r1z);
//                  double n2_y = 0.0;
//                  double n2_z = (r2x - r1x) - (r3x - r1x);       
//                  double n3_x = -(r3y - r1y) + (r2y - r1y);
//                  double n3_y = -(r2x - r1x) + (r3x - r1x);
//                  double n3_z = 0.0;
//
//                  // Compute derivatives of normalized normal vector components (N1, N2, N3) with respect to vertex r1.
//                  double N1_x = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_x - n1*(1.0/norm_n)*(n1*n1_x + n2*n2_x + n3*n3_x));
//                  double N1_y = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_y - n1*(1.0/norm_n)*(n1*n1_y + n2*n2_y + n3*n3_y));
//                  double N1_z = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n1_z - n1*(1.0/norm_n)*(n1*n1_z + n2*n2_z + n3*n3_z));
//                  double N2_x = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n2_x - n2*(1.0/norm_n)*(n1*n1_x + n2*n2_x + n3*n3_x));
//                  double N2_y = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n2_y - n2*(1.0/norm_n)*(n1*n1_y + n2*n2_y + n3*n3_y));
//                  double N2_z = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n2_z - n2*(1.0/norm_n)*(n1*n1_z + n2*n2_z + n3*n3_z));
//                  double N3_x = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n3_x - n3*(1.0/norm_n)*(n1*n1_x + n2*n2_x + n3*n3_x));
//                  double N3_y = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n3_y - n3*(1.0/norm_n)*(n1*n1_y + n2*n2_y + n3*n3_y));
//                  double N3_z = (1.0/norm_n)*(1.0/norm_n)*(norm_n*n3_z - n3*(1.0/norm_n)*(n1*n1_z + n2*n2_z + n3*n3_z));
//                  
//
//                  double r1dN = r1x*N1 + r1y*N2 + r1z*N3;
//        
//                  // Cross products of vectors r1 and r2 for the calculation of derivatives.
//                  double r1cr2x = r1y*r2z - r2y*r1z;
//                  double r1cr2x_x = 0.0;
//                  double r1cr2x_y = r2z;
//                  double r1cr2x_z = -r2y;
//                  
//                  double r1cr2y = -r1x*r2z + r2x*r1z;
//                  double r1cr2y_x = -r2z;
//                  double r1cr2y_y = 0.0;
//                  double r1cr2y_z = r2x;
//                  
//                  double r1cr2z = r1x*r2y - r2x*r1y;
//                  double r1cr2z_x = r2y;
//                  double r1cr2z_y = -r2x;
//                  double r1cr2z_z = 0.0;
//                  
//                  // Cross products of vectors r2 and r3 for the calculation of derivatives.
//                  double r2cr3x = r2y*r3z - r3y*r2z;
//                  double r2cr3x_x = 0.0;
//                  double r2cr3x_y = 0.0;
//                  double r2cr3x_z = 0.0;
//        
//                  double r2cr3y = -r2x*r3z + r3x*r2z;
//                  double r2cr3y_x = 0.0;
//                  double r2cr3y_y = 0.0;
//                  double r2cr3y_z = 0.0;
//                  
//                  double r2cr3z = r2x*r3y - r3x*r2y;
//                  double r2cr3z_x = 0.0;
//                  double r2cr3z_y = 0.0;
//                  double r2cr3z_z = 0.0;
//        
//                  // Cross products of vectors r3 and r1 for the calculation of derivatives.
//                  double r3cr1x = r3y*r1z - r1y*r3z;
//                  double r3cr1x_x = 0.0;
//                  double r3cr1x_y = -r3z;
//                  double r3cr1x_z = r3y;
//                  
//                  double r3cr1y = -r3x*r1z + r1x*r3z;
//                  double r3cr1y_x = r3z;
//                  double r3cr1y_y = 0.0;
//                  double r3cr1y_z = -r3x;
//                  
//                  double r3cr1z = r3x*r1y - r1x*r3y;
//                  double r3cr1z_x = -r3y;
//                  double r3cr1z_y = r3x;
//                  double r3cr1z_z = 0.0;
//          
//                  // Compute the dot product of N1, N2, and N3 with (r1cr2 + r2cr3 + r3cr1) for the calculation of derivatives.
//                  double NN = N1*(r1cr2x + r2cr3x + r3cr1x) + N2*(r1cr2y + r2cr3y + r3cr1y) + N3*(r1cr2z + r2cr3z + r3cr1z);
//
//                  // Compute derivatives of r1_dot_N (r1dN) with respect to r1x, r1y, and r1z.
//                  //Derivatives of r1_dot_N (r1dN) with respect to r1x, y, z
//                  double r1dN_x = r1x*N1_x + 1.0*N1 + r1y*N2_x + 0.0*N2 + r1z*N3_x + 0.0*N3; //N1 + r1y*((r3z - r1z) - (r2z - r1z)) + r1z*(-(r3y - r1y) + (r2y - r1y));
//                  double r1dN_y = r1x*N1_y + 0.0*N1 + r1y*N2_y + 1.0*N2 + r1z*N3_y + 0.0*N3; //r1x*(-(r3z - r1z) + (r2z - r1z)) + N2 + r1z*(-(r2x - r1x) + (r3x - r1x));
//                  double r1dN_z = r1x*N1_z + 0.0*N1 + r1y*N2_z + 0.0*N2 + r1z*N3_z + 1.0*N3; //r1x*(-(r2y - r1y) + (r3y - r1y)) + r1y*((r2x - r1x) - (r3x - r1x)) + N3;
//
//                  //Derivatives of NN
//                  double NN_x = N1_x*(r1cr2x + r2cr3x + r3cr1x) + N1*(r1cr2x_x + r2cr3x_x + r3cr1x_x) +
//                        N2_x*(r1cr2y + r2cr3y + r3cr1y) + N2*(r1cr2y_x + r2cr3y_x + r3cr1y_x) +
//                        N3_x*(r1cr2z + r2cr3z + r3cr1z) + N3*(r1cr2z_x + r2cr3z_x + r3cr1z_x);
//                  double NN_y = N1_y*(r1cr2x + r2cr3x + r3cr1x) + N1*(r1cr2x_y + r2cr3x_y + r3cr1x_y) +
//                        N2_y*(r1cr2y + r2cr3y + r3cr1y) + N2*(r1cr2y_y + r2cr3y_y + r3cr1y_y) +
//                        N3_y*(r1cr2z + r2cr3z + r3cr1z) + N3*(r1cr2z_y + r2cr3z_y + r3cr1z_y);
//                  double NN_z = N1_z*(r1cr2x + r2cr3x + r3cr1x) + N1*(r1cr2x_z + r2cr3x_z + r3cr1x_z) +
//                        N2_z*(r1cr2y + r2cr3y + r3cr1y) + N2*(r1cr2y_z + r2cr3y_z + r3cr1y_z) +
//                        N3_z*(r1cr2z + r2cr3z + r3cr1z) + N3*(r1cr2z_z + r2cr3z_z + r3cr1z_z);
//
//                  // Accumulate the derivatives for calculating the total derivative with respect to r1x, r1y, and r1z.
//                  r1_dx += r1dN_x*sqrt(NN*NN) + r1dN*(NN/sqrt(NN*NN))*NN_x;
//                  r1_dy += r1dN_y*sqrt(NN*NN) + r1dN*(NN/sqrt(NN*NN))*NN_y;
//                  r1_dz += r1dN_z*sqrt(NN*NN) + r1dN*(NN/sqrt(NN*NN))*NN_z;
//
//                  // Calculate the volume derivative components for r1x, r1y, and r1z.
//                  dV_r1_dx = (1.0/12.0)*(2.0*current_total_volume/sqrt(current_total_volume*current_total_volume))*r1_dx;
//                  dV_r1_dy = (1.0/12.0)*(2.0*current_total_volume/sqrt(current_total_volume*current_total_volume))*r1_dy;
//                  dV_r1_dz = (1.0/12.0)*(2.0*current_total_volume/sqrt(current_total_volume*current_total_volume))*r1_dz; 
//                  }
//              else{
//                  // Continue to the next iteration if the triangle is not valid (contains INT_MAX vertices).
//                  continue;
//                  }
//
//              }
//
//          // Calculate the magnitude of the spring force.
//          double magnitude = (spring_constant/(2.0*Rmin*Rmin*Rmin*eq_total_volume))*2.0*(true_current_total_volume - eq_total_volume);
//
//          // Update the force components with the volume spring force contribution.
//          forceX += magnitude*(-dV_r1_dx);
//          forceY += magnitude*(-dV_r1_dy);
//          forceZ += magnitude*(-dV_r1_dz);
//
//
//
//
//          // Return the updated force components as a tuple.
//          return thrust::make_tuple(forceX, forceY, forceZ);
//
//
//          }
//      };

#endif




#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"
#include <math.h>

// Function declaration for computing volume springs.
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

// Definition of the functor used to calculate volume spring forces.
struct VolumeSpringFunctor : public thrust::unary_function<U2CVec3, CVec3> {
    // Member variables representing parameters required by the functor.
    double current_total_volume;
    double true_current_total_volume;
    double eq_total_volume;
    double spring_constant;
    int num_of_triangles;
    double Rmin;

    int* triangles2Nodes_1;
    int* triangles2Nodes_2;
    int* triangles2Nodes_3;
    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    int* id_value_expanded;
    int* keyBegin;
    int* keyEnd;

    // Constructor for the VolumeSpringFunctor, which initializes member variables with the provided arguments.
    __host__ __device__
    VolumeSpringFunctor(
        double& _current_total_volume,
        double& _true_current_total_volume,
        double& _eq_total_volume,
        double& _spring_constant,
        int& _num_of_triangles,
        double& _Rmin,

        int* _triangles2Nodes_1,
        int* _triangles2Nodes_2,
        int* _triangles2Nodes_3,
        double* _membraneNodeXAddr,
        double* _membraneNodeYAddr,
        double* _membraneNodeZAddr,

        int* _id_value_expanded,
        int* _keyBegin,
        int* _keyEnd) :

        // Initialize member variables with constructor arguments.
        current_total_volume(_current_total_volume),
        true_current_total_volume(_true_current_total_volume),
        eq_total_volume(_eq_total_volume),
        spring_constant(_spring_constant),
        num_of_triangles(_num_of_triangles),
        Rmin(_Rmin),
        triangles2Nodes_1(_triangles2Nodes_1),
        triangles2Nodes_2(_triangles2Nodes_2),
        triangles2Nodes_3(_triangles2Nodes_3),
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),
        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

    // Functor operator for the VolumeSpringFunctor, which calculates the volume spring forces for each node.
    __device__
    CVec3 operator()(const U2CVec3& u2d3) {
        // The functor's operator() is called for each node and calculates the volume spring forces.

        int node_id = thrust::get<0>(u2d3);
        int bucketId = thrust::get<1>(u2d3); // Bucket containing the nodeId

        // Initialize variables to accumulate force components for the current node.
        double forceX = thrust::get<2>(u2d3);
        double forceY = thrust::get<3>(u2d3);
        double forceZ = thrust::get<4>(u2d3);

        // Loop through all triangles, calculating forces only for the triangles that contain the current node.
        for (int i = 0; i < num_of_triangles; i++) {
            // Determine the triangle's vertices (r1, r2, and r3) based on whether the current node belongs to the triangle.
            // Here, we always set the "counter" vertex to be node 1 (or r1 in this case) for derivative calculations.
            int r1, r2, r3;

            // Check if the current node is the first vertex of the triangle.
            if (triangles2Nodes_1[i] == node_id) {
                r1 = triangles2Nodes_1[i];
                r2 = triangles2Nodes_2[i];
                r3 = triangles2Nodes_3[i];
            }
            // Check if the current node is the second vertex of the triangle.
            else if (triangles2Nodes_2[i] == node_id) {
                r1 = triangles2Nodes_2[i];
                r2 = triangles2Nodes_3[i];
                r3 = triangles2Nodes_1[i];
            }
            // Check if the current node is the third vertex of the triangle.
            else if (triangles2Nodes_3[i] == node_id) {
                r1 = triangles2Nodes_3[i];
                r2 = triangles2Nodes_1[i];
                r3 = triangles2Nodes_2[i];
            } else {
                // If the current node does not belong to the current triangle, skip the computation for this triangle.
                forceX += 0.0;
                forceY += 0.0;
                forceZ += 0.0;
                continue;
            }

            // Check if all three vertices (r1, r2, and r3) of the triangle are valid (not INT_MAX).
            if (r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX) {
                // Extract coordinates of the vertices of the triangle.
                double r1x = membraneNodeXAddr[r1];
                double r1y = membraneNodeYAddr[r1];
                double r1z = membraneNodeZAddr[r1];
                double r2x = membraneNodeXAddr[r2];
                double r2y = membraneNodeYAddr[r2];
                double r2z = membraneNodeZAddr[r2];
                double r3x = membraneNodeXAddr[r3];
                double r3y = membraneNodeYAddr[r3];
                double r3z = membraneNodeZAddr[r3];

                // Compute the normal vector components of the triangle (un-normalized).
                double n1 = (r2y - r1y) * (r3z - r1z) - (r3y - r1y) * (r2z - r1z);
                double n2 = -(r2x - r1x) * (r3z - r1z) + (r3x - r1x) * (r2z - r1z);
                double n3 = (r2x - r1x) * (r3y - r1y) - (r3x - r1x) * (r2y - r1y);
                double norm_n = sqrt(n1 * n1 + n2 * n2 + n3 * n3);

                // Normalize the components of the normal vector.
                double N1 = n1 / norm_n;
                double N2 = n2 / norm_n;
                double N3 = n3 / norm_n;

                // Compute derivatives of the components of the normal vector with respect to vertex r1 (r1x, r1y, r1z).
                double n1_x = 0.0;
                double n1_y = -(r3z - r1z) + (r2z - r1z);
                double n1_z = -(r2y - r1y) + (r3y - r1y);
                double n2_x = (r3z - r1z) - (r2z - r1z);
                double n2_y = 0.0;
                double n2_z = (r2x - r1x) - (r3x - r1x);
                double n3_x = -(r3y - r1y) + (r2y - r1y);
                double n3_y = -(r2x - r1x) + (r3x - r1x);
                double n3_z = 0.0;

                // Compute derivatives of normalized normal vector components (N1, N2, N3) with respect to vertex r1.
                double N1_x = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n1_x - n1 * (1.0 / norm_n) * (n1 * n1_x + n2 * n2_x + n3 * n3_x));
                double N1_y = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n1_y - n1 * (1.0 / norm_n) * (n1 * n1_y + n2 * n2_y + n3 * n3_y));
                double N1_z = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n1_z - n1 * (1.0 / norm_n) * (n1 * n1_z + n2 * n2_z + n3 * n3_z));
                double N2_x = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n2_x - n2 * (1.0 / norm_n) * (n1 * n1_x + n2 * n2_x + n3 * n3_x));
                double N2_y = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n2_y - n2 * (1.0 / norm_n) * (n1 * n1_y + n2 * n2_y + n3 * n3_y));
                double N2_z = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n2_z - n2 * (1.0 / norm_n) * (n1 * n1_z + n2 * n2_z + n3 * n3_z));
                double N3_x = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n3_x - n3 * (1.0 / norm_n) * (n1 * n1_x + n2 * n2_x + n3 * n3_x));
                double N3_y = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n3_y - n3 * (1.0 / norm_n) * (n1 * n1_y + n2 * n2_y + n3 * n3_y));
                double N3_z = (1.0 / norm_n) * (1.0 / norm_n) * (norm_n * n3_z - n3 * (1.0 / norm_n) * (n1 * n1_z + n2 * n2_z + n3 * n3_z));

                double r1dN = r1x * N1 + r1y * N2 + r1z * N3;

                // Cross products of vectors r1 and r2 for the calculation of derivatives.
                double r1cr2x = r1y * r2z - r2y * r1z;
                double r1cr2x_x = 0.0;
                double r1cr2x_y = r2z;
                double r1cr2x_z = -r2y;

                double r1cr2y = -r1x * r2z + r2x * r1z;
                double r1cr2y_x = -r2z;
                double r1cr2y_y = 0.0;
                double r1cr2y_z = r2x;

                double r1cr2z = r1x * r2y - r2x * r1y;
                double r1cr2z_x = r2y;
                double r1cr2z_y = -r2x;
                double r1cr2z_z = 0.0;

                // Cross products of vectors r2 and r3 for the calculation of derivatives.
                double r2cr3x = r2y * r3z - r3y * r2z;
                double r2cr3x_x = 0.0;
                double r2cr3x_y = 0.0;
                double r2cr3x_z = 0.0;

                double r2cr3y = -r2x * r3z + r3x * r2z;
                double r2cr3y_x = 0.0;
                double r2cr3y_y = 0.0;
                double r2cr3y_z = 0.0;

                double r2cr3z = r2x * r3y - r3x * r2y;
                double r2cr3z_x = 0.0;
                double r2cr3z_y = 0.0;
                double r2cr3z_z = 0.0;

                // Cross products of vectors r3 and r1 for the calculation of derivatives.
                double r3cr1x = r3y * r1z - r1y * r3z;
                double r3cr1x_x = 0.0;
                double r3cr1x_y = -r3z;
                double r3cr1x_z = r3y;

                double r3cr1y = -r3x * r1z + r1x * r3z;
                double r3cr1y_x = r3z;
                double r3cr1y_y = 0.0;
                double r3cr1y_z = -r3x;

                double r3cr1z = r3x * r1y - r1x * r3y;
                double r3cr1z_x = -r3y;
                double r3cr1z_y = r3x;
                double r3cr1z_z = 0.0;

                // Compute the dot product of N1, N2, and N3 with (r1cr2 + r2cr3 + r3cr1) for the calculation of derivatives.
                double NN = N1 * (r1cr2x + r2cr3x + r3cr1x) + N2 * (r1cr2y + r2cr3y + r3cr1y) + N3 * (r1cr2z + r2cr3z + r3cr1z);

                // Compute derivatives of r1_dot_N (r1dN) with respect to r1x, r1y, and r1z.
                double r1dN_x = r1x * N1_x + 1.0 * N1 + r1y * N2_x + 0.0 * N2 + r1z * N3_x + 0.0 * N3;
                double r1dN_y = r1x * N1_y + 0.0 * N1 + r1y * N2_y + 1.0 * N2 + r1z * N3_y + 0.0 * N3;
                double r1dN_z = r1x * N1_z + 0.0 * N1 + r1y * N2_z + 0.0 * N2 + r1z * N3_z + 1.0 * N3;

                // Accumulate the derivatives for calculating the total derivative with respect to r1x, r1y, and r1z.
                r1_dx += r1dN_x * sqrt(NN * NN) + r1dN * (NN / sqrt(NN * NN)) * NN_x;
                r1_dy += r1dN_y * sqrt(NN * NN) + r1dN * (NN / sqrt(NN * NN)) * NN_y;
                r1_dz += r1dN_z * sqrt(NN * NN) + r1dN * (NN / sqrt(NN * NN)) * NN_z;

                // Calculate the volume derivative components for r1x, r1y, and r1z.
                dV_r1_dx = (1.0 / 12.0) * (2.0 * current_total_volume / sqrt(current_total_volume * current_total_volume)) * r1_dx;
                dV_r1_dy = (1.0 / 12.0) * (2.0 * current_total_volume / sqrt(current_total_volume * current_total_volume)) * r1_dy;
                dV_r1_dz = (1.0 / 12.0) * (2.0 * current_total_volume / sqrt(current_total_volume * current_total_volume)) * r1_dz;
            } else {
                // Continue to the next iteration if the triangle is not valid (contains INT_MAX vertices).
                continue;
            }
        }

        // Calculate the magnitude of the spring force.
        double magnitude = (spring_constant / (2.0 * Rmin * Rmin * Rmin * eq_total_volume)) * 2.0 * (true_current_total_volume - eq_total_volume);

        // Update the force components with the volume spring force contribution.
        forceX += magnitude * (-dV_r1_dx);
        forceY += magnitude * (-dV_r1_dy);
        forceZ += magnitude * (-dV_r1_dz);

        // Return the updated force components as a tuple.
        return thrust::make_tuple(forceX, forceY, forceZ);
    }
};

#endif

