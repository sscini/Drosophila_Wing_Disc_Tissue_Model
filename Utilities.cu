#include "System.h"
#include <random>
#include "Utilities.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <omp.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//REMEMBER TO CHANGE THE NEXT LINE IF YOU CHANGE THE ACCEPTANCE RULE!
//CURRENT SETUP: swap is always accepted for boltzmann.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Utilities::Utilities(CoordInfoVecs& coordInfoVecs,
GeneralParams& generalParams) {
            
};


int Utilities::surfaceNormal_device_vecs(
    int inode,
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams
){
    double current_normalX = 0.0;
    double current_normalY = 0.0;
    double current_normalZ = 0.0;
    int neighbor;
    double AREA;
    double scaled_forceX, scaled_forceY, scaled_forceZ;
    for (int j = 0; j < 9; j ++){
        if (j == 0 && (coordInfoVecs.nodes2Triangles_1[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_1[inode];}
        //else{continue;}
        else if (j == 1 && (coordInfoVecs.nodes2Triangles_2[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_2[inode];}
        //else{continue;}
        else if (j == 2 && (coordInfoVecs.nodes2Triangles_3[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_3[inode];}
        //else{continue;}
        else if (j == 3 && (coordInfoVecs.nodes2Triangles_4[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_4[inode];}
        //else{continue;}
        else if (j == 4 && (coordInfoVecs.nodes2Triangles_5[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_5[inode];}
        //else{continue;}
        else if (j == 5 && (coordInfoVecs.nodes2Triangles_6[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_6[inode];}
        //else{continue;}
        else if (j == 6 && (coordInfoVecs.nodes2Triangles_7[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_7[inode];}
        //else{continue;}
        else if (j == 7 && (coordInfoVecs.nodes2Triangles_8[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_8[inode];}
        //else{continue;}
        else if (j == 8 && (coordInfoVecs.nodes2Triangles_9[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_9[inode];}
        //else{continue;}
       // else if (j == 9 && (coordInfoVecs.nodes2Triangles_10[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_10[inode];}
        //else{continue;}
       // else if (j == 10 && (coordInfoVecs.nodes2Triangles_11[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_11[inode];}
        //else{continue;}
        //else if (j == 11 && (coordInfoVecs.nodes2Triangles_12[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_12[inode];}
        else{continue;}
       //for (int i = 0; i < 12; i++){
        
        int node1 = coordInfoVecs. triangles2Nodes_1[neighbor];
        int node2 = coordInfoVecs. triangles2Nodes_2[neighbor];
        int node3 = coordInfoVecs. triangles2Nodes_3[neighbor];

        double vec1x = coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1];
        double vec1y = coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1];
        double vec1z = coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1];
        double vec2x = coordInfoVecs.nodeLocX[node3] - coordInfoVecs.nodeLocX[node1];
        double vec2y = coordInfoVecs.nodeLocY[node3] - coordInfoVecs.nodeLocY[node1];
        double vec2z = coordInfoVecs.nodeLocZ[node3] - coordInfoVecs.nodeLocZ[node1];

        double normalX = vec1y*vec2z - vec2y*vec1z;
        double normalY = -(vec1x*vec2z - vec2x*vec1z);
        double normalZ = vec1x*vec2y - vec2x*vec1y;
        double normalize = sqrt((normalX*normalX)+(normalY*normalY)+(normalZ*normalZ));
        normalX = normalX/normalize;
        normalY = normalY/normalize;
        normalZ = normalZ/normalize;
        AREA = normalize/2.0;
        scaled_forceX = ((1.0/3.0)*AREA*normalX);
        scaled_forceY = ((1.0/3.0)*AREA*normalY);
        scaled_forceZ = ((1.0/3.0)*AREA*normalZ);
        //current_normalX += normalX;
        //current_normalY += normalY;
        //current_normalZ += normalZ;
        coordInfoVecs.nodeForceX[inode] += generalParams.volume_spring_constant*scaled_forceX;
        coordInfoVecs.nodeForceY[inode] += generalParams.volume_spring_constant*scaled_forceY;
        coordInfoVecs.nodeForceZ[inode] += generalParams.volume_spring_constant*scaled_forceZ;
    }

   
};


int Utilities::nodes2Triangles_host_vecs(
    int inode,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
    AuxVecs& auxVecs
){

    hostSetInfoVecs.nodes2Triangles_1[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_2[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_3[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_4[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_5[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_6[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_7[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_8[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_9[inode] = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_10 = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_11 = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_12 = -INT_MAX;

    //for now iterate through all membrane id's
    // int begin = keyBegin[bucketId];
    // int end = keyEnd[bucketId];

    //First, create a new vector detailing all the neighboring nodes of node i//
    
        int COUNT = 0;
        for (int k = 0; k < coordInfoVecs.num_triangles; k++){
            int v1 = hostSetInfoVecs.triangles2Nodes_1[k];
            int v2 = hostSetInfoVecs.triangles2Nodes_2[k];
            int v3 = hostSetInfoVecs.triangles2Nodes_3[k];
            if (v1 < 0 || v2 < 0 || v3 < 0){
                continue;
            }
            else if (v1 == INT_MAX || v2 == INT_MAX || v3 == INT_MAX){
                continue;
            }

            if (v1 == inode || v2 == inode || v3 == inode){
                
                COUNT += 1;
                if (COUNT == 1){
                    hostSetInfoVecs.nodes2Triangles_1[inode] = k;
                }
                else if (COUNT == 2){
                    hostSetInfoVecs.nodes2Triangles_2[inode] = k;
                }
                else if (COUNT == 3){
                    hostSetInfoVecs.nodes2Triangles_3[inode] = k;
                }
                else if (COUNT == 4){
                    hostSetInfoVecs.nodes2Triangles_4[inode] = k;
                }
                else if (COUNT == 5){
                    hostSetInfoVecs.nodes2Triangles_5[inode] = k;
                }
                else if (COUNT == 6){
                    hostSetInfoVecs.nodes2Triangles_6[inode] = k;
                }
                else if (COUNT == 7){
                    hostSetInfoVecs.nodes2Triangles_7[inode] = k;
                }
                else if (COUNT == 8){
                    hostSetInfoVecs.nodes2Triangles_8[inode] = k;
                }
                else if (COUNT == 9){
                    hostSetInfoVecs.nodes2Triangles_9[inode] = k;
                }
                /*else if (COUNT == 10){
                    coordInfoVecs.nodes2Triangles_10 = k;
                }*/
                /* else if (COUNT == 11){
                    coordInfoVecs.nodes2Triangles_11 = k;
                }
                else if (COUNT == 12){
                    coordInfoVecs.nodes2Triangles_12 = k;
                } */

            }
        }

}
void Utilities::triangles2Triangles_host_vecs(
    int elem,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
    AuxVecs& auxVecs){
        if (hostSetInfoVecs.triangles2Edges_1[elem] >= (INT_MAX-1000) ||
            hostSetInfoVecs.triangles2Edges_1[elem] <= (-INT_MAX+1000) ||
            hostSetInfoVecs.triangles2Edges_1[elem] < 0 ||
            hostSetInfoVecs.triangles2Edges_2[elem] >= (INT_MAX-1000)  ||
            hostSetInfoVecs.triangles2Edges_2[elem] <= (-INT_MAX+1000) ||
            hostSetInfoVecs.triangles2Edges_2[elem] < 0 || 
            hostSetInfoVecs.triangles2Edges_3[elem] >= (INT_MAX-1000) || // nav changed these to _3 from _1 03/10/2025
            hostSetInfoVecs.triangles2Edges_3[elem] <= (-INT_MAX+1000)|| // nav changed these to _3 from _1 03/10/2025
            hostSetInfoVecs.triangles2Edges_3[elem] < 0 // nav changed these to _3 from _1 03/10/2025
        ){
            hostSetInfoVecs.triangles2Triangles_1[elem] = -1.0;//-INT_MAX;
            hostSetInfoVecs.triangles2Triangles_2[elem] = -1.0;//-INT_MAX;
            hostSetInfoVecs.triangles2Triangles_3[elem] = -1.0;//-INT_MAX;
        }
        else{

            int edge1_n1 = hostSetInfoVecs.triangles2Nodes_1[elem];
            int edge1_n2 = hostSetInfoVecs.triangles2Nodes_2[elem];
            int edge2_n1 = hostSetInfoVecs.triangles2Nodes_2[elem];
            int edge2_n2 = hostSetInfoVecs.triangles2Nodes_3[elem];
            int edge3_n1 = hostSetInfoVecs.triangles2Nodes_3[elem];
            int edge3_n2 = hostSetInfoVecs.triangles2Nodes_1[elem];
            
            if(1 == elem) {
                std::cout<<"edge1_n = ("<<edge1_n1<<","<<edge1_n2<<")"<<std::endl;
                std::cout<<"edge2_n = ("<<edge2_n1<<","<<edge2_n2<<")"<<std::endl;
                std::cout<<"edge3_n = ("<<edge3_n1<<","<<edge3_n2<<")"<<std::endl;
            }
            
            if (edge1_n1 != edge3_n2){
                std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING UP THE NODE ORDER OF EDGES!"<<std::endl;
            }
            int edge;
            for (int i = 0; i < 3; i++){
                //std::cout<<"i = "<<i<<std::endl;
                if (i == 0){
                    edge = hostSetInfoVecs.triangles2Edges_1[elem];
                }
                else if (i == 1){
                    edge = hostSetInfoVecs.triangles2Edges_2[elem];
                }
                else if (i == 2){
                    edge = hostSetInfoVecs.triangles2Edges_3[elem];
                }
                
                if(1 == elem) {
                    std::cout<<"tri = "<<elem<<"; i = "<<i<<"; edge = "<<edge<<std::endl;
                    //std::cout<<"edge1_n1="<<edge1_n1<<"; edge1_n2 = "<<edge1_n2<<std::endl; 
                    std::cout<<"edges2Tri_1 = "<<hostSetInfoVecs.edges2Triangles_1[edge]<<std::endl;
                    std::cout<<"edges2Tri_2 = "<<hostSetInfoVecs.edges2Triangles_2[edge]<<std::endl;
                    //std::cout<<"edges2Tri_3 = "<<hostSetInfoVecs.edges2Triangles_3[edge]<<std::endl;
                    std::cout<<std::endl;
                }                
                if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge1_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge1_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge1_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge1_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_1[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_1[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                            std::cout<<"error at here, exit(-1)"<<std::endl;
                            exit(-1);
                            
                            //return;//nav changed this 10/26/24 - this was not it. so commenting it out. 
                        }
                    }
                else if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge2_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge2_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge2_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge2_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_2[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_2[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                            std::cout<<"error at second place, exit(-1)"<<std::endl;
                            exit(-1);
                            
                        }
                    }
                else if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge3_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge3_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge3_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge3_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_3[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_3[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                        }
                    }
                else{
                    std::cout<<"error: triangles2Triangles_host_vecs() => edges2Nodes info did not agree with any nodes on the triangles listing"<<std::endl;
                }
            }
 /*           if (hostSetInfoVecs.triangles2Triangles_1[elem] == hostSetInfoVecs.triangles2Triangles_2[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            }
            else if (hostSetInfoVecs.triangles2Triangles_1[elem] == hostSetInfoVecs.triangles2Triangles_3[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            }
            else if (hostSetInfoVecs.triangles2Triangles_2[elem] == hostSetInfoVecs.triangles2Triangles_3[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            } */
        }
       // std::cout<<"triangles2Triangles["<<elem<<"] -> "<<hostSetInfoVecs.triangles2Triangles_1[elem]<<" "<<hostSetInfoVecs.triangles2Triangles_2[elem]<<" "<<hostSetInfoVecs.triangles2Triangles_3[elem]<<std::endl;
    }



//copy configuration from device to host
void Utilities::transferDtoH(GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end(),hostSetInfoVecs.nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end(),hostSetInfoVecs.nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end(),hostSetInfoVecs.nodeLocZ.begin());

    thrust::copy(generalParams.nodes_in_upperhem.begin(),generalParams.nodes_in_upperhem.end(),hostSetInfoVecs.nodes_in_upperhem.begin());
    thrust::copy(generalParams.triangles_in_upperhem.begin(),generalParams.triangles_in_upperhem.end(),hostSetInfoVecs.triangles_in_upperhem.begin());
    thrust::copy(generalParams.edges_in_upperhem.begin(),generalParams.edges_in_upperhem.end(),hostSetInfoVecs.edges_in_upperhem.begin());
    thrust::copy(generalParams.edges_in_upperhem_list.begin(),generalParams.edges_in_upperhem_list.end(),hostSetInfoVecs.edges_in_upperhem_list.begin());
    thrust::copy(generalParams.boundaries_in_upperhem.begin(),generalParams.boundaries_in_upperhem.end(),hostSetInfoVecs.boundaries_in_upperhem.begin());

    thrust::copy(coordInfoVecs.triangles2Nodes_1.begin(),coordInfoVecs.triangles2Nodes_1.end(),hostSetInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_2.begin(),coordInfoVecs.triangles2Nodes_2.end(),hostSetInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_3.begin(),coordInfoVecs.triangles2Nodes_3.end(),hostSetInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(coordInfoVecs.edges2Nodes_1.begin(),coordInfoVecs.edges2Nodes_1.end(),hostSetInfoVecs.edges2Nodes_1.begin());
    thrust::copy(coordInfoVecs.edges2Nodes_2.begin(),coordInfoVecs.edges2Nodes_2.end(),hostSetInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(coordInfoVecs.edges2Triangles_1.begin(),coordInfoVecs.edges2Triangles_1.end(),hostSetInfoVecs.edges2Triangles_1.begin());
    thrust::copy(coordInfoVecs.edges2Triangles_2.begin(),coordInfoVecs.edges2Triangles_2.end(),hostSetInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(coordInfoVecs.triangles2Edges_1.begin(),coordInfoVecs.triangles2Edges_1.end(),hostSetInfoVecs.triangles2Edges_1.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_2.begin(),coordInfoVecs.triangles2Edges_2.end(),hostSetInfoVecs.triangles2Edges_2.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_3.begin(),coordInfoVecs.triangles2Edges_3.end(),hostSetInfoVecs.triangles2Edges_3.begin());

    thrust::copy(coordInfoVecs.nndata1.begin(),coordInfoVecs.nndata1.end(),hostSetInfoVecs.nndata1.begin());
    thrust::copy(coordInfoVecs.nndata2.begin(),coordInfoVecs.nndata2.end(),hostSetInfoVecs.nndata2.begin());
    thrust::copy(coordInfoVecs.nndata3.begin(),coordInfoVecs.nndata3.end(),hostSetInfoVecs.nndata3.begin());
    thrust::copy(coordInfoVecs.nndata4.begin(),coordInfoVecs.nndata4.end(),hostSetInfoVecs.nndata4.begin());
    thrust::copy(coordInfoVecs.nndata5.begin(),coordInfoVecs.nndata5.end(),hostSetInfoVecs.nndata5.begin());
    thrust::copy(coordInfoVecs.nndata6.begin(),coordInfoVecs.nndata6.end(),hostSetInfoVecs.nndata6.begin());
    thrust::copy(coordInfoVecs.nndata7.begin(),coordInfoVecs.nndata7.end(),hostSetInfoVecs.nndata7.begin());
    thrust::copy(coordInfoVecs.nndata8.begin(),coordInfoVecs.nndata8.end(),hostSetInfoVecs.nndata8.begin());
    thrust::copy(coordInfoVecs.nndata9.begin(),coordInfoVecs.nndata9.end(),hostSetInfoVecs.nndata9.begin());
    //thrust::copy(coordInfoVecs.nndata10.begin(),coordInfoVecs.nndata10.end(),hostSetInfoVecs.nndata10.begin());
    //thrust::copy(coordInfoVecs.nndata11.begin(),coordInfoVecs.nndata11.end(),hostSetInfoVecs.nndata11.begin());
    //thrust::copy(coordInfoVecs.nndata12.begin(),coordInfoVecs.nndata12.end(),hostSetInfoVecs.nndata12.begin());

    thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(),coordInfoVecs.nodes2Triangles_1.end(),hostSetInfoVecs.nodes2Triangles_1.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(),coordInfoVecs.nodes2Triangles_2.end(),hostSetInfoVecs.nodes2Triangles_2.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(),coordInfoVecs.nodes2Triangles_3.end(),hostSetInfoVecs.nodes2Triangles_3.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(),coordInfoVecs.nodes2Triangles_4.end(),hostSetInfoVecs.nodes2Triangles_4.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(),coordInfoVecs.nodes2Triangles_5.end(),hostSetInfoVecs.nodes2Triangles_5.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(),coordInfoVecs.nodes2Triangles_6.end(),hostSetInfoVecs.nodes2Triangles_6.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(),coordInfoVecs.nodes2Triangles_7.end(),hostSetInfoVecs.nodes2Triangles_7.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(),coordInfoVecs.nodes2Triangles_8.end(),hostSetInfoVecs.nodes2Triangles_8.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(),coordInfoVecs.nodes2Triangles_9.end(),hostSetInfoVecs.nodes2Triangles_9.begin());

    thrust::copy(coordInfoVecs.triangles2Triangles_1.begin(),coordInfoVecs.triangles2Triangles_1.end(),hostSetInfoVecs.triangles2Triangles_1.begin());
    thrust::copy(coordInfoVecs.triangles2Triangles_2.begin(),coordInfoVecs.triangles2Triangles_2.end(),hostSetInfoVecs.triangles2Triangles_2.begin());
    thrust::copy(coordInfoVecs.triangles2Triangles_3.begin(),coordInfoVecs.triangles2Triangles_3.end(),hostSetInfoVecs.triangles2Triangles_3.begin());
};

//copy configuration from host to device
void Utilities::transferHtoD(GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(),hostSetInfoVecs.nodeLocX.end(),coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(),hostSetInfoVecs.nodeLocY.end(),coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(),hostSetInfoVecs.nodeLocZ.end(),coordInfoVecs.nodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.nodes_in_upperhem.begin(),hostSetInfoVecs.nodes_in_upperhem.end(),generalParams.nodes_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.triangles_in_upperhem.begin(),hostSetInfoVecs.triangles_in_upperhem.end(),generalParams.triangles_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.edges_in_upperhem.begin(),hostSetInfoVecs.edges_in_upperhem.end(),generalParams.edges_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.edges_in_upperhem_list.begin(),hostSetInfoVecs.edges_in_upperhem_list.end(),generalParams.edges_in_upperhem_list.begin());
    thrust::copy(hostSetInfoVecs.boundaries_in_upperhem.begin(),hostSetInfoVecs.boundaries_in_upperhem.end(),generalParams.boundaries_in_upperhem.begin());
    
    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(),hostSetInfoVecs.triangles2Nodes_1.end(),coordInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(),hostSetInfoVecs.triangles2Nodes_2.end(),coordInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(),hostSetInfoVecs.triangles2Nodes_3.end(),coordInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(),hostSetInfoVecs.edges2Nodes_1.end(),coordInfoVecs.edges2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(),hostSetInfoVecs.edges2Nodes_2.end(),coordInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(),hostSetInfoVecs.edges2Triangles_1.end(),coordInfoVecs.edges2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(),hostSetInfoVecs.edges2Triangles_2.end(),coordInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(),hostSetInfoVecs.triangles2Edges_1.end(),coordInfoVecs.triangles2Edges_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(),hostSetInfoVecs.triangles2Edges_2.end(),coordInfoVecs.triangles2Edges_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(),hostSetInfoVecs.triangles2Edges_3.end(),coordInfoVecs.triangles2Edges_3.begin());

    thrust::copy(hostSetInfoVecs.nndata1.begin(),hostSetInfoVecs.nndata1.end(),coordInfoVecs.nndata1.begin());
    thrust::copy(hostSetInfoVecs.nndata2.begin(),hostSetInfoVecs.nndata2.end(),coordInfoVecs.nndata2.begin());
    thrust::copy(hostSetInfoVecs.nndata3.begin(),hostSetInfoVecs.nndata3.end(),coordInfoVecs.nndata3.begin());
    thrust::copy(hostSetInfoVecs.nndata4.begin(),hostSetInfoVecs.nndata4.end(),coordInfoVecs.nndata4.begin());
    thrust::copy(hostSetInfoVecs.nndata5.begin(),hostSetInfoVecs.nndata5.end(),coordInfoVecs.nndata5.begin());
    thrust::copy(hostSetInfoVecs.nndata6.begin(),hostSetInfoVecs.nndata6.end(),coordInfoVecs.nndata6.begin());
    thrust::copy(hostSetInfoVecs.nndata7.begin(),hostSetInfoVecs.nndata7.end(),coordInfoVecs.nndata7.begin());
    thrust::copy(hostSetInfoVecs.nndata8.begin(),hostSetInfoVecs.nndata8.end(),coordInfoVecs.nndata8.begin());
    thrust::copy(hostSetInfoVecs.nndata9.begin(),hostSetInfoVecs.nndata9.end(),coordInfoVecs.nndata9.begin());
    //thrust::copy(hostSetInfoVecs.nndata10.begin(),hostSetInfoVecs.nndata10.end(),coordInfoVecs.nndata10.begin());
    //thrust::copy(hostSetInfoVecs.nndata11.begin(),hostSetInfoVecs.nndata11.end(),coordInfoVecs.nndata11.begin());
    //thrust::copy(hostSetInfoVecs.nndata12.begin(),hostSetInfoVecs.nndata12.end(),coordInfoVecs.nndata12.begin());

    thrust::copy(hostSetInfoVecs.nodes2Triangles_1.begin(),hostSetInfoVecs.nodes2Triangles_1.end(),coordInfoVecs.nodes2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_2.begin(),hostSetInfoVecs.nodes2Triangles_2.end(),coordInfoVecs.nodes2Triangles_2.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_3.begin(),hostSetInfoVecs.nodes2Triangles_3.end(),coordInfoVecs.nodes2Triangles_3.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_4.begin(),hostSetInfoVecs.nodes2Triangles_4.end(),coordInfoVecs.nodes2Triangles_4.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_5.begin(),hostSetInfoVecs.nodes2Triangles_5.end(),coordInfoVecs.nodes2Triangles_5.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_6.begin(),hostSetInfoVecs.nodes2Triangles_6.end(),coordInfoVecs.nodes2Triangles_6.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_7.begin(),hostSetInfoVecs.nodes2Triangles_7.end(),coordInfoVecs.nodes2Triangles_7.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_8.begin(),hostSetInfoVecs.nodes2Triangles_8.end(),coordInfoVecs.nodes2Triangles_8.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_9.begin(),hostSetInfoVecs.nodes2Triangles_9.end(),coordInfoVecs.nodes2Triangles_9.begin());

    thrust::copy(hostSetInfoVecs.triangles2Triangles_1.begin(),hostSetInfoVecs.triangles2Triangles_1.end(),coordInfoVecs.triangles2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Triangles_2.begin(),hostSetInfoVecs.triangles2Triangles_2.end(),coordInfoVecs.triangles2Triangles_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Triangles_3.begin(),hostSetInfoVecs.triangles2Triangles_3.end(),coordInfoVecs.triangles2Triangles_3.begin());
};

