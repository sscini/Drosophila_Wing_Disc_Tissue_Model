#include "System.h"
#include <random>
#include "Utilities.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
// #include "cdecs.h"
// #include "vmalloc.h"
// #include "cleanup.c"
// #include "vectormalloc.c"
// #include "screen.c"
// #include "ppsub.c"
// #include "error.c"
// #include "debug.c"
// #include "output.c"
// #include "simpleio.c"
// #include "times.c"
// #include "other.c"
// #include "machine.c"
// #include "fgetstrin.c"
//#include "General.h"
//#include <gintrp>
// #include "gas/gdecs/gdecs.h"
// #include "driver/ddecs.h"
#include <omp.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
//
///* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// *      CSV  ?  std::vector<StrainKeyframe>
// * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
//std::size_t
//readStrainScheduleCSV(const std::string&           csvPath,
//                      std::vector<StrainKeyframe>& keyframes)
//{
//    std::ifstream fp(csvPath);
//    if(!fp.good())
//        throw std::runtime_error("Utilities::readStrainScheduleCSV – cannot open "
//                                 + csvPath);
//
//    keyframes.clear();
//    std::string line;
//    bool firstLine = true;
//    while(std::getline(fp,line)) {
//        if(line.empty()) continue;
//
//        std::stringstream ss(line);
//        std::string tok;
//
//        // ignore header
//        if(firstLine) { firstLine=false; continue; }
//
//        StrainKeyframe kf{};
//        std::getline(ss,tok,',');  kf.t      = std::stod(tok);
//        std::getline(ss,tok,',');  kf.eps_r  = std::stod(tok);
//        std::getline(ss,tok,',');  kf.eps_t  = std::stod(tok);
//
//        keyframes.push_back(kf);
//    }
//    /* ensure chronological order (good habit if the CSV isn’t guaranteed) */
//    std::sort(keyframes.begin(), keyframes.end(),
//              [](const StrainKeyframe& a,const StrainKeyframe& b){return a.t<b.t;});
//    return keyframes.size();
//}


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
            
    /*int nnode = generalParams.maxNodeCount;
    std::vector<bool> boundary_node_temp(nnode,false);
    for (int i = 0; i < nnode; i++){
        if (coordInfoVecs.edges2Triangles_1[i] == coordInfoVecs.edges2Triangles_2[i]){
            boundary_node_temp[coordInfoVecs.edges2Nodes_1[i]] = true;
            boundary_node_temp[coordInfoVecs.edges2Nodes_2[i]] = true;
        }
    }
    
    //This creates a int vector whose length equals to number of nodes.
    //The initial mesh has every node paired with 6 neighboring nodes. 
    //During the simulation, the number will be updated accordingly. Therefore this has to be moved
    //to another location to avoid re-initialization every time Edgeswap is called.

    std::vector<int> nndata_temp(nnode, 6);
    
    boundary_node = boundary_node_temp;
    nndata = nndata_temp;*/
};

//This function aims to have a fixed range or neighborhood around the tip be weakened uniformly
void Utilities::gradient_weakening_update_host_vecs_tip(double sigma,
    //double max_height_index,
    double max_height_x,
    double max_height_y,
    double max_height_z,
    double distance_to_boundary,
    double distance_uniform_weak,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){

        double pi = 3.1415927;
    /* Scaling by gaussian distribution in the form of 
    (1/sqrt(2*pi*sigma^2))*Exp(-(d/|d|)^2/(sigma^2))/(1/sqrt(2*pi*sigma^2))
    */
   double scale;
   //double tip_threshold = 2.05*generalParams.Rmin;
   bool scale_need = false;
    if (pi < 0){
    //if (sigma == INT_MAX){
        //  for (int i = 0; i < coordInfoVecs.num_edges; i++){
        //     if (hostSetInfoVecs.edges_in_upperhem[i] != 1){
        //         hostSetInfoVecs.scaling_per_edge[i] = 0.0;
        //         continue;
        //     }
        //     else{
        //         hostSetInfoVecs.scaling_per_edge[i] = 1.0;
        //     }
        //  }
    }
    else{
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            int v1 = hostSetInfoVecs.edges2Nodes_1[i];
            int v2 = hostSetInfoVecs.edges2Nodes_2[i];
            if (hostSetInfoVecs.edges2Nodes_1[i] > (INT_MAX-100) || hostSetInfoVecs.edges2Nodes_1[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            else if (hostSetInfoVecs.edges2Nodes_2[i] > (INT_MAX-100) || hostSetInfoVecs.edges2Nodes_2[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
            if (hostSetInfoVecs.edges_in_upperhem[i] != 1 && hostSetInfoVecs.edges_in_upperhem[i] != 0){
                if (avg_z < generalParams.boundary_z){
                    hostSetInfoVecs.scaling_per_edge[i] = 1.0;
                    continue;
                }
                else{
                    double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                    double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                    double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                    double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                      //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                        //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                    //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                    scale = ((dtt-distance_uniform_weak)/(distance_to_boundary-distance_uniform_weak));
                    
                    if ((distance_to_boundary-distance_uniform_weak)<0.0){
                        scale = 0.0;
                    }
                    else if (scale > 1.0){
                        scale = 1.0;
                    }
                    else if (scale < 0.0){
                        scale = 0.0;
                    }
                    else{}
                    scale_need = true;
                }
            }
           
            
           

            if (generalParams.edges_in_upperhem[i] == 1 || generalParams.edges_in_upperhem[i] == 0){//(generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                  //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                    //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                scale = ((dtt-distance_uniform_weak)/(distance_to_boundary-distance_uniform_weak));
                
                
                //if (dtt < tip_threshold){
                //    scale = 1.0;
                //}
                //else {
                //    scale = ((dtt - tip_threshold)/(distance_to_boundary - tip_threshold));
                //}
                if ((distance_to_boundary-distance_uniform_weak) < 0.0){
                    scale = 0.0;
                }
                else if (scale > 1.0){
                    scale = 1.0;
                }
                else if (scale < 0.0){
                    scale = 0.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
            // coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            /*else if (generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
               // scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = 1.0;//scale;
                scale_need = true;
            }*/
            else{
                //coordInfoVecs.scaling_per_edge[i] = 0.0;
            }

            if (scale_need == true){
                hostSetInfoVecs.scaling_per_edge[i] = scale;
            }
            else{
                hostSetInfoVecs.scaling_per_edge[i] = 1.0;
            }
        }
    }
};

//This function is for hill function type weakening purpose
void Utilities::gradient_weakening_update_host_vecs(double sigma,
    //double max_height_index,
    double max_height_x,
    double max_height_y,
    double max_height_z,
    double distance_to_boundary,
    double distance_to_boundary_max,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){

        double pi = 3.1415927;
    /* Scaling by gaussian distribution in the form of 
    (1/sqrt(2*pi*sigma^2))*Exp(-(d/|d|)^2/(sigma^2))/(1/sqrt(2*pi*sigma^2))
    */
   double scale;
   //double tip_threshold = 2.05*generalParams.Rmin;
   bool scale_need = false;
    if (pi < 0){
    //if (sigma == INT_MAX){
        //  for (int i = 0; i < coordInfoVecs.num_edges; i++){
        //     if (hostSetInfoVecs.edges_in_upperhem[i] != 1){
        //         hostSetInfoVecs.scaling_per_edge[i] = 0.0;
        //         continue;
        //     }
        //     else{
        //         hostSetInfoVecs.scaling_per_edge[i] = 1.0;
        //     }
        //  }
    }
    else{
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            int v1 = hostSetInfoVecs.edges2Nodes_1[i];
            int v2 = hostSetInfoVecs.edges2Nodes_2[i];
            if (hostSetInfoVecs.edges2Nodes_1[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_1[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            else if (hostSetInfoVecs.edges2Nodes_2[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_2[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
            if (hostSetInfoVecs.edges_in_upperhem[i] != 1 || hostSetInfoVecs.edges_in_upperhem[i] != 0){
                if (avg_z < generalParams.boundary_z){
                    hostSetInfoVecs.scaling_per_edge[i] = 1.0;
                    continue;
                }
                else{
                    double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                    double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                    double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                    double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                      //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                        //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                    //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                    scale = (dtt/distance_to_boundary_max);
                    
                    if (scale > 1.0){
                        scale = 1.0;
                    }
                    else{}
                    scale_need = true;
                }
            }
           
            
           

            if (generalParams.edges_in_upperhem[i] == 1 || generalParams.edges_in_upperhem[i] == 0){//(generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                  //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                    //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                scale = (dtt/distance_to_boundary_max);
                
                
                //if (dtt < tip_threshold){
                //    scale = 1.0;
                //}
                //else {
                //    scale = ((dtt - tip_threshold)/(distance_to_boundary - tip_threshold));
                //}

                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
            // coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            /*else if (generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
               // scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = 1.0;//scale;
                scale_need = true;
            }*/
            else{
                //coordInfoVecs.scaling_per_edge[i] = 0.0;
            }

            if (scale_need == true){
                hostSetInfoVecs.scaling_per_edge[i] = scale;
            }
            else{
                hostSetInfoVecs.scaling_per_edge[i] = 1.0;
            }
        }
    }
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

    /*double UN = sqrt((current_normalX*current_normalX) +
                     (current_normalY*current_normalY) + 
                     (current_normalZ*current_normalZ));
    current_normalX = current_normalX/UN;
    current_normalY = current_normalY/UN;
    current_normalZ = current_normalZ/UN;*/

    /*for (int j = 0; j < AREA.size(); j++){
        coordInfoVecs.nodeForceX[inode] += generalParams.volume_spring_constant*current_normalX;
        coordInfoVecs.nodeForceY[inode] += generalParams.volume_spring_constant*current_normalY;
        coordInfoVecs.nodeForceZ[inode] += generalParams.volume_spring_constant*current_normalZ;
    }*/
};

double Utilities::find_suitable_location_to_grow(
    int iedge, 
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs
){
    int alpha = -1;
        
    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
    double vol_0, vol_1;
    double P0x_vol1, P0y_vol1, P0z_vol1, P0x_vol2, P0y_vol2, P0z_vol2;
    double N1x_vol, N1y_vol, N1z_vol, N2x_vol, N2y_vol, N2z_vol;
    bool GROWTH_ACCEPTED = false;
    double avg_area;
      
    if ( hostSetInfoVecs.edges2Triangles_1[iedge] != hostSetInfoVecs.edges2Triangles_2[iedge]){
        H0 = hostSetInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        //std::cout<<"H0 = "<<H0<<std::endl;
        T0 = hostSetInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        //std::cout<<"T0 = "<<T0<<std::endl;
        edge_start = hostSetInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = hostSetInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = hostSetInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        //std::cout<<"a1 = "<<a1<<std::endl;
        b1 = hostSetInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        //std::cout<<"b1 = "<<b1<<std::endl;
        c1 = hostSetInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        //std::cout<<"c1 = "<<c1<<std::endl;
        
        a2 = hostSetInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        //std::cout<<"a2 = "<<a2<<std::endl;
        b2 = hostSetInfoVecs.triangles2Edges_2[T0];
        //std::cout<<"b2 = "<<b2<<std::endl;
        c2 = hostSetInfoVecs.triangles2Edges_3[T0];
        //std::cout<<"c2 = "<<c2<<std::endl;
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //std::cout<<"H1 = "<<H1<<std::endl;
        //std::cout<<"H2 = "<<H2<<std::endl;
        //std::cout<<"T1 = "<<T1<<std::endl;
        //std::cout<<"T2 = "<<T2<<std::endl;

        //Now search for the associated 

        int CANDIDATE1_1 = hostSetInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = hostSetInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = hostSetInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = hostSetInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = hostSetInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = hostSetInfoVecs.triangles2Nodes_3[T0];
        
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}

        bool BAD_CHOICE = false;
        for (int q = 0; q < 4; q++){
            double qq;
            if (q == 0){
                qq = edge_start;
            }
            else if (q == 1){
                qq = edge_end;
            }
            else if (q == 2){
                qq = HEAD;
            }
            else if (q == 3){
                qq = TAIL;
            }
            int safe_flip1 = 0;
            if (hostSetInfoVecs.nndata1[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata2[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata3[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata4[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata5[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata6[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata7[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata8[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata9[qq] >= 0){safe_flip1 += 1;        }

            if (q == 0){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
                //std::cout<<"SAFE_FLIP_start = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 1){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_end = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 2 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_H = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 3 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_T = "<<safe_flip1<<std::endl;
                break;
            }
        }
        

        if (BAD_CHOICE == false){//(safe_flip1 < generalParams.safeguardthreshold && safe_flip2 < generalParams.safeguardthreshold){

            int H0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[H0];
            int H0n2 = edge_end;//hostSetInfoVecs.triangles2Nodes_2[H0];
            int H0n3 = HEAD;//hostSetInfoVecs.triangles2Nodes_3[H0];
            int T0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[T0];
            int T0n2 = TAIL;//hostSetInfoVecs.triangles2Nodes_2[T0];
            int T0n3 = edge_end;//hostSetInfoVecs.triangles2Nodes_3[T0];
            double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2.0;
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            avg_area = (area_H0 + area_T0)/2.0;
            if ((avg_area - areaTriangleInfoVecs.initial_area)/areaTriangleInfoVecs.initial_area >= generalParams.strain_threshold){
            GROWTH_ACCEPTED = true;
            }
        }
    }
    if (GROWTH_ACCEPTED == true){
        return (avg_area - areaTriangleInfoVecs.initial_area)/areaTriangleInfoVecs.initial_area;
    }
    else{
        return -1.0;
    }
}

// int Utilities::growth_host_vecs(
int Utilities::growth_host_vecs(
    int iedge, 
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs
){
    int alpha;
    // std::vector<int> alpha;

    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
    double vol_0, vol_1;
    double P0x_vol1, P0y_vol1, P0z_vol1, P0x_vol2, P0y_vol2, P0z_vol2;
    double N1x_vol, N1y_vol, N1z_vol, N2x_vol, N2y_vol, N2z_vol;
    bool GROWTH_ACCEPTED = false;
      
    if ( hostSetInfoVecs.edges2Triangles_1[iedge] != hostSetInfoVecs.edges2Triangles_2[iedge]){
        H0 = hostSetInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        generalParams.triangle_undergoing_growth.push_back(H0);
        //std::cout<<"H0 = "<<H0<<std::endl;
        T0 = hostSetInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        generalParams.triangle_undergoing_growth.push_back(T0);
        //std::cout<<"T0 = "<<T0<<std::endl;
        edge_start = hostSetInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = hostSetInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = hostSetInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        //std::cout<<"a1 = "<<a1<<std::endl;
        b1 = hostSetInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        //std::cout<<"b1 = "<<b1<<std::endl;
        c1 = hostSetInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        //std::cout<<"c1 = "<<c1<<std::endl;
        
        a2 = hostSetInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        //std::cout<<"a2 = "<<a2<<std::endl;
        b2 = hostSetInfoVecs.triangles2Edges_2[T0];
        //std::cout<<"b2 = "<<b2<<std::endl;
        c2 = hostSetInfoVecs.triangles2Edges_3[T0];
        //std::cout<<"c2 = "<<c2<<std::endl;
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //std::cout<<"H1 = "<<H1<<std::endl;
        //std::cout<<"H2 = "<<H2<<std::endl;
        //std::cout<<"T1 = "<<T1<<std::endl;
        //std::cout<<"T2 = "<<T2<<std::endl;

        //Now search for the associated 

        int CANDIDATE1_1 = hostSetInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = hostSetInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = hostSetInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = hostSetInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = hostSetInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = hostSetInfoVecs.triangles2Nodes_3[T0];
        
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}

        bool BAD_CHOICE = false;
        for (int q = 0; q < 4; q++){
            double qq;
            if (q == 0){
                qq = edge_start;
            }
            else if (q == 1){
                qq = edge_end;
            }
            else if (q == 2){
                qq = HEAD;
            }
            else if (q == 3){
                qq = TAIL;
            }
            
            int safe_flip1 = 0;
            if (hostSetInfoVecs.nndata1[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata2[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata3[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata4[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata5[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata6[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata7[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata8[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata9[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata10[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata11[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata12[qq] >= 0){safe_flip1 += 1;        }

            if (q == 0){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
                //std::cout<<"SAFE_FLIP_start = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 1){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_end = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 2 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_H = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 3 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_T = "<<safe_flip1<<std::endl;
                break;
            }
        }
        

        if (BAD_CHOICE == false){//(safe_flip1 < generalParams.safeguardthreshold && safe_flip2 < generalParams.safeguardthreshold){

            int H0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[H0];
            int H0n2 = edge_end;//hostSetInfoVecs.triangles2Nodes_2[H0];
            int H0n3 = HEAD;//hostSetInfoVecs.triangles2Nodes_3[H0];
            int T0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[T0];
            int T0n2 = TAIL;//hostSetInfoVecs.triangles2Nodes_2[T0];
            int T0n3 = edge_end;//hostSetInfoVecs.triangles2Nodes_3[T0];
            double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2.0;
            // double area_spring_constant_1, area_spring_constant_2;
            // if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){
            //     area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
            // }
            // else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){
            //     area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            // }
            // else{
            //     area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
            // }
            // if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){
            //     area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
            // }
            // else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){
            //     area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            // }
            // else{
            //     area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
            // }
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            double avg_area = (area_H0 + area_T0)/2.0;
            // std::random_device rd;  //Will be used to obtain a seed for the random number engine
            // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            // std::uniform_real_distribution<> dis(0.0, 1.0);
            // double random_number = dis(gen);
            // //double random_number = 0;
            // double Edif = generalParams.insertion_energy_cost*(1.0 - (1.0/generalParams.strain_threshold)*((avg_area/areaTriangleInfoVecs.initial_area) - 1));
            // //std::cout<<Edif<<std::endl;
            // double prob = generalParams.tau*exp(-(Edif*generalParams.growth_energy_scaling)/generalParams.kT);
            // if (prob >= 1){
            //     prob = 1;
            // }
            // //std::cout<<"prob = "<<prob<<std::endl;
            // if (random_number < prob){
            //     GROWTH_ACCEPTED = true;
            // }
            if ((avg_area - areaTriangleInfoVecs.initial_area)/areaTriangleInfoVecs.initial_area >= generalParams.strain_threshold){
            GROWTH_ACCEPTED = true;
        }
      }
    //if (H0 > 20000.0){
      if (GROWTH_ACCEPTED == true){
        /*  int k = iedge;
		if (hostSetInfoVecs.edges2Nodes_1[k] == INT_MAX || hostSetInfoVecs.edges2Nodes_2[k] == INT_MAX){
			continue;
		}
		if (generalParams.boundaries_in_upperhem[k] == 1){
			continue;
		}
		if (generalParams.edges_in_upperhem[k] < 0 || generalParams.edges_in_upperhem[k] == INT_MAX){
			continue;
		}
		int iedge = k;*/
		
		
		////std::cout<<"GROWTH ERROR 2"<<std::endl;	
		int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;
        int elem1 = H0;
        int elem2 = T0;
		if (hostSetInfoVecs.triangles2Edges_1[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_2[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_3[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_1[elem1];
		}
		else if (hostSetInfoVecs.triangles2Edges_2[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_3[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_1[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_2[elem1];
		} 
		else if (hostSetInfoVecs.triangles2Edges_3[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_1[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_2[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_3[elem1];
		}
		////std::cout<<"GROWTH ERROR 3"<<std::endl;	

		if (hostSetInfoVecs.triangles2Edges_1[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_2[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_3[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_1[elem2];
		}
		else if (hostSetInfoVecs.triangles2Edges_2[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_3[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_1[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_2[elem2];
		} 
		else if (hostSetInfoVecs.triangles2Edges_3[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_1[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_2[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_3[elem2];
		}
		//Note that in the above assignment, t1e3 and t2e3 are the edges shared by the neighboring triangles T1 and T2.
		////std::cout<<"GROWTH ERROR 4"<<std::endl;	

		//VectorShuffleForEdgeswapLoop.push_back(t1e1);
		//VectorShuffleForEdgeswapLoop.push_back(t1e2);
		//VectorShuffleForEdgeswapLoop.push_back(t2e1);
		//VectorShuffleForEdgeswapLoop.push_back(t2e2);
		int n1, n2, n3, n4;
		
		if ((hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]) || (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]) ){
			n1 = hostSetInfoVecs.edges2Nodes_1[t1e1];
			n2 = hostSetInfoVecs.edges2Nodes_2[t1e1];
			if (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_2[iedge];
			}
			else if (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_1[iedge];
			}
		}
		else if ((hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]) || (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]) ){
			n1 = hostSetInfoVecs.edges2Nodes_2[t1e1];
			n2 = hostSetInfoVecs.edges2Nodes_1[t1e1];
			if (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_2[iedge];
			}
			else if (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_1[iedge];
			}
		}
		////std::cout<<"GROWTH ERROR 5"<<std::endl;	

		if (hostSetInfoVecs.edges2Nodes_1[t2e1] == hostSetInfoVecs.edges2Nodes_1[iedge] || hostSetInfoVecs.edges2Nodes_1[t2e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
			n4 = hostSetInfoVecs.edges2Nodes_2[t2e1];
		}
		else if (hostSetInfoVecs.edges2Nodes_2[t2e1] == hostSetInfoVecs.edges2Nodes_1[iedge] || hostSetInfoVecs.edges2Nodes_2[t2e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
			n4 = hostSetInfoVecs.edges2Nodes_1[t2e1];
		}
		

		//std::cout<<"n1 = "<<n1<<std::endl;
		//std::cout<<"n2 = "<<n2<<std::endl;
		//std::cout<<"n3 = "<<n3<<std::endl;
		//std::cout<<"n4 = "<<n4<<std::endl;
		//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).

		////std::cout<<"GROWTH ERROR 6"<<std::endl;	
		int edgeindex, a, a1, a2, a3, temp1, temp2;
		//std::cout<<"maxNodeCount = "<< generalParams.maxNodeCount<<std::endl;
		double newx = (hostSetInfoVecs.nodeLocX[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocX[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.1"<<std::endl;	
		hostSetInfoVecs.nodeLocX[generalParams. maxNodeCount] = newx;
		////std::cout<<"GROWTH ERROR 6.2"<<std::endl;	
		double newy = (hostSetInfoVecs.nodeLocY[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocY[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.3"<<std::endl;	
		hostSetInfoVecs.nodeLocY[generalParams. maxNodeCount] = newy;
		////std::cout<<"GROWTH ERROR 6.4"<<std::endl;	
		double newz = (hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.5"<<std::endl;	
		hostSetInfoVecs.nodeLocZ[generalParams. maxNodeCount] = newz;
		//These are the coordinate of the new vertex. Its index is "hostSetInfoVecs.nodeLocX.size()-1"

		//Before editing major data structures, we will update the nndata here since it is only affected by the addition of new nodes.

		//int NODESIZE= generalParams.maxNodeCount;//hostSetInfoVecs.nodeLocX.size();
		////std::cout<<"GROWTH ERROR 7"<<std::endl;			
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] = n1;
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] = n2;
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//NOTE: What this +1 actually does is that it specifies the location to write
		//any new data. Here it points to the location to write new triangles information.
		//This is a new triangle associated with (tn1, tn2, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-4".
		////std::cout<<"GROWTH ERROR 8"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n2);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n3);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn2, tn3, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-3".
		////std::cout<<"GROWTH ERROR 9"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n3);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n4);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-2".
		////std::cout<<"GROWTH ERROR 10"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n4);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n1);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-1".
		////std::cout<<"GROWTH ERROR 11"<<std::endl;	
		//Now we add new edges formed by the addition of the new node.
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n1);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 12"<<std::endl;	
		//This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-4".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n2);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 13"<<std::endl;	
		//This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-3".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n3);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 14"<<std::endl;	
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-2".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n4);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 15"<<std::endl;	
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 16"<<std::endl;				
			//Now we check to see if the order of update is correct, i.e. are edges2Triangles data in correct orientation.
			//This is crucial in the bendingspring computation.
			edgeindex = (coordInfoVecs.num_edges - (4-j));
			a = hostSetInfoVecs.edges2Triangles_1[edgeindex];
			if ((hostSetInfoVecs.triangles2Nodes_1[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_2[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a1 = 1;
			}
			else{
				a1 = 0;
			}
			if ((hostSetInfoVecs.triangles2Nodes_2[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_3[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a2 = 1;
			}
			else{
				a2 = 0;
			}
			if ((hostSetInfoVecs.triangles2Nodes_3[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_1[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a3 = 1;
			}
			else{
				a3 = 0;
			}

			if ((a1+a2+a3) == 0){
				temp1 = hostSetInfoVecs.edges2Triangles_1[edgeindex];
				temp2 = hostSetInfoVecs.edges2Triangles_2[edgeindex];
				hostSetInfoVecs.edges2Triangles_1[edgeindex] = temp2;
				hostSetInfoVecs.edges2Triangles_2[edgeindex] = temp1;
			}
			else{}
			//This checks if the orientation is correct or not, if not, flip the ordering.
		}
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-1".
		generalParams.maxNodeCount += 1;

		hostSetInfoVecs.nndata1[generalParams.maxNodeCount-1] =  (n1);
		hostSetInfoVecs.nndata2[generalParams.maxNodeCount-1] =  (n2);
		hostSetInfoVecs.nndata3[generalParams.maxNodeCount-1] =  (n3);
		hostSetInfoVecs.nndata4[generalParams.maxNodeCount-1] =  (n4);
		hostSetInfoVecs.nndata5[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata6[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata7[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata8[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata9[generalParams.maxNodeCount-1] =  (-2);
		//hostSetInfoVecs.nndata10[generalParams.maxNodeCount-1] = (-2);
		//hostSetInfoVecs.nndata11[generalParams.maxNodeCount-1] = (-2);
		//hostSetInfoVecs.nndata12[generalParams.maxNodeCount-1] = (-2);
        hostSetInfoVecs.nodes2Triangles_1[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 4);
        hostSetInfoVecs.nodes2Triangles_2[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 3);
        hostSetInfoVecs.nodes2Triangles_3[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 2);
        hostSetInfoVecs.nodes2Triangles_4[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 1);
        hostSetInfoVecs.nodes2Triangles_5[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_6[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_7[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_8[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_9[generalParams.maxNodeCount-1] =  (-INT_MAX);
		for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n1;
				nnn = n3;
				nnnn = generalParams.maxNodeCount-1;
			}
			else if (j == 1){
				nn = n3;
				nnn = n1;
				nnnn = generalParams.maxNodeCount-1;
			}
			if (hostSetInfoVecs.nndata1[nn] == nnn){
				hostSetInfoVecs.nndata1[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata2[nn] == nnn){
				hostSetInfoVecs.nndata2[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata3[nn] == nnn){
				hostSetInfoVecs.nndata3[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata4[nn] == nnn){
				hostSetInfoVecs.nndata4[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata5[nn] == nnn){
				hostSetInfoVecs.nndata5[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata6[nn] == nnn){
				hostSetInfoVecs.nndata6[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata7[nn] == nnn){
				hostSetInfoVecs.nndata7[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata8[nn] == nnn){
				hostSetInfoVecs.nndata8[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata9[nn] == nnn){
				hostSetInfoVecs.nndata9[nn] = nnnn;
			}
			/*else if (hostSetInfoVecs.nndata10[nn] == nnn){
				hostSetInfoVecs.nndata10[nn] = nnnn;
			}
			 else if (hostSetInfoVecs.nndata11[nn] == nnn){
				hostSetInfoVecs.nndata11[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata12[nn] == nnn){
				hostSetInfoVecs.nndata12[nn] = nnnn;
			} */
		}

        for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n1;
			}
			else if (j == 1){
				nn = n3;
			}

			nnn = elem1;
			nnnn = elem2;

			if (hostSetInfoVecs.nodes2Triangles_1[nn] == nnn || hostSetInfoVecs.nodes2Triangles_1[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_2[nn] == nnn || hostSetInfoVecs.nodes2Triangles_2[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_3[nn] == nnn || hostSetInfoVecs.nodes2Triangles_3[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_4[nn] == nnn || hostSetInfoVecs.nodes2Triangles_4[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_5[nn] == nnn || hostSetInfoVecs.nodes2Triangles_5[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_6[nn] == nnn || hostSetInfoVecs.nodes2Triangles_6[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_7[nn] == nnn || hostSetInfoVecs.nodes2Triangles_7[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_8[nn] == nnn || hostSetInfoVecs.nodes2Triangles_8[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_9[nn] == nnn || hostSetInfoVecs.nodes2Triangles_9[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
			}

			for (int k = 0; k < 2; k++){
				if (j == 0 && k == 0){
					nnn = coordInfoVecs.num_triangles - 4;
				}
				else if (j == 0 && k == 1){
					nnn = coordInfoVecs.num_triangles - 1;
				}
				else if (j == 1 && k == 0){
					nnn = coordInfoVecs.num_triangles - 3;
				}
				else if (j == 1 && k == 1){
					nnn = coordInfoVecs.num_triangles - 2;
				}
				if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_1[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_2[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_2[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_3[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_4[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_5[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_6[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_6[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_7[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_7[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_8[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_9[nn] = nnn;
				}
			}
		}

		for (int j = 0; j < 2; j++){
			int nn, nnn;
			if (j == 0){
				nn = n2;
				nnn = generalParams.maxNodeCount-1;
			}
			else if (j == 1){
				nn = n4;
				nnn = generalParams.maxNodeCount-1;
			}
			if (hostSetInfoVecs.nndata1[nn] < 0){
				hostSetInfoVecs.nndata1[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata2[nn] < 0){
				hostSetInfoVecs.nndata2[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata3[nn] < 0){
				hostSetInfoVecs.nndata3[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata4[nn] < 0){
				hostSetInfoVecs.nndata4[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata5[nn] < 0){
				hostSetInfoVecs.nndata5[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata6[nn] < 0){
				hostSetInfoVecs.nndata6[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata7[nn] < 0){
				hostSetInfoVecs.nndata7[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata8[nn] < 0){
				hostSetInfoVecs.nndata8[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata9[nn] < 0){
				hostSetInfoVecs.nndata9[nn] = nnn;
			}
			/*else if (hostSetInfoVecs.nndata10[nn] < 0){
				hostSetInfoVecs.nndata10[nn] = nnn;
			}
			 else if (hostSetInfoVecs.nndata11[nn] < 0){
				hostSetInfoVecs.nndata11[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata12[nn] < 0){
				hostSetInfoVecs.nndata12[nn] = nnn;
			} */
		}

        for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n2;
				nnn = elem1;
                //nnnn = elem2;
			}
			else if (j == 1){
				nn = n4;
				nnn = elem2;
                //nnnn = elem1;
			}
			if (hostSetInfoVecs.nodes2Triangles_1[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_1[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_2[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_2[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_3[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_3[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_4[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_4[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_5[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_5[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_6[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_6[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_7[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_7[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_8[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_8[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_9[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_9[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
			}
		}
		
		for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n2;
			}
			else if (j == 1){
				nn = n4;
			}
			for (int k = 0; k < 2; k++){
				if (j == 0 && k == 0){
					nnn = coordInfoVecs.num_triangles - 4;
				}
				else if (j == 0 && k == 1){
					nnn = coordInfoVecs.num_triangles - 3;
				}
				else if (j == 1 && k == 0){
					nnn = coordInfoVecs.num_triangles - 2;
				}
				else if (j == 1 && k == 1){
					nnn = coordInfoVecs.num_triangles - 1;
				}
				if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_1[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_2[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_2[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_3[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_4[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_5[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_6[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_6[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_7[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_7[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_8[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_9[nn] = nnn;
				}
			}
		}
		//generalParams.num_of_nodes += 1;

		


		////std::cout<<"GROWTH ERROR 17"<<std::endl;	
		//Now we update the edges2Triangles data structure with new edges.
		//std::cout<<"elem 1 = "<<elem1<<std::endl;
		//std::cout<<"elem 2 = "<<elem2<<std::endl;
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
		}
		int TRIANGLESIZE = coordInfoVecs.num_triangles;//hostSetInfoVecs.triangles2Nodes_1.size();
		if (hostSetInfoVecs.edges2Triangles_1[t1e1] == elem1){
			hostSetInfoVecs.edges2Triangles_1[t1e1] = TRIANGLESIZE-4;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t1e1] == elem1){
			hostSetInfoVecs.edges2Triangles_2[t1e1] = TRIANGLESIZE-4;
		}
		else{}
		////std::cout<<"GROWTH ERROR 18"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t1e2] == elem1){
			hostSetInfoVecs.edges2Triangles_1[t1e2] = TRIANGLESIZE-3;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t1e2] == elem1){
			hostSetInfoVecs.edges2Triangles_2[t1e2] = TRIANGLESIZE-3;
		}
		else{}
		////std::cout<<"GROWTH ERROR 19"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t2e1] == elem2){
			hostSetInfoVecs.edges2Triangles_1[t2e1] = TRIANGLESIZE-2;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t2e1] == elem2){
			hostSetInfoVecs.edges2Triangles_2[t2e1] = TRIANGLESIZE-2;
		}
		else{}
		////std::cout<<"GROWTH ERROR 20"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t2e2] == elem2){
			hostSetInfoVecs.edges2Triangles_1[t2e2] = TRIANGLESIZE-1;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t2e2] == elem2){
			hostSetInfoVecs.edges2Triangles_2[t2e2] = TRIANGLESIZE-1;
		}
		else{}
		//std::cout<<"t1e1 "<<t1e1<<std::endl;
		//std::cout<<"t1e2 "<<t1e2<<std::endl;
		//std::cout<<"t1e3 "<<t1e3<<std::endl;
		//std::cout<<"t2e1 "<<t2e1<<std::endl;
		//std::cout<<"t2e2 "<<t2e2<<std::endl;
		//std::cout<<"t2e3 "<<t2e3<<std::endl;

		//for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
		//}
		//The above change the existing edges2Triangles data structure to accomodate new triangles added.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Now we will take care of the last unedited data structure "triangles2Edges".
		//int aa, bb;
		int EDGESIZE = coordInfoVecs.num_edges;//hostSetInfoVecs.edges2Nodes_1.size();
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 21"<<std::endl;	
			if (j == 0){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 4] = (EDGESIZE-4);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 4] = (t1e1);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 4] = (EDGESIZE-3);   
			}
			else if (j == 1){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 3] = (EDGESIZE-3);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 3] = (t1e2);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 3] = (EDGESIZE-2);   
			}
			else if (j ==2){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 2] = (EDGESIZE-2);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 2] = (t2e1);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 2] = (EDGESIZE-1);   
			}
			else if (j ==3){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 1] = (EDGESIZE-1);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 1] = (t2e2);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 1] = (EDGESIZE-4);   
			}
			
		}
	
 // nav commented the below lines on 04/20/25 because she wants to introduce the upperhem vector stuff in the XML parser and systemBuilder. 
		// nav look here. This is an interesting way to go about the computation. Will keep this scheme definitely. Now here they already have (somehow, the nodes_in_upperhem populated with upperhem stuff.
//		if (hostSetInfoVecs.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_1[iedge]] == 1 && generalParams.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_2[iedge]] == 1){
//			//generalParams.nodes_in_upperhem.push_back(1);
//            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = 1;
//
//		}
//		else if (hostSetInfoVecs.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_1[iedge]] == 1 || generalParams.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_2[iedge]] == 1){
//			//generalParams.nodes_in_upperhem.push_back(0);
//            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = 0;
//		}
//		else{
//			//generalParams.nodes_in_upperhem.push_back(-1);
//            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = -1;
	//	}
		//Finally, we will fill the edge data chosen for growth (expansion) with INT_MAX so its data is no longer relevant to the computation
		////std::cout<<"GROWTH ERROR 22"<<std::endl;	
		hostSetInfoVecs.edges2Nodes_1[iedge] = INT_MAX;
		hostSetInfoVecs.edges2Nodes_2[iedge] = INT_MAX;
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//	//std::cout<<"GROWTH ERROR 23"<<std::endl;	
			if (hostSetInfoVecs.triangles2Edges_1[i] == iedge){
				hostSetInfoVecs.triangles2Edges_1[i] = INT_MAX;
			}
			if (hostSetInfoVecs.triangles2Edges_2[i] == iedge){
				hostSetInfoVecs.triangles2Edges_2[i] = INT_MAX;
			}
			if (hostSetInfoVecs.triangles2Edges_3[i] == iedge){
				hostSetInfoVecs.triangles2Edges_3[i] = INT_MAX;
			}
		}
		hostSetInfoVecs.edges2Triangles_1[iedge] = INT_MAX;
		hostSetInfoVecs.edges2Triangles_2[iedge] = INT_MAX;
		
		////std::cout<<"GROWTH ERROR 24"<<std::endl;	
		
			hostSetInfoVecs.triangles2Nodes_1[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_2[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_3[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_1[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_2[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_3[elem2] = INT_MAX;
			
			//Delete the associated vertices information of the selected triangle.
			//Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
			//Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
			for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//		//std::cout<<"GROWTH ERROR 25"<<std::endl;	
				if (hostSetInfoVecs.edges2Triangles_1[i] == elem1 || hostSetInfoVecs.edges2Triangles_1[i] == elem2){
					hostSetInfoVecs.edges2Triangles_1[i] = INT_MAX;
				}
				if (hostSetInfoVecs.edges2Triangles_2[i] == elem1 || hostSetInfoVecs.edges2Triangles_2[i] == elem2){
					hostSetInfoVecs.edges2Triangles_2[i] = INT_MAX;
				}
			if (hostSetInfoVecs.edges2Triangles_1[i] != INT_MAX && hostSetInfoVecs.edges2Triangles_2[i] == INT_MAX){
				std::cout<<"modified edges2Triangles "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
				}
				else if (hostSetInfoVecs.edges2Triangles_1[i] == INT_MAX && hostSetInfoVecs.edges2Triangles_2[i] != INT_MAX){
					std::cout<<"modified edges2Triangles "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
					}
			}
			//This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
		//	//std::cout<<"GROWTH ERROR 26"<<std::endl;	
			hostSetInfoVecs.triangles2Edges_1[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_2[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_3[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_1[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_2[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_3[elem2] = INT_MAX;
			for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//		//std::cout<<"GROWTH ERROR 27"<<std::endl;	
				if (hostSetInfoVecs.triangles2Edges_1[i] == iedge){
					hostSetInfoVecs.triangles2Edges_1[i] = INT_MAX;
				}
				if (hostSetInfoVecs.triangles2Edges_2[i] == iedge ){
					hostSetInfoVecs.triangles2Edges_2[i] = INT_MAX;
				}
				if (hostSetInfoVecs.triangles2Edges_3[i] == iedge ){
					hostSetInfoVecs.triangles2Edges_3[i] = INT_MAX;
				}
			}
		
		//Erase the edge infomation related to the deleted triangle. Note the deletion should always start with the largest index.

		//Before we delete the edge, determine whether the newly added node is part of nodes_in_upperhem or not.
		

		
						//Erase the edge infomation related to the deleted triangle.

						//Now we update the nodes_in_upperhem and edges_in_upperhem data structures.
						//This ensures that the newly created edges will have the correct associated spring constant.
//std::cout<<"ERROR HERE?"<<std::endl;
		generalParams.edges_in_upperhem[iedge] = INT_MAX;
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
			if (generalParams.edges_in_upperhem_list[i] == iedge){
				generalParams.edges_in_upperhem_list[i] = INT_MAX;
				//break;
			}
		}

		for (int q = 0; q < 4; q++){
		//	//std::cout<<"GROWTH ERROR 30"<<std::endl;	
			int nodeP = hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles - (4-q)]; 
			int nodeQ = hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles - (4-q)];
			int nodeR = hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles - (4-q)];

			if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(1);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 1;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else{
				//generalParams.triangles_in_upperhem.push_back(INT_MAX);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = INT_MAX;
			}
		}
		//std::cout<<"edges2Triangles size"<<""<<hostSetInfoVecs.edges2Triangles_1.size()<<" "<<hostSetInfoVecs.edges2Triangles_2.size()<<std::endl;
		//std::cout<<"triangles_in_upperhem size "<<generalParams.triangles_in_upperhem.size()<<std::endl;	
		//std::cout<<"GROWTH ERROR 29"<<std::endl;	
		for (int q = 0; q < 4; q++){
			//std::cout<<"GROWTH ERROR 31"<<std::endl;	
			int elem_1 = hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<coordInfoVecs.num_edges-(4 - q)<<std::endl;
			//std::cout<<"elem_1 "<<elem_1<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeP]<<std::endl;
			int elem_2 = hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<"elem_2"<<elem_2<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeQ]<<std::endl;
			//std::cout<<"GROWTH ERROR 31.5"<<std::endl;
			if (hostSetInfoVecs.triangles_in_upperhem[elem_1] == 1 && hostSetInfoVecs.triangles_in_upperhem[elem_2] == 1){
				
				//generalParams.edges_in_upperhem.push_back(1);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = 1;
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
				//generalParams.edges_in_upperhem_list.push_back(coordInfoVecs.num_edges - (4 - q));
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = coordInfoVecs.num_edges - (4-q);
                generalParams.edges_in_upperhem_list_length += 1;
			}
			
			else if (generalParams.triangles_in_upperhem[elem_1] == 1 || generalParams.triangles_in_upperhem[elem_2] == 1){
				
				//generalParams.edges_in_upperhem.push_back(0);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = 0;
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = -INT_MAX;
                generalParams.edges_in_upperhem_list_length += 1;
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
			}
			else{
				
				//generalParams.edges_in_upperhem.push_back(-1);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = -1;
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = -INT_MAX;
                generalParams.edges_in_upperhem_list_length += 1;
			}

			//generalParams.boundaries_in_upperhem.push_back(-1);
            hostSetInfoVecs.boundaries_in_upperhem[coordInfoVecs.num_edges - (4-q)] = -1;		
			
		}
		generalParams.triangles_in_upperhem[elem1] = INT_MAX;
		generalParams.triangles_in_upperhem[elem2] = INT_MAX;
      }
    }  

    if (GROWTH_ACCEPTED == true){
        double midpoint_x = (hostSetInfoVecs.nodeLocX[edge_start] + hostSetInfoVecs.nodeLocX[edge_end])/2.0;
        // std::cout<<"wwwww"<<std::endl;
        double midpoint_y = (hostSetInfoVecs.nodeLocY[edge_start] + hostSetInfoVecs.nodeLocY[edge_end])/2.0;
        // std::cout<<"wwwwww"<<std::endl;
        double midpoint_z = (hostSetInfoVecs.nodeLocZ[edge_start] + hostSetInfoVecs.nodeLocZ[edge_end])/2.0;
        std::cout<<"midpoint of edge chosen for growth : "<<midpoint_x<<", "<<midpoint_y<<", "<<midpoint_z<<std::endl;
        alpha = iedge;//1;
    }
    else{
        // alpha = -1;
        // generalParams.triangle_undergoing_growth.clear();
        // alpha.push_back(-1);
        alpha = -1;
    }
    // return alpha;
    return alpha;
}

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

// void Utilities::triangles2Triangles_host_vecs(
//     int elem,
//     HostSetInfoVecs& hostSetInfoVecs,
//     CoordInfoVecs& coordInfoVecs,
// 	GeneralParams& generalParams,
//     AuxVecs& auxVecs){
//         if (hostSetInfoVecs.triangles2Edges_1[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] < 0 ||
//             hostSetInfoVecs.triangles2Edges_2[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_2[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_2[elem] < 0 ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] < 0
//         ){
//             hostSetInfoVecs.triangles2Triangles_1[elem] = -1.0;//-INT_MAX;
//             hostSetInfoVecs.triangles2Triangles_2[elem] = -1.0;//-INT_MAX;
//             hostSetInfoVecs.triangles2Triangles_3[elem] = -1.0;//-INT_MAX;
//         }
//         else{

//             int edge1 = hostSetInfoVecs.triangles2Edges_1[elem];
//             int edge2 = hostSetInfoVecs.triangles2Edges_2[elem];
//             int edge3 = hostSetInfoVecs.triangles2Edges_3[elem];
//             int edge;
//             for (int i = 0; i < 3; i++){
//                 //std::cout<<"i = "<<i<<std::endl;
//                 if (i == 0){
//                     edge = edge1;
//                 }
//                 else if (i == 1){
//                     edge = edge2;
//                 }
//                 else if (i == 2){
//                     edge = edge3;
//                 }
//                 int triangle1 = hostSetInfoVecs.edges2Triangles_1[edge];
//                 int triangle2 = hostSetInfoVecs.edges2Triangles_2[edge];
//                 if (i == 0 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_1[elem] = triangle2;
//                 }
//                 else if (i == 0 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle2<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_1[elem] = triangle1;
//                 }
//                 else if (i == 0 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }

//                 if (i == 1 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     //std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_2[elem] = triangle2;
//                 }
//                 else if (i == 1 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // //std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_2[elem] = triangle1;
//                 }
//                 else if (i == 1 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }

//                 if (i == 2 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     //std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_3[elem] = triangle2;
//                 }
//                 else if (i == 2 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // //std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_3[elem] = triangle1;
//                 }
//                 else if (i == 2 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }
//             }
//         }
//     }

int Utilities::edge_swap_host_vecs(
    int iedge, 
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

    // double   generalParams.scaling_pow = 4.0;
    int alpha = 0;
        
    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
    double preferred_angle;
    double vol_0, vol_1;
    double P0x_vol1, P0y_vol1, P0z_vol1, P0x_vol2, P0y_vol2, P0z_vol2;
    double N1x_vol, N1y_vol, N1z_vol, N2x_vol, N2y_vol, N2z_vol;
    // std::cout<<"iedge = "<<iedge<<std::endl;
    //std::cout<<hostSetInfoVecs.edges2Triangles_1[iedge] <<" "<< hostSetInfoVecs.edges2Triangles_2[iedge]<<std::endl;   
    if ( hostSetInfoVecs.edges2Triangles_1[iedge] != hostSetInfoVecs.edges2Triangles_2[iedge]){
        H0 = hostSetInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        //std::cout<<"H0 = "<<H0<<std::endl;
        T0 = hostSetInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        //std::cout<<"T0 = "<<T0<<std::endl;
        edge_start = hostSetInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = hostSetInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = hostSetInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        //std::cout<<"a1 = "<<a1<<std::endl;
        b1 = hostSetInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        //std::cout<<"b1 = "<<b1<<std::endl;
        c1 = hostSetInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        //std::cout<<"c1 = "<<c1<<std::endl;
        
        a2 = hostSetInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        //std::cout<<"a2 = "<<a2<<std::endl;
        b2 = hostSetInfoVecs.triangles2Edges_2[T0];
        //std::cout<<"b2 = "<<b2<<std::endl;
        c2 = hostSetInfoVecs.triangles2Edges_3[T0];
        //std::cout<<"c2 = "<<c2<<std::endl;
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //std::cout<<"H1 = "<<H1<<std::endl;
        //std::cout<<"H2 = "<<H2<<std::endl;
        //std::cout<<"T1 = "<<T1<<std::endl;
        //std::cout<<"T2 = "<<T2<<std::endl;

        //Now search for the associated 

        int CANDIDATE1_1 = hostSetInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = hostSetInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = hostSetInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = hostSetInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = hostSetInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = hostSetInfoVecs.triangles2Nodes_3[T0];
        
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}

        bool BAD_CHOICE = false;
        for (int q = 0; q < 4; q++){
            double qq;
            if (q == 0){
                qq = edge_start;
            }
            else if (q == 1){
                qq = edge_end;
            }
            else if (q == 2){
                qq = HEAD;
            }
            else if (q == 3){
                qq = TAIL;
            }
            int safe_flip1 = 0;
            if (hostSetInfoVecs.nndata1[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata2[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata3[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata4[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata5[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata6[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata7[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata8[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata9[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata10[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata11[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata12[qq] >= 0){safe_flip1 += 1;        }

            if (q == 0 && safe_flip1 == 4){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
               // std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 1 && safe_flip1 == 4){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 2 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 3 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }

            //WE NEED ONE LAST CHECK TO AVOID OVERLAPPING EDGES TO OCCUR, THIS IS EXCLUSIVELY OBSERVED WITH A NARROW BUDDING NECK SHOWS UP IN THE SYSTEM
            if (hostSetInfoVecs.nndata1[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata2[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata3[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata4[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata5[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata6[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata7[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata8[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata9[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            if (hostSetInfoVecs.nndata1[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata2[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata3[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata4[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata5[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata6[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata7[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata8[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata9[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}

        }
        

        if (BAD_CHOICE == false){//(safe_flip1 < generalParams.safeguardthreshold && safe_flip2 < generalParams.safeguardthreshold){



            //int temp_edges2Nodes_2 = HEAD;
            ////std::cout<<"head tail in loop: "<< HEAD << " "<< TAIL <<std::endl;
            //The small subsystem we will be working with is
            //          
            //           edge_start
            //    T10    *   |    *     H10
            //         T1    |     H1
            //        *      |       *
            //    TAIL   T0  |  H0    HEAD
            //        *      |       *
            //         T2    |     H2
            //    T20    *   v    *     H20
            //            edge_end
            //
            //H10 is the triangle sharing the same edge H1 with triangle H0.

            //energy E_0 calculation
            //Since nodes are NOT moved and area is not changed, we only need to calculate 
            //linear spring energy and bending energy.
            //Furthermore, since linear spring energy will only be nontrivial for the edge swapped,
            //we can condense the linear spring energy computation to only one edge.
            //Bending energy is more complicated due to the need of unit normals.
            
            

            std::vector<int> edges_iteration(5);
            edges_iteration[0] = iedge;
            edges_iteration[1] = H1;
            edges_iteration[2] = H2;
            edges_iteration[3] = T1;
            edges_iteration[4] = T2;

            double rep_0;
            double rep_energy_TAIL_HEAD;
            double R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD])*(hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD]) +
                                    (hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD])*(hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD]) +
                                    (hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD])*(hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD]));
            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                if (hostSetInfoVecs.nodes_in_upperhem[TAIL] == 1 && hostSetInfoVecs.nodes_in_upperhem[HEAD] == 1){
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                else{
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                }
            else{
                rep_energy_TAIL_HEAD = 0.0;
                }
            int node_index_H1, node_index_H2, node_index_T1, node_index_T2;
            rep_0 = rep_energy_TAIL_HEAD;
            /*if (generalParams.edges_in_upperhem[edges_iteration[0]] == 1){
                        linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[0]] == 0){
                        linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak + linearSpringInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        linear_spring_constant = linearSpringInfoVecs.spring_constant;
            }*/
            if (generalParams.SCALE_TYPE == 0){
                linear_spring_constant = linearSpringInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],2.0)/generalParams.gausssigma)));
                if (linear_spring_constant < linearSpringInfoVecs.spring_constant_weak){linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;};
            }
            else if (generalParams.SCALE_TYPE == 1){
                // linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.scaling_pow) + 
                //                         linearSpringInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],generalParams.scaling_pow));
                linear_spring_constant = linearSpringInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.scaling_pow) + 
                                        linearSpringInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],generalParams.scaling_pow));
            }
            else if (generalParams.SCALE_TYPE == 2){
                linear_spring_constant = linearSpringInfoVecs.spring_constant - 
                                    (linearSpringInfoVecs.spring_constant - linearSpringInfoVecs.spring_constant_weak)*
                                    hostSetInfoVecs.scaling_per_edge[edges_iteration[0]];
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 1){
                        linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 0){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 0){
                        // linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak*2.0);
                        linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak + linearSpringInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        linear_spring_constant = linearSpringInfoVecs.spring_constant;
                    }
            }
            else if (generalParams.SCALE_TYPE == 4){
                if (generalParams.nonuniform_wall_weakening_linear==true){
                    //double scaling = 0.0;//linearSpringInfoVecs.spring_constant_weak/linearSpringInfoVecs.spring_constant;
                    //linear_spring_constant = linearSpringInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
                    double spectrum = generalParams.maxSpringScaler_linear*linearSpringInfoVecs.spring_constant - linearSpringInfoVecs.spring_constant_weak;
                    linear_spring_constant = linearSpringInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.hilleqnpow)))*spectrum);
                    if (linear_spring_constant < linearSpringInfoVecs.spring_constant_weak){linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;}
                }
                else{
                    if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 1){
                        linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 0){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 0){
                        // linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak*2.0);
                        linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak + linearSpringInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        linear_spring_constant = linearSpringInfoVecs.spring_constant;
                    }
                }
		    }
            
                int wrong1, wrong2, wrong3;
                for (int j = 0; j < 5; j++){
                   /* if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                    }*/
                    if (generalParams.SCALE_TYPE == 0){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],2.0)/generalParams.gausssigma)));
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;};
                       /* if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                        }*/
                    }
                    else if (generalParams.SCALE_TYPE == 1){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*4.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                                            bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                        // bend_spring_constant = bendingTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                        //                     bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                    }   
                    else if (generalParams.SCALE_TYPE == 2){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak - 
                        (bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak)*
                                            hostSetInfoVecs.scaling_per_edge[edges_iteration[j]];
                    }
                    else if (generalParams.SCALE_TYPE == 3){
                        if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1 ){//&& hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                        }
                        else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){// 1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                            //bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                            preferred_angle = (bendingTriangleInfoVecs.initial_angle_bud + bendingTriangleInfoVecs.initial_angle[edges_iteration[j]])/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];
                        }
                    }
                    else if (generalParams.SCALE_TYPE == 4){
                        if (generalParams.nonuniform_wall_weakening_bend==true){
                            //double scaling = 0.0;//bendingTriangleInfoVecs.spring_constant_weak/bendingTriangleInfoVecs.spring_constant;
                            //bend_spring_constant = bendingTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
                            double spectrum = generalParams.maxSpringScaler_bend*bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak;
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*spectrum);
                            if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                            if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1 ){
                                preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                            }
                            else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){
                                preferred_angle = (bendingTriangleInfoVecs.initial_angle_bud + bendingTriangleInfoVecs.initial_angle[edges_iteration[j]])/2.0;
                            }
                            else{
                                preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];
                            }
                        }
                        else{
                            if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1 ){//&& hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                                bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                                preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                            }
                            else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){// 1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                                //bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                                bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                                preferred_angle = (bendingTriangleInfoVecs.initial_angle_bud + bendingTriangleInfoVecs.initial_angle[edges_iteration[j]])/2.0;
                            }
                            else{
                                bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                                preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];
                            }
                        }
                    }
                    int Tri1 = hostSetInfoVecs.edges2Triangles_1[edges_iteration[j]];//index of the 1st triangle
                    int Tri2 = hostSetInfoVecs.edges2Triangles_2[edges_iteration[j]];
                    //int id_k = hostSetInfoVecs.edges2Nodes_1[edges_iteration[j]];
                    //int id_i = hostSetInfoVecs.edges2Nodes_2[edges_iteration[j]];

                    double N1vec1x, N1vec1y, N1vec1z, N1vec2x, N1vec2y, N1vec2z;
                    double N2vec1x, N2vec1y, N2vec1z, N2vec2x, N2vec2y, N2vec2z;
                    double N1_x, N1_y, N1_z, N2_x, N2_y, N2_z;
                    
                    if (Tri1 != Tri2) {
                        int Tri1_n1 = hostSetInfoVecs.triangles2Nodes_1[Tri1];
                        if (j == 0){
                            P0x_vol1 = hostSetInfoVecs.nodeLocX[Tri1_n1];
                            P0y_vol1 = hostSetInfoVecs.nodeLocY[Tri1_n1];
                            P0z_vol1 = hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        }
                        int Tri1_n2 = hostSetInfoVecs.triangles2Nodes_2[Tri1];
                        int Tri1_n3 = hostSetInfoVecs.triangles2Nodes_3[Tri1];
                        N1vec1x = hostSetInfoVecs.nodeLocX[Tri1_n2] - hostSetInfoVecs.nodeLocX[Tri1_n1];
                        N1vec1y = hostSetInfoVecs.nodeLocY[Tri1_n2] - hostSetInfoVecs.nodeLocY[Tri1_n1];
                        N1vec1z = hostSetInfoVecs.nodeLocZ[Tri1_n2] - hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        N1vec2x = hostSetInfoVecs.nodeLocX[Tri1_n3] - hostSetInfoVecs.nodeLocX[Tri1_n1];
                        N1vec2y = hostSetInfoVecs.nodeLocY[Tri1_n3] - hostSetInfoVecs.nodeLocY[Tri1_n1];
                        N1vec2z = hostSetInfoVecs.nodeLocZ[Tri1_n3] - hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        //std::vector<double> N1(3);
                        N1_x = N1vec1y*N1vec2z - N1vec2y*N1vec1z;
                        N1_y = -(N1vec1x*N1vec2z - N1vec2x*N1vec1z);
                        N1_z = N1vec1x*N1vec2y - N1vec2x*N1vec1y;
                        double nN1 = sqrt(pow(N1_x,2)+pow(N1_y,2)+pow(N1_z,2));
                        ////std::cout<<"nN1 = "<<nN1<<std::endl;

                        int Tri2_n1 = hostSetInfoVecs.triangles2Nodes_1[Tri2];
                        int Tri2_n2 = hostSetInfoVecs.triangles2Nodes_2[Tri2];
                        int Tri2_n3 = hostSetInfoVecs.triangles2Nodes_3[Tri2];
                        if (j == 0){
                            P0x_vol2 = hostSetInfoVecs.nodeLocX[Tri2_n1];
                            P0y_vol2 = hostSetInfoVecs.nodeLocY[Tri2_n1];
                            P0z_vol2 = hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        }
                        N2vec1x = hostSetInfoVecs.nodeLocX[Tri2_n2] - hostSetInfoVecs.nodeLocX[Tri2_n1];
                        N2vec1y = hostSetInfoVecs.nodeLocY[Tri2_n2] - hostSetInfoVecs.nodeLocY[Tri2_n1];
                        N2vec1z = hostSetInfoVecs.nodeLocZ[Tri2_n2] - hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        N2vec2x = hostSetInfoVecs.nodeLocX[Tri2_n3] - hostSetInfoVecs.nodeLocX[Tri2_n1];
                        N2vec2y = hostSetInfoVecs.nodeLocY[Tri2_n3] - hostSetInfoVecs.nodeLocY[Tri2_n1];
                        N2vec2z = hostSetInfoVecs.nodeLocZ[Tri2_n3] - hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        //std::vector<double> N2(3);
                        N2_x = N2vec1y*N2vec2z - N2vec2y*N2vec1z;
                        N2_y = -(N2vec1x*N2vec2z - N2vec2x*N2vec1z);
                        N2_z = N2vec1x*N2vec2y - N2vec2x*N2vec1y; 
                        double nN2 = sqrt(pow(N2_x,2)+pow(N2_y,2)+pow(N2_z,2));;
                        ////std::cout<<"nN2 = "<<nN2<<std::endl;

                        double cosAngle = (N1_x*N2_x + N1_y*N2_y + N1_z*N2_z)/(nN1*nN2);
                        ////std::cout<<"dotproduct = "<<cosAngle<<std::endl;
                    
        
                        if (cosAngle > 1.0) {
                            cosAngle = 1.0;
                        }
                        else if (cosAngle < -1.0){
                            cosAngle = -1.0;
                        }

                        double theta_current = acos( cosAngle );
                        
                        
                        double local_energy = bend_spring_constant * (1 - cos(theta_current - preferred_angle) );
                        
                        //bendingTriangleInfoVecs.spring_constant * (1 - cos(theta_current - bendingTriangleInfoVecs.initial_angle) );
                        temp_bend = temp_bend + local_energy;
                        if (j == 1){
                            wrong1 = HEAD;
                            wrong2 = edge_start;
                            wrong3 = edge_end;
                        }
                        if (j == 2){
                            wrong1 = edge_end;
                            wrong2 = HEAD;
                            wrong3 = edge_start;
                        }
                        if (j == 3){
                            wrong1 = edge_start;
                            wrong2 = TAIL;
                            wrong3 = edge_end;
                        }
                        if (j == 4){
                            wrong1 = TAIL;
                            wrong2 = edge_end;
                            wrong3 = edge_start;
                        }
                        if (Tri1_n1 != wrong1 && Tri1_n1 != wrong2 && Tri1_n1 != wrong3){
                            if (j == 1){
                                node_index_H1 = Tri1_n1;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n1;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n1;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n1;
                            }
                        }
                        else if (Tri1_n2 != wrong1 && Tri1_n2 != wrong2 && Tri1_n2 != wrong3){
                            if (j == 1){
                                node_index_H1 = Tri1_n2;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n2;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n2;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n2;
                            }
                        }
                        else if (Tri1_n3 != wrong1 && Tri1_n3 != wrong2 && Tri1_n3 != wrong3){

                            if (j == 1){
                                node_index_H1 = Tri1_n3;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n3;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n3;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n3;
                            }
                        }
                        else if (Tri2_n1 != wrong1 && Tri2_n1 != wrong2 && Tri2_n1 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n1;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n1;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n1;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n1;
                            }
                        }
                        else if (Tri2_n2 != wrong1 && Tri2_n2 != wrong2 && Tri2_n2 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n2;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n2;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n2;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n2;
                            }
                        }
                        else if (Tri2_n3 != wrong1 && Tri2_n3 != wrong2 && Tri2_n3 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n3;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n3;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n3;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n3;
                            }
                        }

                        
                        /*//std::cout<<"bending energy "<<local_energy<<std::endl;
                        for (int COUNT = 0; COUNT < 3; COUNT++){
                        //std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
                        //std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
                        //std::cout<<"angle "<<theta_current<<std::endl;*/
                        if (j == 0){
                            N1x_vol = N1_x/nN1;
                            N1y_vol = N1_y/nN1;
                            N1z_vol = N1_z/nN1;
                            N2x_vol = N2_x/nN2;
                            N2y_vol = N2_y/nN2;
                            N2z_vol = N2_z/nN2;
                        }
                    }
                }

            int nb1, nb2, nb3, nb4, nb5, nb6, nb7, nb8, nb9;
            int midpt, endpt1, endpt2, endpt3, falsept1, falsept2;
            int true_endpoints;
            for (int t = 0; t < 4; t++){    

                if (t == 0){
                    midpt = edge_start;
                    endpt1 = HEAD;
                    endpt2 = TAIL;
                    endpt3 = edge_end;
                    true_endpoints = 3;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_T1;
                }
                else if (t == 1){
                    midpt = HEAD;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = TAIL;
                    true_endpoints = 2;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_H2;
                }
                else if (t == 2){
                    midpt = edge_end;
                    endpt1 = HEAD;
                    endpt2 = edge_start;
                    endpt3 = TAIL;
                    true_endpoints = 3;
                    falsept1 = node_index_H2;
                    falsept2 = node_index_T2;
                }
                else if (t == 3){
                    midpt = TAIL;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = HEAD;
                    true_endpoints = 2;
                    falsept1 = node_index_T1;
                    falsept2 = node_index_T2;
                }
                nb1 = hostSetInfoVecs.nndata1[midpt];
                nb2 = hostSetInfoVecs.nndata2[midpt];
                nb3 = hostSetInfoVecs.nndata3[midpt];
                nb4 = hostSetInfoVecs.nndata4[midpt];
                nb5 = hostSetInfoVecs.nndata5[midpt];
                nb6 = hostSetInfoVecs.nndata6[midpt];
                nb7 = hostSetInfoVecs.nndata7[midpt];
                nb8 = hostSetInfoVecs.nndata8[midpt];
                nb9 = hostSetInfoVecs.nndata9[midpt];
                for (int y = 0; y < 9; y++){
                    int startpt;
                    if (y == 0 && nb1>= 0 && nb1 != endpt1 && nb1 != endpt2 && nb1 != endpt3 && nb1 != falsept1 && nb1 != falsept2){startpt = nb1;}
                    else if (y == 1 && nb2>= 0 && nb2 != endpt1 && nb2 != endpt2 && nb2 != endpt3 && nb2 != falsept1 && nb2 != falsept2){startpt = nb2;}
                    else if (y == 2 && nb3>= 0 && nb3 != endpt1 && nb3 != endpt2 && nb3 != endpt3 && nb3 != falsept1 && nb3 != falsept2){startpt = nb3;}
                    else if (y == 3 && nb4>= 0 && nb4 != endpt1 && nb4 != endpt2 && nb4 != endpt3 && nb4 != falsept1 && nb4 != falsept2){startpt = nb4;}
                    else if (y == 4 && nb5>= 0 && nb5 != endpt1 && nb5 != endpt2 && nb5 != endpt3 && nb5 != falsept1 && nb5 != falsept2){startpt = nb5;}
                    else if (y == 5 && nb6>= 0 && nb6 != endpt1 && nb6 != endpt2 && nb6 != endpt3 && nb6 != falsept1 && nb6 != falsept2){startpt = nb6;}
                    else if (y == 6 && nb7>= 0 && nb7 != endpt1 && nb7 != endpt2 && nb7 != endpt3 && nb7 != falsept1 && nb7 != falsept2){startpt = nb7;}
                    else if (y == 7 && nb8>= 0 && nb8 != endpt1 && nb8 != endpt2 && nb8 != endpt3 && nb8 != falsept1 && nb8 != falsept2){startpt = nb8;}
                    else if (y == 8 && nb9>= 0 && nb9 != endpt1 && nb9 != endpt2 && nb9 != endpt3 && nb9 != falsept1 && nb9 != falsept2){startpt = nb9;}
                    else{continue;}

                    if (true_endpoints == 3){
                        for (int h = 0; h < 3; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            else if (h==2){aa = endpt3;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_0 += rep_energy_TAIL_HEAD;
                        }
                    }
                    else if (true_endpoints == 2){
                        for (int h = 0; h < 2; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                               // rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                               //                         (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                               if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_0 += rep_energy_TAIL_HEAD;
                        }
                    }
                }
                
                //NOW WE UTILIZE THIS LOOP TO CALCULATE THE REPULSION ENERGY ASSOCIATED WITH THE FALSEPTS (I.E. THE ONES ASSOCIATED WITH H10, T10, H20, T20 BUT NOT IN THE SUBSYSTEM DEPECTED ABOVE)
                //int oiu1, oiu2, oiu3;
                // if (t == 0){oiu1 = node_index_H1; oiu2 = edge_end; oiu3 = TAIL;}
                // else if (t == 1){oiu1 = node_index_H2; oiu2 = edge_start; oiu3 = TAIL;}
                // else if (t == 2){oiu1 = node_index_T1; oiu2 = edge_end; oiu3 = HEAD;}
                // else if (t == 3){oiu1 = node_index_T2; oiu2 = edge_start; oiu3 = HEAD;}
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_0 += rep_energy_TAIL_HEAD;

                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_0 += rep_energy_TAIL_HEAD;
            }

            

            double bend_0 = temp_bend;
            //
            double linear_0;
            double DISTANCE = sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0));
            //if (DISTANCE < generalParams.abs_Rmin){
            //    linear_0 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
            //	(DISTANCE - generalParams.Rmin);// + 
                //(linearSpringInfoVecs.spring_constant_rep1/2)*(DISTANCE - generalParams.abs_Rmin)*
                //(DISTANCE - generalParams.abs_Rmin);
                //linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)))*
                //(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)));
            //}
            //else if (DISTANCE != generalParams.Rmin){
            //    linear_0 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
            //	(DISTANCE - generalParams.Rmin);
            //}
            //else{
                // linear_0 = (linear_spring_constant/(2.0*generalParams.Rmin*generalParams.Rmin))*(DISTANCE - generalParams.Rmin)*
                // (DISTANCE - generalParams.Rmin);
                // if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1 && hostSetInfoVecs.boundaries_in_upperhem[edges_iteration[0]] == 0){
                //     linear_0 = (linear_spring_constant/(2.0))*(DISTANCE - generalParams.Rmin_growth)*
                //     (DISTANCE - generalParams.Rmin_growth);
                // }
                // else{
                    linear_0 = (linear_spring_constant/(2.0))*(DISTANCE - generalParams.Rmin)*
                    (DISTANCE - generalParams.Rmin);
                //}
            //}
            
            //else if (DISTANCE < generalParams.Rmin ){
            //    linear_0 = linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
            //}
            
            /*double linear_0 = (linearSpringInfoVecs.spring_constant/2)*(sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0])*
                (sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0]);*/
                ////std::cout<<"the energy of this edge is = "<<linear_0<<std::endl;
            
            int H0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[H0];
            int H0n2 = edge_end;//hostSetInfoVecs.triangles2Nodes_2[H0];
            int H0n3 = HEAD;//hostSetInfoVecs.triangles2Nodes_3[H0];
            int T0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[T0];
            int T0n2 = TAIL;//hostSetInfoVecs.triangles2Nodes_2[T0];
            int T0n3 = edge_end;//hostSetInfoVecs.triangles2Nodes_3[T0];
            double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2.0;
            double area_spring_constant_1, area_spring_constant_2;
            /*if (generalParams.triangles_in_upperhem[H0] == 1){
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
            }
            else if (generalParams.triangles_in_upperhem[H0] == 0){
                area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            }
            else{
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
            }*/
            if (generalParams.SCALE_TYPE == 0){
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if(generalParams.SCALE_TYPE == 1){
                area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
            //  area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
            //                             areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
            //                             areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 2){
             area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H1]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H2]))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                    // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                    area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                }
            }
            else if (generalParams.SCALE_TYPE == 4){
                    if (generalParams.nonuniform_wall_weakening_area==true){
                    // double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
                        // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //                    areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //                    areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
                        double spectrum = generalParams.maxSpringScaler_area*areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                        area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                        areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*spectrum) +
                        areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*spectrum))/3.0;
                        if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
                    }
                    else{
                        if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                            area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                            // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                            area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                        }
                    }
		       }
            /*if (generalParams.triangles_in_upperhem[T0] == 1){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
            }
            else if (generalParams.triangles_in_upperhem[T0] == 0){
                area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            }
            else{
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
            }*/

            if (generalParams.SCALE_TYPE == 0){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if (generalParams.SCALE_TYPE == 1){
                area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)))/3.0;
                // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                //                         areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)) +
                //                         areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 2){
                area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T1])+
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T2]))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                    // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                }
            }
            else if (generalParams.SCALE_TYPE == 4){
                if (generalParams.nonuniform_wall_weakening_area==true){
                    //double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
                    //area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                    //                          areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                    //                          areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
                    double spectrum = generalParams.maxSpringScaler_area*areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                            areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*spectrum) +
                                            areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*spectrum))/3.0;
                    if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else{
                    if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                        area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                        // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                        area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                    }
                }
		   }
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            double area_0_energy = area_spring_constant_1*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                                area_spring_constant_2*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);
            
            //double vol_H0 = (1.0/3.0)*(P0x_vol1*N1x_vol + P0y_vol1*N1y_vol + P0z_vol1*N1z_vol)*area_H0;
            //double vol_T0 = (1.0/3.0)*(P0x_vol2*N2x_vol + P0y_vol2*N2y_vol + P0z_vol2*N2z_vol)*area_T0;
            //vol_0 = vol_H0 + vol_T0;
            double E_0 = linear_0 + bend_0 + area_0_energy + rep_0;// + generalParams.volume_energy;
            // std::cout<<"old linear energy: "<<linear_0<<" , old length = "<<DISTANCE<<std::endl;
            // std::cout<<"old bend energy: "<<bend_0<<std::endl;
            // std::cout<<"old area energy: "<<area_0_energy<<std::endl;
            // std::cout<<"old total energy: "<<E_0<<std::endl;

            
            //Flip the edge, build the data structure for the smaller system.
            /*bool BAD_CHOICE = false;
            int temp_edges2Nodes_1 = TAIL;
            int temp_edges2Nodes_2 = HEAD;

            int temp_nndata_HEAD = nndata[HEAD] + 1;
            int temp_nndata_TAIL = nndata[TAIL] + 1;
            int temp_nndata_edge_start = nndata[edge_start] - 1;
            
            int temp_nndata_edge_end = nndata[edge_end] - 1;
            
            if (boundary_node[HEAD] == false && temp_nndata_HEAD < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[TAIL] == false && temp_nndata_TAIL < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_start] == false && temp_nndata_edge_start < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_end] == false && temp_nndata_edge_end < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[HEAD] == false && temp_nndata_HEAD > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[TAIL] == false && temp_nndata_TAIL > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_start] == false && temp_nndata_edge_start > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_end] == false && temp_nndata_edge_end > 12){
                BAD_CHOICE = true;
            }
            else {
                BAD_CHOICE = false;
            }*/


            if (BAD_CHOICE == false) {
                /*temp_edges2Nodes_1[iedge] = TAIL;
                temp_edges2Nodes_2[iedge] = HEAD;
                temp_nndata[HEAD] = temp_nndata[HEAD] + 1;
                temp_nndata[TAIL] = temp_nndata[TAIL] + 1;
                temp_nndata[edge_start] = temp_nndata[edge_start] - 1;
                temp_nndata[edge_end] = temp_nndata[edge_end] - 1;*/

                //The labeling of neighboring edge will as follows after swap:
                //          
                //           edge_start
                //           *        *
                //         T1    H0     H1
                //        *              *
                //      TAIL ----------> HEAD
                //        *              *
                //         T2    T0    H2
                //           *        *
                //            edge_end
                //
                //Now we will update the temporary data structure to accomodate the edgeswap
                
                //Update the new triangles2Nodes information
                /*temp_triangles2Nodes_1[H0] = HEAD;
                temp_triangles2Nodes_2[H0] = edge_start;
                temp_triangles2Nodes_3[H0] = TAIL;
                temp_triangles2Nodes_1[T0] = HEAD;
                temp_triangles2Nodes_2[T0] = TAIL;
                temp_triangles2Nodes_3[T0] = edge_end;*/

                
                //Creating vectors to compute the normal vectors under the swapped configuration.
                int H1t1 = hostSetInfoVecs.edges2Triangles_1[H1];
                //std::cout<<"H1t1 = "<<H1t1<<std::endl;
                int H1t2 = hostSetInfoVecs.edges2Triangles_2[H1]; 
                //std::cout<<"H1t2 = "<<H1t2<<std::endl;//These are the associated triangles to edge H1
                //For the following if statement, we identify the triangles that are affected by the edge-swap.
                //Since we do not know the exact index of the affected triangle, we use the if statement to consider possible cases.
                //This gives us the vectors necessary to compute unit normal vectors required for bending energy.
                

        
                double H1t1_vec1x;
                double H1t1_vec1y;
                double H1t1_vec1z;
                double H1t1_vec2x;
                double H1t1_vec2y;
                double H1t1_vec2z;
                double H1t2_vec1x;
                double H1t2_vec1y;
                double H1t2_vec1z;
                double H1t2_vec2x;
                double H1t2_vec2y;
                double H1t2_vec2z;

                double H2t1_vec1x;
                double H2t1_vec1y;
                double H2t1_vec1z;
                double H2t1_vec2x;
                double H2t1_vec2y;
                double H2t1_vec2z;
                double H2t2_vec1x;
                double H2t2_vec1y;
                double H2t2_vec1z;
                double H2t2_vec2x;
                double H2t2_vec2y;
                double H2t2_vec2z;

                double T1t2_vec1x;
                double T1t2_vec1y;
                double T1t2_vec1z;
                double T1t2_vec2x;
                double T1t2_vec2y;
                double T1t2_vec2z;
                double T1t1_vec1x;
                double T1t1_vec1y;
                double T1t1_vec1z;
                double T1t1_vec2x;
                double T1t1_vec2y;
                double T1t1_vec2z;

                double T2t2_vec1x;
                double T2t2_vec1y;
                double T2t2_vec1z;
                double T2t2_vec2x;
                double T2t2_vec2y;
                double T2t2_vec2z;
                double T2t1_vec1x;
                double T2t1_vec1y;
                double T2t1_vec1z;
                double T2t1_vec2x;
                double T2t1_vec2y;
                double T2t1_vec2z;
                
                //           edge_start
                //           *        *
                //         T1    H0     H1
                //        *              *
                //      TAIL ----------> HEAD
                //        *              *
                //         T2    T0    H2
                //           *        *
                //            edge_end
                
                if (H1t1 == H0){H1t1_vec1x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t1_vec1y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t1_vec1z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t1_vec2x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t1_vec2y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t1_vec2z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                //std::cout<<"H1t1 = H0"<<std::endl;
                                }
                else if (H1t2 == H0){H1t2_vec1x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t2_vec1y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t2_vec1z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t2_vec2x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t2_vec2y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t2_vec2z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                //std::cout<<"H1t2 = H0"<<std::endl;
                                }
                int H2t1 = hostSetInfoVecs.edges2Triangles_1[H2];
                int H2t2 = hostSetInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
                if (H2t1 == H0){//In this case H2t1 turns into T0.
                                H2t1_vec1x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t1_vec1y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t1_vec1z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t1_vec2x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t1_vec2y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t1_vec2z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                //std::cout<<"H2t1 = H0"<<std::endl;
                                }
                else if (H2t2 == H0){//In this case H2t2 tunrs into T0
                                H2t2_vec1x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t2_vec1y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t2_vec1z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t2_vec2x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t2_vec2y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t2_vec2z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                //std::cout<<"H2t2 = H0"<<std::endl;
                                }
                int T1t1 = hostSetInfoVecs.edges2Triangles_1[T1];
                //std::cout<<"T1t1 = "<<T1t1<<std::endl;
                int T1t2 = hostSetInfoVecs.edges2Triangles_2[T1];
                //std::cout<<"T1t2 = "<<T1t2<<std::endl;
                if (T1t1 == T0){//In this case T1t1 turns into H0.
                                T1t1_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t1_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t1_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t1_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t1_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t1_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                //std::cout<<"T1t1 = T0"<<std::endl;
                                }
                else if (T1t2 == T0){//In this case T1t2 turns into H0.
                                T1t2_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t2_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t2_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t2_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t2_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t2_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                //std::cout<<"T1t2 = T0"<<std::endl;
                                }
                int T2t1 = hostSetInfoVecs.edges2Triangles_1[T2];
                int T2t2 = hostSetInfoVecs.edges2Triangles_2[T2];
                if (T2t1 == T0){T2t1_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t1_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t1_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t1_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t1_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t1_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                //std::cout<<"T2t1 = T0"<<std::endl;
                                }
                else if (T2t2 == T0){
                                T2t2_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t2_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t2_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t2_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t2_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t2_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                //std::cout<<"T2t2 = T0"<<std::endl;
                                }
                
                //First calculate the linear spring energy due to edge-swap.
            double linear_1;
            double DISTANCE = sqrt(
            pow(hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL], 2.0) + 
            pow(hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL], 2.0) + 
            pow(hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL], 2.0));
            


            /*if (DISTANCE < generalParams.abs_Rmin){
                linear_1 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
                (DISTANCE - generalParams.Rmin) + 
                //(linearSpringInfoVecs.spring_constant_rep1/2)*(DISTANCE - generalParams.abs_Rmin)*
                //(DISTANCE - generalParams.abs_Rmin);
                linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)))*
                (1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)));
            }
            else if (DISTANCE != generalParams.Rmin){
                linear_1 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
                (DISTANCE - generalParams.Rmin);
            }*/
            //else{
                // linear_1 = (linearSpringInfoVecs.spring_constant/(2.0*generalParams.Rmin*generalParams.Rmin))*(DISTANCE - generalParams.Rmin)*
                // (DISTANCE - generalParams.Rmin);
                // if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1 && hostSetInfoVecs.boundaries_in_upperhem[edges_iteration[0]] == 0){
                //     linear_1 = (linearSpringInfoVecs.spring_constant/(2.0))*(DISTANCE - generalParams.Rmin_growth)*
                //     (DISTANCE - generalParams.Rmin_growth);
                // }
                //else{
                    linear_1 = (linearSpringInfoVecs.spring_constant/(2.0))*(DISTANCE - generalParams.Rmin)*
                    (DISTANCE - generalParams.Rmin);
                //}
            //}
            
            //else if (DISTANCE < generalParams.Rmin){
            //   linear_1 =   linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
            //}

            double rep_1;
            rep_energy_TAIL_HEAD = 0.0;
            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[edge_end])*(hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[edge_end]) +
                                    (hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[edge_end])*(hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[edge_end]) +
                                    (hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[edge_end])*(hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[edge_end]));
            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                if (hostSetInfoVecs.nodes_in_upperhem[edge_start] == 1 && hostSetInfoVecs.nodes_in_upperhem[edge_end] == 1){
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                else{
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                }
            else{
                rep_energy_TAIL_HEAD = 0.0;
                }
            //double rep_energy_T1, rep_energy_T2, rep_energy_H1, rep_energy_H2;
            rep_1 = rep_energy_TAIL_HEAD;
                
            double prob;
            double random_number;
            double Edif;
            if (DISTANCE >= 0.0){
                //WARNING: RESET BENDING COUNTER
                temp_bend = 0.0;


        
                double N1_vec1x, N1_vec1y, N1_vec1z, N1_vec2x, N1_vec2y, N1_vec2z, N2_vec1x, N2_vec1y, N2_vec1z, N2_vec2x, N2_vec2y, N2_vec2z;
                double N1_x, N1_y, N1_z, N2_x, N2_y, N2_z;
                bool THIS_SHOULD_NOT_HAPPEN = false;
                for (int j = 0; j < 5; j++){
                  /*  if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                    }*/
                    if (generalParams.SCALE_TYPE == 0){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],2.0)/generalParams.gausssigma)));
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;};
                       /*if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                        }*/
                    }
                    else if (generalParams.SCALE_TYPE == 1){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*4.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                                            (bendingTriangleInfoVecs.spring_constant_weak)*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                        // bend_spring_constant = bendingTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                        //                     bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                    }
                    else if (generalParams.SCALE_TYPE == 2){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak - 
                                            (bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak)*
                                            hostSetInfoVecs.scaling_per_edge[edges_iteration[j]];
                    }
                    else if (generalParams.SCALE_TYPE == 3){
                        if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                        }
                        else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){//1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                            // bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                            preferred_angle = (bendingTriangleInfoVecs.initial_angle[edges_iteration[j]] + bendingTriangleInfoVecs.initial_angle_bud)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];
                        }
                    }
                    else if (generalParams.SCALE_TYPE == 4){
                        if (generalParams.nonuniform_wall_weakening_bend==true){
                        //double scaling = 0.0;//bendingTriangleInfoVecs.spring_constant_weak/bendingTriangleInfoVecs.spring_constant;
                            //bend_spring_constant = bendingTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
                            double spectrum = generalParams.maxSpringScaler_bend*bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak;
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*spectrum);
                            if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                        }
                        else{
                            if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                                bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                                preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];//bendingTriangleInfoVecs.initial_angle_bud;
                            }
                            else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){//1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                                // bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                                bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                                preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];//(bendingTriangleInfoVecs.initial_angle + bendingTriangleInfoVecs.initial_angle_bud)/2.0;
                            }
                            else{
                                bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                                preferred_angle = bendingTriangleInfoVecs.initial_angle[edges_iteration[j]];
                            }
                        }
                      }

                        if (j == 0){
                            N1_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];//x component of the 1st vector to calculate N1
                            N1_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                            N1_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N1_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];//x component of the 2nd vector to calculate N1
                            N1_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                            N1_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N2_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                            N2_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                            N2_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N2_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                            N2_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                            N2_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];

                            P0x_vol1 = hostSetInfoVecs.nodeLocX[HEAD];
                            P0y_vol1 = hostSetInfoVecs.nodeLocY[HEAD];
                            P0z_vol1 = hostSetInfoVecs.nodeLocZ[HEAD];
                            P0x_vol2 = hostSetInfoVecs.nodeLocX[HEAD];
                            P0y_vol2 = hostSetInfoVecs.nodeLocY[HEAD];
                            P0z_vol2 = hostSetInfoVecs.nodeLocZ[HEAD];
                        }
                        else if (j == 1){
                            N1_vec1x = H1t1_vec1x;
                            //std::cout<<"N1_vec1x = "<<N1_vec1x<<std::endl;
                            N1_vec1y = H1t1_vec1y;
                            //std::cout<<"N1_vec1y = "<<N1_vec1y<<std::endl;
                            N1_vec1z = H1t1_vec1z;
                            //std::cout<<"N1_vec1z = "<<N1_vec1z<<std::endl;
                            N1_vec2x = H1t1_vec2x;
                            //std::cout<<"N1_vec2x = "<<N1_vec2x<<std::endl;
                            N1_vec2y = H1t1_vec2y;
                            //std::cout<<"N1_vec2y = "<<N1_vec2y<<std::endl;
                            N1_vec2z = H1t1_vec2z;
                            //std::cout<<"N1_vec2z = "<<N1_vec2z<<std::endl;
                            N2_vec1x = H1t2_vec1x;
                            //std::cout<<"N2_vec1x = "<<N2_vec1x<<std::endl;
                            N2_vec1y = H1t2_vec1y;
                            //std::cout<<"N2_vec1y = "<<N2_vec1y<<std::endl;
                            N2_vec1z = H1t2_vec1z;
                            //std::cout<<"N2_vec1z = "<<N2_vec1z<<std::endl;
                            N2_vec2x = H1t2_vec2x;
                            //std::cout<<"N2_vec2x = "<<N2_vec2x<<std::endl;
                            N2_vec2y = H1t2_vec2y;
                            //std::cout<<"N2_vec2y = "<<N2_vec2y<<std::endl;
                            N2_vec2z = H1t2_vec2z;
                            //std::cout<<"N2_vec2z = "<<N2_vec2z<<std::endl;
                        }
                        else if (j == 2){
                            N1_vec1x = H2t1_vec1x;
                            N1_vec1y = H2t1_vec1y;
                            N1_vec1z = H2t1_vec1z;
                            N1_vec2x = H2t1_vec2x;
                            N1_vec2y = H2t1_vec2y;
                            N1_vec2z = H2t1_vec2z;
                            N2_vec1x = H2t2_vec1x;
                            N2_vec1y = H2t2_vec1y;
                            N2_vec1z = H2t2_vec1z;
                            N2_vec2x = H2t2_vec2x;
                            N2_vec2y = H2t2_vec2y;
                            N2_vec2z = H2t2_vec2z;
                        }
                        else if (j == 3){
                            N1_vec1x = T1t1_vec1x;
                            N1_vec1y = T1t1_vec1y;
                            N1_vec1z = T1t1_vec1z;
                            N1_vec2x = T1t1_vec2x;
                            N1_vec2y = T1t1_vec2y;
                            N1_vec2z = T1t1_vec2z;
                            N2_vec1x = T1t2_vec1x;
                            N2_vec1y = T1t2_vec1y;
                            N2_vec1z = T1t2_vec1z;
                            N2_vec2x = T1t2_vec2x;
                            N2_vec2y = T1t2_vec2y;
                            N2_vec2z = T1t2_vec2z;
                        }
                        else if (j == 4){
                            N1_vec1x = T2t1_vec1x;
                            N1_vec1y = T2t1_vec1y;
                            N1_vec1z = T2t1_vec1z;
                            N1_vec2x = T2t1_vec2x;
                            N1_vec2y = T2t1_vec2y;
                            N1_vec2z = T2t1_vec2z;
                            N2_vec1x = T2t2_vec1x;
                            N2_vec1y = T2t2_vec1y;
                            N2_vec1z = T2t2_vec1z;
                            N2_vec2x = T2t2_vec2x;
                            N2_vec2y = T2t2_vec2y;
                            N2_vec2z = T2t2_vec2z;
                        }                        
                        //std::vector<double> N1(3);
                        N1_x = N1_vec1y*N1_vec2z - N1_vec2y*N1_vec1z;
                        N1_y = -(N1_vec1x*N1_vec2z - N1_vec2x*N1_vec1z);
                        N1_z = N1_vec1x*N1_vec2y - N1_vec2x*N1_vec1y;
                        double nN1 = sqrt(pow(N1_x,2)+pow(N1_y,2)+pow(N1_z,2));
                        ////std::cout<<"newnN1 = "<<nN1<<std::endl;

        
                        //std::vector<double> N2(3);
                        N2_x = N2_vec1y*N2_vec2z - N2_vec2y*N2_vec1z;
                        N2_y = -(N2_vec1x*N2_vec2z - N2_vec2x*N2_vec1z);
                        N2_z = N2_vec1x*N2_vec2y - N2_vec2x*N2_vec1y; 
                        double nN2 = sqrt(pow(N2_x,2)+pow(N2_y,2)+pow(N2_z,2));;
                        ////std::cout<<"newnN2 = "<<nN2<<std::endl;

                        if (j == 0){
                            N1x_vol = N1_x/nN1;
                            N1y_vol = N1_y/nN1;
                            N1z_vol = N1_z/nN1;
                            N2x_vol = N2_x/nN2;
                            N2y_vol = N2_y/nN2;
                            N2z_vol = N2_z/nN2;
                        }

                        double cosAngle = (N1_x*N2_x + N1_y*N2_y + N1_z*N2_z)/(nN1*nN2);
                        ////std::cout<<"cosAngle = "<<cosAngle<<std::endl;
                        
                        if (cosAngle > 1.0) {
                            cosAngle = 1.0;
                        }
                        else if (cosAngle < -1.0){
                            cosAngle = -1.0;
                        }
                        if (cosAngle == -1.0){
                            THIS_SHOULD_NOT_HAPPEN = true;
                        }

                        double theta_current = acos( cosAngle );
                        
                        double local_energy = bend_spring_constant * (1 - cos(theta_current - preferred_angle) );
                        temp_bend = temp_bend + local_energy;
                        //std::cout<<"bending energy "<<local_energy<<std::endl;
                        //for (int COUNT = 0; COUNT < 3; COUNT++){
                        ////std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
                        ////std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
                        //std::cout<<"angle "<<theta_current<<std::endl;
                    }
                double bend_1 = temp_bend;

                for (int t = 0; t < 4; t++){    

                if (t == 0){
                    midpt = edge_start;
                    endpt1 = HEAD;
                    endpt2 = TAIL;
                    endpt3 = edge_end;
                    true_endpoints = 2;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_T1;
                }
                else if (t == 1){
                    midpt = HEAD;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = TAIL;
                    true_endpoints = 3;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_H2;
                }
                else if (t == 2){
                    midpt = edge_end;
                    endpt1 = HEAD;
                    endpt2 = edge_start;
                    endpt3 = TAIL;
                    true_endpoints = 2;
                    falsept1 = node_index_H2;
                    falsept2 = node_index_T2;
                }
                else if (t == 3){
                    midpt = TAIL;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = HEAD;
                    true_endpoints = 3;
                    falsept1 = node_index_T1;
                    falsept2 = node_index_T2;
                }
                nb1 = hostSetInfoVecs.nndata1[midpt];
                nb2 = hostSetInfoVecs.nndata2[midpt];
                nb3 = hostSetInfoVecs.nndata3[midpt];
                nb4 = hostSetInfoVecs.nndata4[midpt];
                nb5 = hostSetInfoVecs.nndata5[midpt];
                nb6 = hostSetInfoVecs.nndata6[midpt];
                nb7 = hostSetInfoVecs.nndata7[midpt];
                nb8 = hostSetInfoVecs.nndata8[midpt];
                nb9 = hostSetInfoVecs.nndata9[midpt];
                for (int y = 0; y < 9; y++){
                    int startpt;
                    if (y == 0 && nb1>= 0 && nb1 != endpt1 && nb1 != endpt2 && nb1 != endpt3 && nb1 != falsept1 && nb1 != falsept2){startpt = nb1;}
                    else if (y == 1 && nb2>= 0 && nb2 != endpt1 && nb2 != endpt2 && nb2 != endpt3 && nb2 != falsept1 && nb2 != falsept2){startpt = nb2;}
                    else if (y == 2 && nb3>= 0 && nb3 != endpt1 && nb3 != endpt2 && nb3 != endpt3 && nb3 != falsept1 && nb3 != falsept2){startpt = nb3;}
                    else if (y == 3 && nb4>= 0 && nb4 != endpt1 && nb4 != endpt2 && nb4 != endpt3 && nb4 != falsept1 && nb4 != falsept2){startpt = nb4;}
                    else if (y == 4 && nb5>= 0 && nb5 != endpt1 && nb5 != endpt2 && nb5 != endpt3 && nb5 != falsept1 && nb5 != falsept2){startpt = nb5;}
                    else if (y == 5 && nb6>= 0 && nb6 != endpt1 && nb6 != endpt2 && nb6 != endpt3 && nb6 != falsept1 && nb6 != falsept2){startpt = nb6;}
                    else if (y == 6 && nb7>= 0 && nb7 != endpt1 && nb7 != endpt2 && nb7 != endpt3 && nb7 != falsept1 && nb7 != falsept2){startpt = nb7;}
                    else if (y == 7 && nb8>= 0 && nb8 != endpt1 && nb8 != endpt2 && nb8 != endpt3 && nb8 != falsept1 && nb8 != falsept2){startpt = nb8;}
                    else if (y == 8 && nb9>= 0 && nb9 != endpt1 && nb9 != endpt2 && nb9 != endpt3 && nb9 != falsept1 && nb9 != falsept2){startpt = nb9;}
                    else{continue;}

                    if (true_endpoints == 3){
                        for (int h = 0; h < 3; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            else if (h==2){aa = endpt3;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_1 += rep_energy_TAIL_HEAD;
                        }
                    }
                    else if (true_endpoints == 2){
                        for (int h = 0; h < 2; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_1 += rep_energy_TAIL_HEAD;
                        }
                    }
                }
                
                //NOW WE UTILIZE THIS LOOP TO CALCULATE THE REPULSION ENERGY ASSOCIATED WITH THE FALSEPTS (I.E. THE ONES ASSOCIATED WITH H10, T10, H20, T20 BUT NOT IN THE SUBSYSTEM DEPECTED ABOVE)
                // int oiu1, oiu2, oiu3;
                // if (t == 0){oiu1 = node_index_H1; oiu2 = edge_end; oiu3 = TAIL;}
                // else if (t == 1){oiu1 = node_index_H2; oiu2 = edge_start; oiu3 = TAIL;}
                // else if (t == 2){oiu1 = node_index_T1; oiu2 = edge_end; oiu3 = HEAD;}
                // else if (t == 3){oiu1 = node_index_T2; oiu2 = edge_start; oiu3 = HEAD;}
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_1 += rep_energy_TAIL_HEAD;
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_1 += rep_energy_TAIL_HEAD;
            }

                int H0n1 = HEAD;
                int H0n2 = edge_start;
                int H0n3 = TAIL;
                int T0n1 = TAIL;
                int T0n2 = edge_end;
                int T0n3 = HEAD;
                double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
                double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                            );
                double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                            );
                double mean_abc = (a + b + c)/2;
                double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                            );
                double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                            );
                double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                            );
                double mean_def = (d + e + f)/2.0;
                /*if (generalParams.triangles_in_upperhem[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (generalParams.triangles_in_upperhem[H0] == 0){
                    area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                }*/
                if (generalParams.SCALE_TYPE == 0){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T1],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 1 ){
                    area_spring_constant_1 =  ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                       (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                                       (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)))/3.0;
                    // area_spring_constant_1 =  (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                    //                     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                    //                     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 2){
                    area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H1]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T1]))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 3){
                    if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                        // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                        area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                    }
                }
                else if (generalParams.SCALE_TYPE == 4){
                    if (generalParams.nonuniform_wall_weakening_area==true){
                        //double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
                        //area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //                   areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //                   areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
                        double spectrum = generalParams.maxSpringScaler_area*areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                        area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                                areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*spectrum) +
                                                areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*spectrum))/3.0;
                        if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
                    }
                    else{
                        if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                            area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                            // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                            area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                        }    
                    }
		      }
                /*if (generalParams.triangles_in_upperhem[T0] == 1){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (generalParams.triangles_in_upperhem[T0] == 0){
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                }*/
                if (generalParams.SCALE_TYPE == 0){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T2],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 1){
                    area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)) +
                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
                    // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                    //     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)) +
                    //     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 2){
                    area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T2] )+
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H2]))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 3){
                    if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                        // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                        area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                    }
                }
                else if (generalParams.SCALE_TYPE == 4){
                    if (generalParams.nonuniform_wall_weakening_area==true){
                    //double scaling = 0.0;// areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
                        //area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //               areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                        //               areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
                        double spectrum = generalParams.maxSpringScaler_area*areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                        area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                                areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*spectrum) +
                                                areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*spectrum))/3.0;
                        if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
                    }
                    else{
                         if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                            area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                            // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                            area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                        }
                    }
		   }
                double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
                double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
                double area_1_energy = area_spring_constant_1*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                                    area_spring_constant_2*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);
                //double vol_H0 = (1.0/3.0)*(P0x_vol1*N1x_vol + P0y_vol1*N1y_vol + P0z_vol1*N1z_vol)*area_H0;
                //double vol_T0 = (1.0/3.0)*(P0x_vol2*N2x_vol + P0y_vol2*N2y_vol + P0z_vol2*N2z_vol)*area_T0;
                //vol_1 = vol_H0 + vol_T0;
                //double new_vol = generalParams.true_current_total_volume + (vol_1 - vol_0);
                //double new_vol_energy = generalParams.volume_spring_constant*(new_vol - generalParams.eq_total_volume)*(new_vol - generalParams.eq_total_volume)/
                //                        (2.0*generalParams.Rmin*generalParams.Rmin*generalParams.Rmin*generalParams.eq_total_volume);
                double E_1 = linear_1 + bend_1 + area_1_energy + rep_1;// + new_vol_energy;
                // std::cout<<"new linear energy = "<<linear_1<<" , new length = "<<DISTANCE<<std::endl;
                // std::cout<<"new bend energy = "<<bend_1<<std::endl;
                // std::cout<<"new area energy = "<<area_1_energy<<std::endl;
                // std::cout<<"new total energy: "<<E_1<<std::endl;
            
            //Now compute the Boltzmann factor to determine if a swap occurs.
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, 1.0);
            random_number = dis(gen);
            //double random_number = 0;
            Edif = (E_1 - E_0);
            

            prob = generalParams.tau*exp(-Edif/generalParams.kT);
            
            }
            else{
                prob = -1.0;
            }
            // std::cout<<"P(swap): "<<prob<<std::endl;
            
            bool ACCEPT2;
            if (!isnan(prob)){
                if (prob >= 1){ACCEPT2 = true;}
                else if (prob < 1 && random_number <= prob){ACCEPT2 = true;}
                else if (prob < 1 && random_number > prob){ACCEPT2 = false;}
            }
            else{ACCEPT2 = false;}
            ////std::cout<<"ACCEPT2 = "<<ACCEPT2<<std::endl;
            //Perform real update
            //if (ACCEPT2 == true){
            //if (Edif < 100.0){
            if (ACCEPT2 == true ){//&& THIS_SHOULD_NOT_HAPPEN == false){
                alpha = 1;
                hostSetInfoVecs.triangles2Nodes_1[H0] = HEAD;
                hostSetInfoVecs.triangles2Nodes_2[H0] = edge_start;
                hostSetInfoVecs.triangles2Nodes_3[H0] = TAIL;
                hostSetInfoVecs.triangles2Nodes_1[T0] = HEAD;
                hostSetInfoVecs.triangles2Nodes_2[T0] = TAIL;
                hostSetInfoVecs.triangles2Nodes_3[T0] = edge_end;
                hostSetInfoVecs.triangles2Edges_1[H0] = iedge;
                hostSetInfoVecs.triangles2Edges_2[H0] = H1;
                hostSetInfoVecs.triangles2Edges_3[H0] = T1;
                hostSetInfoVecs.triangles2Edges_1[T0] = iedge;
                hostSetInfoVecs.triangles2Edges_2[T0] = T2;
                hostSetInfoVecs.triangles2Edges_3[T0] = H2;
                hostSetInfoVecs.edges2Nodes_1[iedge] = TAIL;
                hostSetInfoVecs.edges2Nodes_2[iedge] = HEAD;
                H1t1 = hostSetInfoVecs.edges2Triangles_1[H1];
                H1t2 = hostSetInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
                if (H1t1 == H0){hostSetInfoVecs.edges2Triangles_1[H1] = H0;}
                if (H1t2 == H0){hostSetInfoVecs.edges2Triangles_2[H1] = H0;}
                H2t1 = hostSetInfoVecs.edges2Triangles_1[H2];
                H2t2 = hostSetInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
                if (H2t1 == H0){hostSetInfoVecs.edges2Triangles_1[H2] = T0;}
                if (H2t2 == H0){hostSetInfoVecs.edges2Triangles_2[H2] = T0;}
                T1t1 = hostSetInfoVecs.edges2Triangles_1[T1];
                T1t2 = hostSetInfoVecs.edges2Triangles_2[T1];
                if (T1t1 == T0){hostSetInfoVecs.edges2Triangles_1[T1] = H0;}
                if (T1t2 == T0){hostSetInfoVecs.edges2Triangles_2[T1] = H0;}
                T2t1 = hostSetInfoVecs.edges2Triangles_1[T2];
                T2t2 = hostSetInfoVecs.edges2Triangles_2[T2];
                if (T2t1 == T0){hostSetInfoVecs.edges2Triangles_1[T2] = T0;}
                if (T2t2 == T0){hostSetInfoVecs.edges2Triangles_2[T2] = T0;}
                ////std::cout<<"IS THE ERROR HERE 1?"<<std::endl;
                ////////////////////////////////////////////////////////////////////////////////////
                ///////////// UPDATING NEIGHBORING NODE INFO ///////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////

                ///////// DELETING CONNECTIVITY BETWEEN EDGE_START AND EDGE_END ////////////////////
                //int data_id;
                if (hostSetInfoVecs.nndata1[edge_start] == edge_end){
                    //data_id = 0;
                    hostSetInfoVecs.nndata1[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata2[edge_start] == edge_end){
                    //data_id = 1;
                    hostSetInfoVecs.nndata2[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata3[edge_start] == edge_end){
                //    data_id = 2;
                    hostSetInfoVecs.nndata3[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata4[edge_start] == edge_end){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata5[edge_start] == edge_end){
                    //data_id = 4;
                    hostSetInfoVecs.nndata5[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata6[edge_start] == edge_end){
                //    data_id = 5;
                    hostSetInfoVecs.nndata6[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata7[edge_start] == edge_end){
                //    data_id = 6;
                    hostSetInfoVecs.nndata7[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata8[edge_start] == edge_end){
                //   data_id = 7;
                    hostSetInfoVecs.nndata8[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata9[edge_start] == edge_end){
                // data_id = 8;
                    hostSetInfoVecs.nndata9[edge_start] = -2;
                }
                /*else if (hostSetInfoVecs.nndata10[edge_start] == edge_end){
                //   data_id = 9;
                    hostSetInfoVecs.nndata10[edge_start] = -2;
                }
                 else if (hostSetInfoVecs.nndata11[edge_start] == edge_end){
                // data_id = 10;
                    hostSetInfoVecs.nndata11[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata12[edge_start] == edge_end){
                //   data_id = 11;
                    hostSetInfoVecs.nndata12[edge_start] = -2;
                } */
                else {}

                if (hostSetInfoVecs.nndata1[edge_end] == edge_start){
                // data_id = 0;
                    hostSetInfoVecs.nndata1[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata2[edge_end] == edge_start){
                //  data_id = 1;
                    hostSetInfoVecs.nndata2[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata3[edge_end] == edge_start){
                //  data_id = 2;
                    hostSetInfoVecs.nndata3[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata4[edge_end] == edge_start){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata5[edge_end] == edge_start){
                //  data_id = 4;
                    hostSetInfoVecs.nndata5[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata6[edge_end] == edge_start){
                //  data_id = 5;
                    hostSetInfoVecs.nndata6[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata7[edge_end] == edge_start){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata8[edge_end] == edge_start){
                // data_id = 7;
                    hostSetInfoVecs.nndata8[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata9[edge_end] == edge_start){
                // data_id = 8;
                    hostSetInfoVecs.nndata9[edge_end] = -2;
                }
                /*else if (hostSetInfoVecs.nndata10[edge_end] == edge_start){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[edge_end] = -2;
                }
                 else if (hostSetInfoVecs.nndata11[edge_end] == edge_start){
                //   data_id = 10;
                    hostSetInfoVecs.nndata11[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata12[edge_end] == edge_start){
                //   data_id = 11;
                    hostSetInfoVecs.nndata12[edge_end] = -2;
                } */
                else {}
////std::cout<<"IS THE ERROR HERE 2?"<<std::endl;

            
                ///////////// ESTABLISHING NEW CONNECTIVITY ////////////////////
            if (hostSetInfoVecs.nndata1[HEAD] < 0){
                //   data_id = 0;
                    hostSetInfoVecs.nndata1[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata2[HEAD] < 0){
                //   data_id = 1;
                    hostSetInfoVecs.nndata2[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata3[HEAD] < 0){
                //   data_id = 2;
                    hostSetInfoVecs.nndata3[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata4[HEAD] < 0){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata5[HEAD] < 0){
                //   data_id = 4;
                    hostSetInfoVecs.nndata5[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata6[HEAD] < 0){
                //   data_id = 5;
                    hostSetInfoVecs.nndata6[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata7[HEAD] < 0){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata8[HEAD] < 0){
                //  data_id = 7;
                    hostSetInfoVecs.nndata8[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata9[HEAD] < 0){
                //  data_id = 8;
                    hostSetInfoVecs.nndata9[HEAD] = TAIL;
                }
                /*else if (hostSetInfoVecs.nndata10[HEAD] < 0){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[HEAD] = TAIL;
                }
                 else if (hostSetInfoVecs.nndata11[HEAD] < 0){
                //  data_id = 10;
                    hostSetInfoVecs.nndata11[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata12[HEAD] < 0){
                //  data_id = 11;
                    hostSetInfoVecs.nndata12[HEAD] = TAIL;
                } */
                else {}

                ////std::cout<<"IS THE ERROR HERE 3?"<<std::endl;

                if (hostSetInfoVecs.nndata1[TAIL] < 0){
                //  data_id = 0;
                    hostSetInfoVecs.nndata1[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata2[TAIL] < 0){
                //  data_id = 1;
                    hostSetInfoVecs.nndata2[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata3[TAIL] < 0){
                //  data_id = 2;
                    hostSetInfoVecs.nndata3[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata4[TAIL] < 0){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata5[TAIL] < 0){
                //  data_id = 4;
                    hostSetInfoVecs.nndata5[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata6[TAIL] < 0){
                //  data_id = 5;
                    hostSetInfoVecs.nndata6[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata7[TAIL] < 0){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata8[TAIL] < 0){
                //  data_id = 7;
                    hostSetInfoVecs.nndata8[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata9[TAIL] < 0){
                //  data_id = 8;
                    hostSetInfoVecs.nndata9[TAIL] = HEAD;
                }
                /*else if (hostSetInfoVecs.nndata10[TAIL] < 0){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[TAIL] = HEAD;
                }
                 else if (hostSetInfoVecs.nndata11[TAIL] < 0){
                //  data_id = 10;
                    hostSetInfoVecs.nndata11[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata12[TAIL] < 0){
                //  data_id = 11;
                    hostSetInfoVecs.nndata12[TAIL] = HEAD;
                } */
                else {}
                ////std::cout<<"IS THE ERROR HERE 4?"<<std::endl;

                //nndata[HEAD] += 1;
                //nndata[TAIL] += 1;
                //nndata[edge_start] -= 1;
                //nndata[edge_end] -= 1;
                
                for (int j = 0; j < 2; j++){
                    int nn, nnn, nnnn;
                    if (j == 0){
                        nn = edge_start;
                        nnnn = T0;
                    }
                    else if (j == 1){
                        nn = edge_end;
                        nnnn = H0;
                    }

                    if (hostSetInfoVecs.nodes2Triangles_1[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_2[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_3[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_4[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_5[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_6[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_7[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_8[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_9[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
                    }
                }

                for (int j = 0; j < 2; j++){
                    int nn, nnn, nnnn;
                    if (j == 0){
                        nn = HEAD;
                        nnnn = T0;
                    }
                    else if (j == 1){
                        nn = TAIL;
                        nnnn = H0;
                    }

                    if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_1[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_2[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_2[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_3[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_4[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_5[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_6[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_6[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_7[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_7[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_8[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_9[nn] = nnnn;
                    }
                }
        
            }
        } 
    }
    };  
    return alpha;
//This completes the update (if necessary) of the following data structures: triangles2Nodes, edges2Nodes, edges2Triangles.
};



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

/* Make_ele_structure_PDE build the data structure necessary for the PDE computation.
    Not all elements of the structure are introduced as some require another function to 
    compute. 
    Contents not updated here: side_vector, side_norm, aff_conorm, aff_adj_ele_conorm,
                                affine_pt, aff_cent */
/*void Utilities::make_ele_structure_PDE(ELEMENT element,
    int id,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs){
        int i = id;
        element.Id = id;
        int vt1 = coordInfoVecs.triangles2Nodes_1[i];
        int vt2 = coordInfoVecs.triangles2Nodes_2[i];
        int vt3 = coordInfoVecs.triangles2Nodes_3[i];
        int e1 = coordInfoVecs.triangles2Edges_1[i];
        int e2 = coordInfoVecs.triangles2Edges_2[i];
        int e3 = coordInfoVecs.triangles2Edges_3[i];
        std::vector<int> e1v = {coordInfoVecs.edges2Nodes_1[e1], coordInfoVecs.edges2Nodes_2[e1]};
        std::vector<int> e2v = {coordInfoVecs.edges2Nodes_1[e2], coordInfoVecs.edges2Nodes_2[e2]};
        std::vector<int> e3v = {coordInfoVecs.edges2Nodes_1[e3], coordInfoVecs.edges2Nodes_2[e3]};
        element.vertex[0] = vt1;
        element.vertex[1] = vt2;
        element.vertex[2] = vt3;
        element.adj_ele[0] = coordInfoVecs.triangles2Triangles_1[i];
        element.adj_ele[1] = coordInfoVecs.triangles2Triangles_2[i];
        element.adj_ele[2] = coordInfoVecs.triangles2Triangles_3[i];
        element.area0 = areaTriangleInfoVecs.initial_area;
        element.center[0] = (coordInfoVecs.nodeLocX[vt1] + coordInfoVecs.nodeLocX[vt2] + coordInfoVecs.nodeLocX[vt3])/3.0;
        element.center[1] = (coordInfoVecs.nodeLocY[vt1] + coordInfoVecs.nodeLocY[vt2] + coordInfoVecs.nodeLocY[vt3])/3.0;
        element.center[2] = (coordInfoVecs.nodeLocZ[vt1] + coordInfoVecs.nodeLocZ[vt2] + coordInfoVecs.nodeLocZ[vt3])/3.0;
        element.length_side[0] = sqrt((coordInfoVecs.nodeLocX[e1v[1]] - coordInfoVecs.nodeLocX[e1v[0]])*(coordInfoVecs.nodeLocX[e1v[1]] - coordInfoVecs.nodeLocX[e1v[0]]) +
                                    (coordInfoVecs.nodeLocY[e1v[1]] - coordInfoVecs.nodeLocY[e1v[0]])*(coordInfoVecs.nodeLocY[e1v[1]] - coordInfoVecs.nodeLocY[e1v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e1v[1]] - coordInfoVecs.nodeLocZ[e1v[0]])*(coordInfoVecs.nodeLocZ[e1v[1]] - coordInfoVecs.nodeLocZ[e1v[0]]));
        //std::cout<<"length_side[0] = "<<element.length_side[0]<<std::endl;
        element.length_side[1] = sqrt((coordInfoVecs.nodeLocX[e2v[1]] - coordInfoVecs.nodeLocX[e2v[0]])*(coordInfoVecs.nodeLocX[e2v[1]] - coordInfoVecs.nodeLocX[e2v[0]]) +
                                    (coordInfoVecs.nodeLocY[e2v[1]] - coordInfoVecs.nodeLocY[e2v[0]])*(coordInfoVecs.nodeLocY[e2v[1]] - coordInfoVecs.nodeLocY[e2v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e2v[1]] - coordInfoVecs.nodeLocZ[e2v[0]])*(coordInfoVecs.nodeLocZ[e2v[1]] - coordInfoVecs.nodeLocZ[e2v[0]]));
        //std::cout<<"length_side[1] = "<<element.length_side[1]<<std::endl;
        element.length_side[2] = sqrt((coordInfoVecs.nodeLocX[e3v[1]] - coordInfoVecs.nodeLocX[e3v[0]])*(coordInfoVecs.nodeLocX[e3v[1]] - coordInfoVecs.nodeLocX[e3v[0]]) +
                                    (coordInfoVecs.nodeLocY[e3v[1]] - coordInfoVecs.nodeLocY[e3v[0]])*(coordInfoVecs.nodeLocY[e3v[1]] - coordInfoVecs.nodeLocY[e3v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e3v[1]] - coordInfoVecs.nodeLocZ[e3v[0]])*(coordInfoVecs.nodeLocZ[e3v[1]] - coordInfoVecs.nodeLocZ[e3v[0]]));
        //std::cout<<"length_side[2] = "<<element.length_side[2]<<std::endl;
        double dummy_length = (element.length_side[0] + element.length_side[1] + element.length_side[2])/2.0;
        //std::cout<<"dummy_length = "<<dummy_length<<std::endl;
        
        
        element.area = sqrt(dummy_length*(-element.length_side[0] + dummy_length)*(-element.length_side[1] + dummy_length)*(-element.length_side[2] + dummy_length)) ;
        //std::cout<<"element.area ="<<element.area<<std::endl;
        element.out_norm.x[0] = (coordInfoVecs.nodeLocY[vt2] - coordInfoVecs.nodeLocY[vt1])*
                                    (coordInfoVecs.nodeLocZ[vt3] - coordInfoVecs.nodeLocZ[vt1]) - 
                                    (coordInfoVecs.nodeLocY[vt3] - coordInfoVecs.nodeLocY[vt3])*
                                    (coordInfoVecs.nodeLocZ[vt2] - coordInfoVecs.nodeLocZ[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[0]<<std::endl;
        element.out_norm.x[1] = (coordInfoVecs.nodeLocX[vt2] - coordInfoVecs.nodeLocX[vt1])*
                                    (coordInfoVecs.nodeLocZ[vt3] - coordInfoVecs.nodeLocZ[vt1]) - 
                                    (coordInfoVecs.nodeLocX[vt3] - coordInfoVecs.nodeLocX[vt3])*
                                    (coordInfoVecs.nodeLocZ[vt2] - coordInfoVecs.nodeLocZ[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[1]<<std::endl;
        element.out_norm.x[2] = (coordInfoVecs.nodeLocX[vt2] - coordInfoVecs.nodeLocX[vt1])*
                                    (coordInfoVecs.nodeLocY[vt3] - coordInfoVecs.nodeLocY[vt1]) - 
                                    (coordInfoVecs.nodeLocX[vt3] - coordInfoVecs.nodeLocX[vt3])*
                                    (coordInfoVecs.nodeLocY[vt2] - coordInfoVecs.nodeLocY[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[2]<<std::endl;
        element.vertex[0] = coordInfoVecs.triangles2Nodes_1[i];
        element.vertex[1] = coordInfoVecs.triangles2Nodes_2[i];
        element.vertex[2] = coordInfoVecs.triangles2Nodes_3[i];
        element.affine_trans = NULL;
        element.inv_affine_trans = NULL;
    }*/