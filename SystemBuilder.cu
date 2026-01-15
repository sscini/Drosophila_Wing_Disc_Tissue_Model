//#include <curand.h>
#include <set>
#include <list>
#include <vector>
#include <memory>
#include "System.h"
#include "SystemBuilder.h"
#include "SystemStructures.h"
# define M_PI 3.14159265358979323846  /* pi */

// Constructor of SystemBuilder that sets the time step and solve time.
SystemBuilder::SystemBuilder(double timestep_, int solve_time_) {
	dt = timestep_/1000.0;
	solve_time = solve_time_;
};

// Function to add a node to the system with the specified coordinates (x, y, z).
void SystemBuilder::addNode(double x, double y, double z) {
	hostSetInfoVecs.nodeLocX.push_back(x);
	hostSetInfoVecs.nodeLocY.push_back(y);
	hostSetInfoVecs.nodeLocZ.push_back(z);

	hostSetInfoVecs.isNodeFixed.push_back(false);
}

// Function to add non-dimensional node data (x1 - x9) for the node. This is likely related to a specific simulation case.
void SystemBuilder::addNndata(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
	hostSetInfoVecs.nndata1.push_back(x1);
	hostSetInfoVecs.nndata2.push_back(x2);
	hostSetInfoVecs.nndata3.push_back(x3);
	hostSetInfoVecs.nndata4.push_back(x4);
	hostSetInfoVecs.nndata5.push_back(x5);
	hostSetInfoVecs.nndata6.push_back(x6);
	hostSetInfoVecs.nndata7.push_back(x7);
	hostSetInfoVecs.nndata8.push_back(x8);
	hostSetInfoVecs.nndata9.push_back(x9);
}

// Function to add an edge between two nodes with the specified IDs (idL and idR) and calculate its initial length/rest length.
void SystemBuilder::addEdge(int idL, int idR) {
	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
	hostSetInfoVecs.edges2Nodes_2.push_back(idR);

	double xL = hostSetInfoVecs.nodeLocX[idL];
	double yL = hostSetInfoVecs.nodeLocY[idL];
	double zL = hostSetInfoVecs.nodeLocZ[idL];
	double xR = hostSetInfoVecs.nodeLocX[idR];
	double yR = hostSetInfoVecs.nodeLocY[idR];
	double zR = hostSetInfoVecs.nodeLocZ[idR];
	double dist = std::sqrt((xL - xR) * (xL - xR) + (yL - yR) * (yL - yR) + (zL - zR) * (zL - zR));
	hostSetInfoVecs.edge_initial_length.push_back(dist);
  hostSetInfoVecs.edge_rest_length.push_back(dist);
 // std::cout<<"dist = " <<dist<<", from - 1 = "<<idL<< "to - 1  = "<<idR<<std::endl;
}

// Function to add an edge between two nodes with the specified IDs (idL and idR) and specified initial length (edge_initial_length) and rest length (edge_rest_length).
void SystemBuilder::addEdge(int idL, int idR, double restLength) {
	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
	hostSetInfoVecs.edges2Nodes_2.push_back(idR);
	hostSetInfoVecs.edge_initial_length.push_back(restLength);
  hostSetInfoVecs.edge_rest_length.push_back(restLength);
}

void SystemBuilder::addLayerFlag_node(int layerflag){
  hostSetInfoVecs.nodes_in_upperhem.push_back(layerflag);
}

void SystemBuilder::addLayerFlag_edge(int layerflag){
  hostSetInfoVecs.edges_in_upperhem.push_back(layerflag);
}



// Function to add an element (triangle) to the system using three node IDs (idA, idB, and idC).
void SystemBuilder::addElement(int idA, int idB, int idC) {
	hostSetInfoVecs.triangles2Nodes_1.push_back(idA);
	hostSetInfoVecs.triangles2Nodes_2.push_back(idB);
	hostSetInfoVecs.triangles2Nodes_3.push_back(idC);
}

// Function to add an element (triangle) to the system using three edge IDs (idA, idB, and idC).
void SystemBuilder::addElement2Edge(int idA, int idB, int idC) {
	hostSetInfoVecs.triangles2Edges_1.push_back(idA);
	hostSetInfoVecs.triangles2Edges_2.push_back(idB);
	hostSetInfoVecs.triangles2Edges_3.push_back(idC);
}

// Function to add an edge to element mapping by providing element ID (idA) and two edge IDs (idB and idC).
void SystemBuilder::addEdge2Elem(int idA, int idB) {
	hostSetInfoVecs.edges2Triangles_1.push_back(idA);
	hostSetInfoVecs.edges2Triangles_2.push_back(idB);
}

// Function to fix a node with the specified ID (id).
void SystemBuilder::fixNodes(int id) {
	hostSetInfoVecs.isNodeFixed[id] = true;
}

// Function to add a capsid node to the system with the specified coordinates (x, y, z).
void SystemBuilder::addCapsidNode(double x, double y, double z) {
	hostSetInfoVecs.capsidNodeLocX.push_back(x);
	hostSetInfoVecs.capsidNodeLocY.push_back(y);
	hostSetInfoVecs.capsidNodeLocZ.push_back(z);
}

void computeBasisVectors(CoordInfoVecs& coords, HostSetInfoVecs& hostInfo,
                         double theta_DV = 0.0873, double R = 1.0) {
    int N = coords.nodeLocX.size();
    
    hostInfo.pathlength_scaled.resize(N);
    hostInfo.e_R_x.resize(N); hostInfo.e_R_y.resize(N); hostInfo.e_R_z.resize(N);
    hostInfo.e_phi_x.resize(N); hostInfo.e_phi_y.resize(N); hostInfo.e_phi_z.resize(N);
    hostInfo.e_h_x.resize(N); hostInfo.e_h_y.resize(N); hostInfo.e_h_z.resize(N);
    
    for (int i = 0; i < N; i++) {
        double x = coords.nodeLocX[i];
        double y = coords.nodeLocY[i];
        double z = coords.nodeLocZ[i];
        double r = sqrt(x*x + y*y + z*z);
        
        // Surface normal (e_h) - for spherical cap, points radially outward
        hostInfo.e_h_x[i] = x / r;
        hostInfo.e_h_y[i] = y / r;
        hostInfo.e_h_z[i] = z / r;
        
        // DV boundary check
        bool inDV = (fabs(x) <= R * sin(theta_DV / 2.0));
        hostInfo.nodes_in_DV[i] = inDV ? 1 : 0;
        
        // Compute center point for pathlength
        double cx, cy, cz;
        if (inDV) {
            // Center on DV midline
            cx = x;
            cy = 0;
            cz = sqrt(R*R - x*x);
        } else {
            // Center at DV boundary edge
            double sign_x = (x >= 0) ? 1.0 : -1.0;
            double el = sign_x * theta_DV / 2.0;
            cx = R * sin(el);
            cy = 0;
            cz = R * cos(el);
        }
        
        // Pathlength: arc distance from center
        double dot = (x*cx + y*cy + z*cz) / (R*R);
        dot = fmax(-1.0, fmin(1.0, dot));  // Clamp for numerical safety
        double pathlength = R * acos(dot);
        hostInfo.pathlength_scaled[i] = pathlength;  // Will normalize later
        
        // e_R: in-surface direction pointing away from center
        double oax = x - cx, oay = y - cy, oaz = z - cz;
        double oa_len = sqrt(oax*oax + oay*oay + oaz*oaz);
        if (oa_len < 1e-10) oa_len = 1e-10;
        oax /= oa_len; oay /= oa_len; oaz /= oa_len;
        
        // Project onto surface: e_R = e_OA - (e_h · e_OA) * e_h
        double eh_dot_eoa = hostInfo.e_h_x[i]*oax + hostInfo.e_h_y[i]*oay + hostInfo.e_h_z[i]*oaz;
        double er_x = oax - eh_dot_eoa * hostInfo.e_h_x[i];
        double er_y = oay - eh_dot_eoa * hostInfo.e_h_y[i];
        double er_z = oaz - eh_dot_eoa * hostInfo.e_h_z[i];
        double er_len = sqrt(er_x*er_x + er_y*er_y + er_z*er_z);
        if (er_len < 1e-10) er_len = 1e-10;
        hostInfo.e_R_x[i] = er_x / er_len;
        hostInfo.e_R_y[i] = er_y / er_len;
        hostInfo.e_R_z[i] = er_z / er_len;
        
        // e_phi = e_h × e_R
        hostInfo.e_phi_x[i] = hostInfo.e_h_y[i] * hostInfo.e_R_z[i] - hostInfo.e_h_z[i] * hostInfo.e_R_y[i];
        hostInfo.e_phi_y[i] = hostInfo.e_h_z[i] * hostInfo.e_R_x[i] - hostInfo.e_h_x[i] * hostInfo.e_R_z[i];
        hostInfo.e_phi_z[i] = hostInfo.e_h_x[i] * hostInfo.e_R_y[i] - hostInfo.e_h_y[i] * hostInfo.e_R_x[i];
    }
    
    // Normalize pathlength_scaled separately for DV and outDV
    double max_pl_inDV = 0, max_pl_outDV = 0;
    for (int i = 0; i < N; i++) {
        if (hostInfo.nodes_in_DV[i] == 1) {
            max_pl_inDV = fmax(max_pl_inDV, hostInfo.pathlength_scaled[i]);
        } else {
            max_pl_outDV = fmax(max_pl_outDV, hostInfo.pathlength_scaled[i]);
        }
    }
    for (int i = 0; i < N; i++) {
        if (hostInfo.nodes_in_DV[i] == 1 && max_pl_inDV > 0) {
            hostInfo.pathlength_scaled[i] /= max_pl_inDV;
        } else if (max_pl_outDV > 0) {
            hostInfo.pathlength_scaled[i] /= max_pl_outDV;
        }
    }
}

// Function to create the system and return a shared pointer to it.
std::shared_ptr<System> SystemBuilder::createSystem() {
	// Create a shared pointer to the System object.
	std::shared_ptr<System> host_ptr_System = std::make_shared<System>();

	// Set individual parameters in the System object using the values stored in hostSetInfoVecs.
	host_ptr_System->generalParams.dt = defaultdt;
	host_ptr_System->generalParams.solve_time = solve_time;
	host_ptr_System->generalParams.tau = defaultTau;
	host_ptr_System->generalParams.kT = defaultKBT;
  host_ptr_System->generalParams.Tf = defaultTf;
  host_ptr_System->generalParams.tol = defaulttol;
  //host_ptr_System->generalParams.Rmin = defaultRmin;
  host_ptr_System->generalParams.lambda_iso_center_outDV = defaultlambda_iso_out_DV_center;
  host_ptr_System->generalParams.lambda_iso_edge_outDV = defaultlambda_iso_out_DV_edge;
  host_ptr_System->generalParams.lambda_aniso_center_outDV = defaultlambda_aniso_out_DV_center;
  host_ptr_System->generalParams.lambda_aniso_edge_outDV = defaultlambda_aniso_out_DV_edge;
  host_ptr_System->generalParams.lambda_iso_center_inDV = defaultlambda_iso_in_DV_center;
  host_ptr_System->generalParams.lambda_iso_edge_inDV = defaultlambda_iso_in_DV_edge;
  host_ptr_System->generalParams.lambda_aniso_center_inDV = defaultlambda_aniso_in_DV_center;
  host_ptr_System->generalParams.lambda_aniso_edge_inDV = defaultlambda_aniso_in_DV_edge;
  
	host_ptr_System->linearSpringInfoVecs.spring_constant = defaultLinear_Const;
	host_ptr_System->areaTriangleInfoVecs.spring_constant = defaultArea_Const;
	host_ptr_System->bendingTriangleInfoVecs.spring_constant = defaultBending_Const;

	host_ptr_System->ljInfoVecs.epsilon_M = defaultLJ_Eps;
	host_ptr_System->ljInfoVecs.Rmin_M = defaultLJ_Rmin;
	host_ptr_System->ljInfoVecs.Rcutoff_M = defaultLJ_Rmax;
	host_ptr_System->ljInfoVecs.spring_constant = defaultLJ_Const;
	host_ptr_System->ljInfoVecs.LJ_PosX = defaultLJ_X;
	host_ptr_System->ljInfoVecs.LJ_PosY = defaultLJ_Y;
	host_ptr_System->ljInfoVecs.LJ_PosZ = defaultLJ_Z;

	host_ptr_System->linearSpringInfoVecs.scalar_edge_length = defaultEdgeEq;
	host_ptr_System->areaTriangleInfoVecs.initial_area = defaultAreaEq;
  host_ptr_System->generalParams.num_layers = defaultLayers;
	//host_ptr_System->bendingTriangleInfoVecs.initial_angle = defaultAngleEq;

	// Initialize the System object with the data stored in hostSetInfoVecs.
	host_ptr_System->initializeSystem(hostSetInfoVecs);

	// Return the shared pointer to the System object.
	return host_ptr_System;
}


////#include <curand.h>
//#include <set>
//#include <list>					 
//#include <vector>
//#include <memory>
//#include "System.h"
//#include "SystemBuilder.h"
//#include "SystemStructures.h"
//# define M_PI 3.14159265358979323846  /* pi */
//
//// This is where the time step is affected. 
//SystemBuilder::SystemBuilder(double timestep_, int solve_time_){
//	dt = timestep_;
//	solve_time = solve_time_;
//};
//
//void SystemBuilder::addNode(double x,double y, double z ) {
//	hostSetInfoVecs.nodeLocX.push_back(x);
//	hostSetInfoVecs.nodeLocY.push_back(y);
//	hostSetInfoVecs.nodeLocZ.push_back(z);
//
//	hostSetInfoVecs.isNodeFixed.push_back(false);
//
//
//	//std::cout<<"adding node: "<< x << " " << y << " "<< z <<  std::endl;
//}
//
//void SystemBuilder::addNndata(double x1,double x2, double x3, double x4,double x5, double x6, double x7,double x8, double x9) {
////void SystemBuilder::addNndata(double x1,double x2, double x3, double x4,double x5, double x6, double x7,double x8, double x9, double x10,double x11, double x12 ) {
//	hostSetInfoVecs.nndata1.push_back(x1);
//	hostSetInfoVecs.nndata2.push_back(x2);
//	hostSetInfoVecs.nndata3.push_back(x3);
//	hostSetInfoVecs.nndata4.push_back(x4);
//	hostSetInfoVecs.nndata5.push_back(x5);
//	hostSetInfoVecs.nndata6.push_back(x6);
//	hostSetInfoVecs.nndata7.push_back(x7);
//	hostSetInfoVecs.nndata8.push_back(x8);
//	hostSetInfoVecs.nndata9.push_back(x9);
//	//hostSetInfoVecs.nndata10.push_back(x10);
//	//hostSetInfoVecs.nndata11.push_back(x11);
//	//hostSetInfoVecs.nndata12.push_back(x12);
//
//	//std::cout<<"adding node: "<< x << " " << y << " "<< z <<  std::endl;
//}
//
//void SystemBuilder::addEdge(int idL, int idR ) {
//	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
//	hostSetInfoVecs.edges2Nodes_2.push_back(idR);
//
//
//	double xL = hostSetInfoVecs.nodeLocX[idL];
//	double yL = hostSetInfoVecs.nodeLocY[idL];
//	double zL = hostSetInfoVecs.nodeLocZ[idL];
//	double xR = hostSetInfoVecs.nodeLocX[idR];
//	double yR = hostSetInfoVecs.nodeLocY[idR];
//	double zR = hostSetInfoVecs.nodeLocZ[idR];
//	double dist = std::sqrt( (xL-xR)*(xL-xR) + (yL-yR)*(yL-yR) + (zL-zR)*(zL-zR));
//	hostSetInfoVecs.edge_initial_length.push_back(dist);
//	//std::cout<<"adding edge with calc dist: "<< idL << " " << idR << " dist: "<< dist<< std::endl;
//}
//void SystemBuilder::addEdge(int idL, int idR, double edge_initial_length) {
//	
//	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
//	hostSetInfoVecs.edges2Nodes_2.push_back(idR);
//
//	hostSetInfoVecs.edge_initial_length.push_back(edge_initial_length);
//	
//	//std::cout<<"adding edge with fixed dist: "<< idL << " " << idR << " dist: "<< edge_initial_length<< std::endl;
//}
//void SystemBuilder::addElement(int idA, int idB, int idC ) {
//	hostSetInfoVecs.triangles2Nodes_1.push_back(idA);
//	hostSetInfoVecs.triangles2Nodes_2.push_back(idB);	
//	hostSetInfoVecs.triangles2Nodes_3.push_back(idC);
//	
//}
//
//void SystemBuilder::addElement2Edge(int idA, int idB, int idC ) {
//	hostSetInfoVecs.triangles2Edges_1.push_back(idA);
//	hostSetInfoVecs.triangles2Edges_2.push_back(idB);	
//	hostSetInfoVecs.triangles2Edges_3.push_back(idC);
//	
//}
//
//void SystemBuilder::addEdge2Elem(int idA, int idB ) {
//	hostSetInfoVecs.edges2Triangles_1.push_back(idA);
//	hostSetInfoVecs.edges2Triangles_2.push_back(idB);	
//}
//
//void SystemBuilder::fixNodes(int id){
//	hostSetInfoVecs.isNodeFixed[id] = true;
//	//std::cout<<"fixing node "<< id << std::endl;
//}
//void SystemBuilder::addCapsidNode(double x, double y, double z){
//	hostSetInfoVecs.capsidNodeLocX.push_back(x);
//	hostSetInfoVecs.capsidNodeLocY.push_back(y);
//	hostSetInfoVecs.capsidNodeLocZ.push_back(z);
//	//std::cout<<"adding capsid node: "<< x << " " << y << " "<< z <<  std::endl;
//}
//
//
////adds all constraints to the nodesystem model so that it can use the constraints.
//std::shared_ptr<System> SystemBuilder::createSystem() {
//
//
//	//now all the edges and variables are set. 
//	//so set the system and return a pointer.
//	std::shared_ptr<System> host_ptr_System = std::make_shared<System>();
//	//hand in individually set parameters here:
//
//
//	//set individual parameters
//	host_ptr_System->generalParams.dt = dt;
//	host_ptr_System->generalParams.solve_time = solve_time;//int
//	host_ptr_System->generalParams.tau = defaultTau; 
//	host_ptr_System->generalParams.kT = defaultKBT;
//
//	host_ptr_System->linearSpringInfoVecs.spring_constant = defaultLinear_Const;
//	host_ptr_System->areaTriangleInfoVecs.spring_constant = defaultArea_Const;
//	host_ptr_System->bendingTriangleInfoVecs.spring_constant = defaultBending_Const;
//
//	host_ptr_System->ljInfoVecs.epsilon_M = defaultLJ_Eps;
//	host_ptr_System->ljInfoVecs.Rmin_M = defaultLJ_Rmin;
//	host_ptr_System->ljInfoVecs.Rcutoff_M = defaultLJ_Rmax;
//	host_ptr_System->ljInfoVecs.spring_constant = defaultLJ_Const;
//	host_ptr_System->ljInfoVecs.LJ_PosX = defaultLJ_X;
//	host_ptr_System->ljInfoVecs.LJ_PosY = defaultLJ_Y;
//	host_ptr_System->ljInfoVecs.LJ_PosZ = defaultLJ_Z;
//
//
//	host_ptr_System->linearSpringInfoVecs.scalar_edge_length = defaultEdgeEq;
//	host_ptr_System->areaTriangleInfoVecs.initial_area = defaultAreaEq;
//	host_ptr_System->bendingTriangleInfoVecs.initial_angle = defaultAngleEq;
//
//
//	//set vectors and allocate memory here:
//	host_ptr_System->initializeSystem(
//		hostSetInfoVecs
//	);
//	
//	
//	return host_ptr_System;							 
//
//}
//
