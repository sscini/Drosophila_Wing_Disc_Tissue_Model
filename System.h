

#ifndef SYSTEM_H_
#define SYSTEM_H_

#pragma once

#include <fstream>
#include <memory>
#include <math.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <stdint.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

using CVec3   = thrust::tuple<double,double,double>;    // 3-component vector
using Mat_3x3 = thrust::tuple<CVec3,CVec3,CVec3>;       // 3×3 matrix (row-major)

// Struct to store information related to the capsid (virus shell) nodes.
struct CapsidInfoVecs {
	  // Bucket IDs for each capsid node.
    thrust::device_vector<int> id_bucket;	//bucket id
	  // bucket value means global rank of a certain point
	  // Node IDs (global rank) for each capsid node.
    thrust::device_vector<int> id_value;//node id

    // X, Y, and Z coordinates of the capsid nodes.
    thrust::device_vector<double> nodeLocX;
    thrust::device_vector<double> nodeLocY;
    thrust::device_vector<double> nodeLocZ;
  
  
  	// Temporary vectors used for storing information related to linear springs binding capsid nodes to the membrane.
    thrust::device_vector<int> tempMembraneId;
  	thrust::device_vector<int> tempCapsideId;
  
  	thrust::device_vector<double> tempLengthsPairs;
  	thrust::device_vector<double> tempNodeForceX;
  	thrust::device_vector<double> tempNodeForceY;
  	thrust::device_vector<double> tempNodeForceZ;
  
  	// Parameters related to capsid repulsion.
    int factor = 10; // Number of default nodes repulsion can interact with.
    thrust::device_vector<int> tempNodeIdUnreduced;
  	thrust::device_vector<double> tempNodeForceXUnreduced;
  	thrust::device_vector<double> tempNodeForceYUnreduced;
  	thrust::device_vector<double> tempNodeForceZUnreduced;
  
  	// Reduced vectors to store information related to capsid repulsion.
    thrust::device_vector<int> tempNodeIdReduced;
  	thrust::device_vector<double> tempNodeForceXReduced;
  	thrust::device_vector<double> tempNodeForceYReduced;
  	thrust::device_vector<double> tempNodeForceZReduced;
   
    // Parameters related to VerletVelocity node update algorithm
    thrust::device_vector<double> nodeVelX;
    thrust::device_vector<double> nodeVelY;
    thrust::device_vector<double> nodeVelZ;
  
  	// Parameters for the linear spring connecting capsid nodes to the membrane.
    double spring_constant = 5.0;
  	double length_zero = 0.97;
  	double length_cutoff = 1.63;
  	int maxNodeCount;
  	double viscosity = 1.0;
  	double forceX = 0.0;
  	double forceY = 0.0;
  	double forceZ = 0.0;
    
  	int num_connections=0;// Number of connections between capsid nodes.

  
};

struct PrismInfoVecs {
    // for each basic prism element P we store its 6 node indices. 
    thrust::device_vector<int> P1; // P1 -> The center point for all node
    thrust::device_vector<int> P2; // P2
    thrust::device_vector<int> P3; // P3
    thrust::device_vector<int> P4; // P4
    thrust::device_vector<int> P5; // P5
    thrust::device_vector<int> P6; // P6
    
    int num_prisms;
};

// Struct to store node location, velocity, and force information.
struct CoordInfoVecs {
    // Parameters for Foppl-von Karman model.
  	double k_0;
  	double k_1;
  	double k_2;
  	double k_3;
  	double k_4;
  	double k_ss;//10.75; // Here is the ks for Foppl-von Karman number k. Now to find kb.
    double beta;///1.45; // Coefficient for bending energy.
    double gamma; // Coefficient for triangular energy.
    double q1; // Coefficient for triangular energy.
    double h; // Coefficient for triangular energy.
  	thrust::device_vector<double> b_per_triangle; // Unused vector (not annotated).
  
  	// Displacement vector u for Foppl-von Karman model.
    thrust::device_vector<double> u;
  
  	// Solution vector per triangle for Foppl-von Karman model.
    thrust::device_vector<double> soln_per_triangle;
  
  	// Scaling vector per edge.
    thrust::device_vector<double> scaling_per_edge;
    
    // Vector for storing old and new positions of all the nodes. 
    //thrust::device_vector<CVec3> pos_old;
    //thrust::device_vector<CVec3> pos_new; 
    
    
  	// Data structure to store node-to-triangle connectivity (up to 24 triangles per node). 
    thrust::device_vector<int> nodes2Triangles_1;
  	thrust::device_vector<int> nodes2Triangles_2;
  	thrust::device_vector<int> nodes2Triangles_3;
  	thrust::device_vector<int> nodes2Triangles_4;
  	thrust::device_vector<int> nodes2Triangles_5;
  	thrust::device_vector<int> nodes2Triangles_6;
  	thrust::device_vector<int> nodes2Triangles_7;
  	thrust::device_vector<int> nodes2Triangles_8;
  	thrust::device_vector<int> nodes2Triangles_9;
  	thrust::device_vector<int> nodes2Triangles_10;
  	thrust::device_vector<int> nodes2Triangles_11;
  	thrust::device_vector<int> nodes2Triangles_12;
    thrust::device_vector<int> nodes2Triangles_13;
    thrust::device_vector<int> nodes2Triangles_14;
    thrust::device_vector<int> nodes2Triangles_15;
    thrust::device_vector<int> nodes2Triangles_16;
    thrust::device_vector<int> nodes2Triangles_17;
    thrust::device_vector<int> nodes2Triangles_18;
    thrust::device_vector<int> nodes2Triangles_19;
    thrust::device_vector<int> nodes2Triangles_20;
    thrust::device_vector<int> nodes2Triangles_21;
    thrust::device_vector<int> nodes2Triangles_22;
    thrust::device_vector<int> nodes2Triangles_23;
    thrust::device_vector<int> nodes2Triangles_24;
  
  	// Surface normal vectors for each triangle.
    thrust::device_vector<double> SurfaceNormalX;
  	thrust::device_vector<double> SurfaceNormalY;
  	thrust::device_vector<double> SurfaceNormalZ;
  
  	// Non-dimensional node data (up to 20 values per node).
    thrust::device_vector<int> nndata1;
  	thrust::device_vector<int> nndata2;
  	thrust::device_vector<int> nndata3;
  	thrust::device_vector<int> nndata4;
  	thrust::device_vector<int> nndata5;
  	thrust::device_vector<int> nndata6;
  	thrust::device_vector<int> nndata7;
  	thrust::device_vector<int> nndata8;
  	thrust::device_vector<int> nndata9;
  	thrust::device_vector<int> nndata10;
  	thrust::device_vector<int> nndata11;
  	thrust::device_vector<int> nndata12;
  	thrust::device_vector<int> nndata13;
  	thrust::device_vector<int> nndata14;
  	thrust::device_vector<int> nndata15;
  	thrust::device_vector<int> nndata16;
  	thrust::device_vector<int> nndata17;
  	thrust::device_vector<int> nndata18;
  	thrust::device_vector<int> nndata19;
  	thrust::device_vector<int> nndata20;
  
  	// Vector to store whether a node is fixed (true) or not (false).
    thrust::device_vector<bool> isNodeFixed;
  	//GLOBAL COORDS
  	// Previous node location, velocity, and force.
    // X,Y,Z, location, velocity and force of all nodes
  	thrust::device_vector<double> prevNodeLocX;
  	thrust::device_vector<double> prevNodeLocY;
  	thrust::device_vector<double> prevNodeLocZ;
  	thrust::device_vector<double> prevNodeForceX;
  	thrust::device_vector<double> prevNodeForceY;
  	thrust::device_vector<double> prevNodeForceZ;
  
   	// Current node location, velocity, and force.
    thrust::device_vector<double> nodeLocX;
  	thrust::device_vector<double> nodeLocY;
  	thrust::device_vector<double> nodeLocZ;
  	thrust::device_vector<double> nodeForceX;
  	thrust::device_vector<double> nodeForceY;
  	thrust::device_vector<double> nodeForceZ;
  
  	// Temporary vectors for reduced and unreduced node forces during computation.
    thrust::device_vector<double> tempNodeForceXReduced;
  	thrust::device_vector<double> tempNodeForceYReduced;
  	thrust::device_vector<double> tempNodeForceZReduced;
  	thrust::device_vector<double> tempNodeForceXUnreduced;
  	thrust::device_vector<double> tempNodeForceYUnreduced;
  	thrust::device_vector<double> tempNodeForceZUnreduced;
  
  	//LOCAL COORDS
  	// Data structures to store triangle and edge information.
    //indices of each triangle
  	int num_triangles;
  	thrust::device_vector<int> triangles2Nodes_1;
  	thrust::device_vector<int> triangles2Nodes_2;
  	thrust::device_vector<int> triangles2Nodes_3;
  
  
  	//indices of each edge
  	int num_edges;
  	thrust::device_vector<int> edges2Nodes_1;
  	thrust::device_vector<int> edges2Nodes_2;
  
  	//indices of 2 triangle on each edge
  	thrust::device_vector<int> edges2Triangles_1;
  	thrust::device_vector<int> edges2Triangles_2;
  
  	//indices of edges on each triangle.
  	thrust::device_vector<int> triangles2Edges_1;
  	thrust::device_vector<int> triangles2Edges_2;
  	thrust::device_vector<int> triangles2Edges_3;
  
  	thrust::device_vector<int> triangles2Triangles_1;
  	thrust::device_vector<int> triangles2Triangles_2;
  	thrust::device_vector<int> triangles2Triangles_3;
    
    thrust::device_vector<double> pathlength_scaled;
    
    // Radial basis vector (in-surface direction pointing away from center)
    thrust::device_vector<double> e_R_x;
    thrust::device_vector<double> e_R_y;
    thrust::device_vector<double> e_R_z;
    
    // Circumferential basis vector (in-surface, perpendicular to e_R)
    thrust::device_vector<double> e_phi_x;
    thrust::device_vector<double> e_phi_y;
    thrust::device_vector<double> e_phi_z;
    
    // Height/normal basis vector (surface normal, points radially outward)
    thrust::device_vector<double> e_h_x;
    thrust::device_vector<double> e_h_y;
    thrust::device_vector<double> e_h_z;

};


// Struct used for linking of nodes in the network.
struct AuxVecs {

// Bucket IDs for each node, i.e. which point fits into which bucket ID. These are bucket keys. 
    thrust::device_vector<int> id_bucket; // bucket id
    // Node IDs (global rank) for each node. <- Bucket Value
    thrust::device_vector<int> id_value; // node id
    // Bucket keys expanded means the bucket IDs of neighbors for each node.
    thrust::device_vector<int> id_bucket_expanded;
    // Bucket values expanded means each node (represented by its global rank) will have multiple copies.
    thrust::device_vector<int> id_value_expanded;

    // Begin position of keys in id_bucket_expanded and id_value_expanded.
    // Entry keyBegin[bucketKey] returns the start of indices to link.
    thrust::device_vector<int> keyBegin;
    // End position of keys in id_bucket_expanded and id_value_expanded.
    thrust::device_vector<int> keyEnd;
};


// Struct to store domain parameters.
struct DomainParams {
    // Domain boundaries.
    double minX;
    double maxX;
    double minY;
    double maxY;
    double minZ;
    double maxZ;
    // Domain boundaries in the original configuration (used for visualization).
    double originMinX;
    double originMaxX;
    double originMinY;
    double originMaxY;
    double originMinZ;
    double originMaxZ;
    double gridSpacing = 1.5;// Bucket scheme search distance (must be larger than the cutoff for capsid).
    int XBucketCount;
  	int YBucketCount;
  	int ZBucketCount;
  	int totalBucketCount = 0;
};

// Struct to store Lennard-Jones (LJ) potential information.
struct LJInfoVecs{
	  // Position of the Lennard-Jones particle in the X, Y, and Z coordinates.
    double LJ_PosX;
  	double LJ_PosY;
  	double LJ_PosZ;
	  // Vectors to store the positions of all Lennard-Jones particles in the X, Y, and Z coordinates.
    thrust::device_vector<double> LJ_PosX_all;
  	thrust::device_vector<double> LJ_PosY_all;
  	thrust::device_vector<double> LJ_PosZ_all;



  	// Parameters for LJ potential.
    double Rmin_M=0.97;//1.0; // Lennard-Jones minimum energy distance.
    double Rcutoff_M=0.97;//1.0; // Lennard-Jones cutoff distance.
    double Rmin_LJ;
    double Rcutoff_LJ;

  	double epsilon_M=1.0; // Lennard-Jones well depth.
  	double epsilon_M_att1;
  	double epsilon_M_att2;
  	double epsilon_M_rep1;
  	double epsilon_M_rep2;
  	double epsilon_LJ;
  	double epsilon_LJ_rep1;
  	double epsilon_LJ_rep2;
  	double spring_constant;
  
  	// Temporary vectors used for computation.
    thrust::device_vector<int> node_id_close;
  	double lj_energy_M;
  	double lj_energy_LJ;
  	double forceX;
  	double forceY;
  	double forceZ;
  	thrust::device_vector<double> forceX_all;
  	thrust::device_vector<double> forceY_all;
  	thrust::device_vector<double> forceZ_all;

};

// Struct to store information related to area springs (triangular elements).
struct AreaTriangleInfoVecs {
	double dummy;
	int factor = 3; // Used for reduction.
	double initial_area = 0.0048013; // Initial area of triangle. 
	double spring_constant = 5.0;
	double spring_constant_weak = 5.0;

	double area_triangle_energy;
  
  // Temporary vectors used for computation. 
	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;

};

// Struct to store information related to bending springs (triangular elements).
struct BendingTriangleInfoVecs {
	int numBendingSprings=0;

	int factor = 4; // Used for reduction.
	double spring_constant = 5.0;
	double spring_constant_weak;
	double spring_constant_raft;
	double spring_constant_coat;
	//double initial_angle = 0.0087; // Initial angle of the bending triangle (radians).
	double initial_angle_raft;
	double initial_angle_coat;
	double initial_angle_bud;

	double bending_triangle_energy;
  // angle per edge of each neighboring triangle
  thrust::device_vector<double> initial_angle;
  // Temporary vectors used for computation. 
	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};

// Struct to store information related to linear springs (edges).
struct LinearSpringInfoVecs {

	double apical_spring_constant = 5.0;
  double basal_spring_constant = 5.0;
  double spring_constant_vertical = 5.0; // spring constant for the vertical springs in the double layer eversion model. 
  
  int factor = 2; // Used for reduction.
	double spring_constant = 5.0;
	double spring_constant_weak = 5.0;
	double spring_constant_att1 = 5.0;
	double spring_constant_att2;
	double spring_constant_rep1; // This is the "D" in Morse potential
	double spring_constant_rep2; // This is the "a" in Morse potential

	double linear_spring_energy;
	double memrepulsion_energy;
	double scalar_edge_length;
	
	thrust::device_vector<double> edge_initial_length;
  thrust::device_vector<double> edge_rest_length;
  thrust::device_vector<double> edge_final_length;

  // Temporary vectors used for computation. 
	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};

// Struct to store general parameters of the system.
struct GeneralParams{
	std::vector<int> edge_undergoing_growth;
	bool nonuniform_wall_weakening_bend;
	bool nonuniform_wall_weakening_linear;
	bool nonuniform_wall_weakening_area;
	bool nonuniform_wall_weakening;
	double maxSpringScaler_linear = 1.0; 
	double maxSpringScaler_area = 1.0;
	double maxSpringScaler_bend = 1.0;
	double ratio_for_HillFunctionStiffness;
  double current_quasi_step;
  double Tf; //= 15.0;
  double gamma = 0.1; //friction
  // growth / bending options
 // bool useEdgeProjection      = true;  // lambda·e projection instead of lambda_rr average
  bool enableBendAnisotropy   = true;  // spontaneous-curvature coupling
  double thickness        = 0.10;    // t  in  C0 = (lambda ff – lambda rr)/t
  double tol;// = 1e-4;//tol of 1e-5 and 1e-6 are too low. The positions begin to oscillate and there is no convergence. So we stick to 1e-4. I'm now going to 1e-9. This works but takes too long so I am submitting a different sim with 1e-4
  
  int num_layers;
  double theta_DV = 0.1002;
	std::vector<int> triangle_undergoing_growth;
	double chemdiff_time_step_size;
	int current_total_sim_step;
	double chemdiff_max_step;
	double current_bud_area;
	double kT;
	double kT_growth;
	double tau;
	int solve_time=100;
	//double Rmin = 100.0; // Current mesh minimum edge length, subject to change.
	double Rmin_growth;
	double abs_Rmin;
	int iteration = 0;
	int maxNodeCount;
	int maxNodeCountLJ;
	double length_scale;
	//parameters for advancing timestep and determining equilibrium

	double dt;// = 0.01; // nav will change dt to 0.001 later. 
  double dx;
	double nodeMass = 1.0;
	int edges_in_upperhem_list_length;
	double growth_energy_scaling;
	thrust::device_vector<int> edge_to_ljparticle;
	thrust::device_vector<int> nodes_in_tip;
	thrust::device_vector<int> nodes_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem_list;
	thrust::device_vector<int> edges_in_tip;
	thrust::device_vector<int> triangles_in_tip;
	thrust::device_vector<int> triangles_in_upperhem;
	thrust::device_vector<int> boundaries_in_upperhem;
  thrust::device_vector<int> boundaries_in_lowerhem;
	thrust::device_vector<double> angle_per_edge;
  thrust::device_vector<int> nodes_in_DV;
  thrust::device_vector<int> edges_in_DV;
  
 
	double centerX = 0.0;
	double centerY = 0.0;
	double centerZ = 0.0; // has been modified in System.cu to reflect the center of the layer currently being computed. 
 
  double c_dx = 0.0; 
  double c_dy = 0.0; 
  double c_dz = -15.239; // center of sphere that the disc is a part of. (calculated from outer surface) //redundant now 12_22_25
  double apical_Rx;
	double current_total_volume;
	double true_current_total_volume;
	double eq_total_volume;
	double volume_spring_constant = 5.0;
	double volume_energy;
	double eq_total_boundary_length;
	double line_tension_energy;
	double line_tension_constant;
	double safeguardthreshold;
	

	//int num_of_triangles;
	//int num_of_edges;
	int true_num_edges;

	double insertion_energy_cost;
	double strain_threshold;
	double strain_threshold2;

	int SCALE_TYPE;
	double scaling_pow;
	double gausssigma;
	double hilleqnconst;
	double hilleqnpow;

	thrust::device_vector<int> no_weakening;
	double septin_ring_z;
	double boundary_z;
 
  //double epsilon_r_center;
  //double epsilon_r = -0.12406004; // change this name to Isotropic strain
  //double epsilon_t = 1.20789496;// change this name to Anisotropic strain
  
  double DV_ax, DV_ay, DV_az, DV_bx, DV_by, DV_bz;
  double ODx, ODy, ODz, OVx, OVy, OVz;
 
  thrust::device_vector<double> lambda_iso_center_outDV_v;

    
  thrust::device_vector<double> lambda_iso_edge_outDV_v;

  
  thrust::device_vector<double> lambda_aniso_center_outDV_v;

  
  thrust::device_vector<double> lambda_aniso_edge_outDV_v;
  
  
  
  // inDV
  thrust::device_vector<double> lambda_iso_center_inDV_v;

  
  thrust::device_vector<double> lambda_iso_edge_inDV_v;

  
  thrust::device_vector<double> lambda_aniso_center_inDV_v;

  
  thrust::device_vector<double> lambda_aniso_edge_inDV_v;


  // Strain field (lambda) values in the outDV region at different stages of eversion.
  double lambda_iso_center_outDV = -0.12406004; //   wl3-0hAPF (-0.12406004) | wl3-2hAPF (-0.29431527) | wl3-4hAPF (-0.20050286)
  double lambda_iso_edge_outDV = 1.20789496; //      wl3-0hAPF ( 1.20789496) | wl3-2hAPF ( 1.44256030) | wl3-4hAPF ( 1.43479468)
  double lambda_aniso_center_outDV = -0.06172103;//  wl3-0hAPF (-0.06172103) | wl3-2hAPF ( 0.21823128) | wl3-4hAPF ( 0.29444448)
  double lambda_aniso_edge_outDV = 1.01053997; //    wl3-0hAPF ( 1.01053997) | wl3-2hAPF ( 0.98874841) | wl3-4hAPF ( 0.97652462)
  
  // Strain field (lambda) values in inDV region at different stages of eversion. 
  double lambda_iso_center_inDV = -0.09848994; //    wl3-0hAPF (-0.09848994) | wl3-2hAPF (-0.11692544) | wl3-4hAPF (-0.06151876)
  double lambda_iso_edge_inDV = 1.18401136;//       wl3-0hAPF ( 1.18401136) | wl3-2hAPF ( 1.21007540) | wl3-4hAPF ( 1.47472744)
  double lambda_aniso_center_inDV = -0.12904887; //  wl3-0hAPF (-0.12904887) | wl3-2hAPF (-0.21271504) | wl3-4hAPF (-0.30567972)
  double lambda_aniso_edge_inDV = 1.03128453; //    wl3-0hAPF ( 1.03128453) | wl3-2hAPF ( 1.24178074) | wl3-4hAPF ( 1.29370391)
  
  double disc_radius;//= 77.66; // manually input radius for each of the discs you're computing. 
  double sphere_R;
  double theta_max=0.87226;
  thrust::device_vector<double> rho;
  
 // double disc_radius;
  

};

__host__ __device__
inline Mat_3x3 makeIdentity3x3()
{
    return Mat_3x3(
        CVec3(1.0, 0.0, 0.0),
        CVec3(0.0, 1.0, 0.0),
        CVec3(0.0, 0.0, 1.0)   // ? no trailing comma
    );
}

__host__ __device__
inline CVec3 makeOnes3() { return CVec3(1.0,1.0,1.0); }



class Storage;
class SystemBuilder;
struct HostSetInfoVecs;

// Class representing the entire system.
class System {
public:
	std::weak_ptr<SystemBuilder> weak_bld_ptr;
	GeneralParams generalParams;
	DomainParams domainParams;
	AuxVecs auxVecs;
	CoordInfoVecs coordInfoVecs;
  PrismInfoVecs prismInfoVecs;

	CapsidInfoVecs capsidInfoVecs;
	LinearSpringInfoVecs linearSpringInfoVecs;
	BendingTriangleInfoVecs bendingTriangleInfoVecs;
	AreaTriangleInfoVecs areaTriangleInfoVecs;
	LJInfoVecs ljInfoVecs;
	// RBC rbc;
	// ELEMENT element;

	std::shared_ptr<Storage> storage;

	//gsl_vector* df;
	//gsl_vector* locations;
    



public:

	System();

	void set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr);
 
  void printEdges();
  
  void printTriangles();

	void PrintForce();

	void initializeSystem(HostSetInfoVecs& hostSetInfoVecs);

	void assignStorage(std::shared_ptr<Storage> _storage);
 
	void solveSystem();

	void setExtras();

	void setBucketScheme();
  
	//void Solve_Forces(const gsl_vector* temp_locations);
	void Solve_Forces();
	//double Solve_Energy(const gsl_vector* temp_locations);
  
	//void dev_to_gsl_loc_update(gsl_vector* temp_locations);
	//void gsl_to_dev_loc_update(const gsl_vector* temp_locations);
	
  //void relaxUntilConverged(int &k);

};



#endif /*POLYMERSYSTEM_H_*/
