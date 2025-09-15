/*

Well you didnt save our last rant did you. Good on you!!
Sarcasm aside..HOW COULD YOU. WE WORKED HARD ON THAT. DAMMIT NAV.

SO here is where I am rn. adding a verlet velocity advancer for energy conservation to my current node advance file. Here is what I have done.

I have added a new vector nodeVelX/Y/Z in system.h and systemstructures.h. coordInfoVecs and hostSetInfoVecs respectively. I have then gone into system.cu and initialized it. Idk if there are other places that this parameter needs to be updated. I'll change it as it goes along. Might be happening in utilities.cpp. Lez see

Okay so here's where we're at. We have the double layered spring model being fed into the tensor and computed but for some reason all the springs are being pulled into each other.
Need to figure out why this is happening right after relaxation. Should not be the case but it is.

Next thing to do is to separate the vertical spring constants from the rest of the mesh computation. Maybe treat them like the septin ring ones? Let's start work on that.

Things to keep in mind:

Linear Springs - the initial length is the edge_initial_length array and the final length it has to go to is the edge_rest_length. So when we're starting out we dont want any linear energy acting
the springs. Once the strain tensor is applied, we will then have the lengths be perturbed. So modify that function accordingly.

In the Science Advances paper they have used a CSV to input the strain parameters. Maybe we need to be doing the same thing.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// NEW WORK FLOW ////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


1. Start by putting vertical spring lengths into the septin ring nodes.


2. MAIN ALGORITH FLOW:
  a. edge_initial_lengths -> initial lengths of the edges. First start by calculating the final FINAL resting lengths through application of the strain tensor OR
     Start by setting a total time for simulation. Use that to get the total number of steps. Then use that to figure out final rest length from the strain tensor. Store this in a new vector called
     edge_final_length.
     NB - TOP AND BOTTOM LAYERS HAVE THEIR REST LENGTHS UPDATED BY THE SAME AMOUNT IN THE ACTIVE SHAPE PROGRAMMING PAPER.

     Now when the strain tensor runs it'll be broken down into steps and incrementally run so it'll be - edge_rest_length(T) = edge_initial_length + (edge_final_length - edge_initial_length) *T/Tf
     Here Tf is the quasi static number of steps the simulation takes place in.

  b. Here's how the energy is going to change: linear_spring_energy = 1/2 sum_over_all_edges ( k * (edge_current_length - edge_rest_length)^2. Overdamped dynamics is where edge_current_length comes
     into play.

  c. Energy must now be minimized. Look at notes for equations.

     Apply first strain increment - solve forces - quasi statically search for the minimum energy using tau - advance node positions - repeat.

  Tf - edge_final_length
  T - edge_rest_length
  tau - edge_current_length
  T_0 - edge_initial length

  okay most of the above things are done.

  Need to:
  1. Look into LinearSprings files
  2. Recenter mesh initially so that it does not have the initial mesh in a separate position.
  3. Modify the strain tensor so its only working on one layer at a time.
  4. Figure out why the vertical springs have so much tension and collapse.



  Okay so here's where we're at right now:

  I'm trying to figure out why my vertical springs are not being calculated properly. I think I'll need to verify the data structure and see which edges are connected to which ones then come up with an algorithm that would find those nodes and calculate the distance between them. This would be added in systembuilder.

  Oh and also currently all my constants are being set as the default ones and not taken from the data structure so that's also a bummer.

*/
/*

Well you didnt save our last rant did you. Good on you!!
Sarcasm aside..HOW COULD YOU. WE WORKED HARD ON THAT. DAMMIT NAV.

SO here is where I am rn. adding a verlet velocity advancer for energy conservation to my current node advance file. Here is what I have done.

I have added a new vector nodeVelX/Y/Z in system.h and systemstructures.h. coordInfoVecs and hostSetInfoVecs respectively. I have then gone into system.cu and initialized it. Idk if there are other places that this parameter needs to be updated. I'll change it as it goes along. Might be happening in utilities.cpp. Lez see

Okay so here's where we're at. We have the double layered spring model being fed into the tensor and computed but for some reason all the springs are being pulled into each other.
Need to figure out why this is happening right after relaxation. Should not be the case but it is.

Next thing to do is to separate the vertical spring constants from the rest of the mesh computation. Maybe treat them like the septin ring ones? Let's start work on that.

Things to keep in mind:

Linear Springs - the initial length is the edge_initial_length array and the final length it has to go to is the edge_rest_length. So when we're starting out we dont want any linear energy acting
the springs. Once the strain tensor is applied, we will then have the lengths be perturbed. So modify that function accordingly.

In the Science Advances paper they have used a CSV to input the strain parameters. Maybe we need to be doing the same thing.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// NEW WORK FLOW ////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


1. Start by putting vertical spring lengths into the septin ring nodes.


2. MAIN ALGORITH FLOW:
a. edge_initial_lengths -> initial lengths of the edges. First start by calculating the final FINAL resting lengths through application of the strain tensor OR
   Start by setting a total time for simulation. Use that to get the total number of steps. Then use that to figure out final rest length from the strain tensor. Store this in a new vector called
   edge_final_length.
   NB - TOP AND BOTTOM LAYERS HAVE THEIR REST LENGTHS UPDATED BY THE SAME AMOUNT IN THE ACTIVE SHAPE PROGRAMMING PAPER.

   Now when the strain tensor runs it'll be broken down into steps and incrementally run so it'll be - edge_rest_length(T) = edge_initial_length + (edge_final_length - edge_initial_length) *T/Tf
   Here Tf is the quasi static number of steps the simulation takes place in.

b. Here's how the energy is going to change: linear_spring_energy = 1/2 sum_over_all_edges ( k * (edge_current_length - edge_rest_length)^2. Overdamped dynamics is where edge_current_length comes
   into play.

c. Energy must now be minimized. Look at notes for equations.

   Apply first strain increment - solve forces - quasi statically search for the minimum energy using tau - advance node positions - repeat.

Tf - edge_final_length
T - edge_rest_length
tau - edge_current_length
T_0 - edge_initial length

okay most of the above things are done.

Need to:
1. Look into LinearSprings files
2. Recenter mesh initially so that it does not have the initial mesh in a separate position.
3. Modify the strain tensor so its only working on one layer at a time.
4. Figure out why the vertical springs have so much tension and collapse.



Okay so here's where we're at right now:

I'm trying to figure out why my vertical springs are not being calculated properly. I think I'll need to verify the data structure and see which edges are connected to which ones then come up with an algorithm that would find those nodes and calculate the distance between them. This would be added in systembuilder.

Oh and also currently all my constants are being set as the default ones and not taken from the data structure so that's also a bummer.

06/06/25:

Okay so quals are 6 days away. less than that actually. You can do it! Good luck

In the meantime let's look at the plan beyond quals. 

1. Make sure strain tensor is correct. I have my doubts about it right now. 
2. Gradient descent takes fricking AGES. Figure out why this is. It really should not be. 
3. Remove the scale parameter settings for hill const diff for wall weakening from your code. 
4. This is the last one. And a big one. START MAKING YOUR 3D SCE model. We'll figure out couplings later. But for now start reading the papers on it.



06/23/25 

Okay so the month of June is pretty much over now. Passed the quals and now I'm just drifting. Where do I go from here. I guess a lot of people feel this way after a taxing run where they 
do a long stretch and burn themselves out. I know I'm not a significant part of this lab and so I will keep working quietly on the things that are assigned me. I wont let it break my spirit tho. 


From tomorrow I will start working on transfering this code to python and seeing if we can have a real time updating screen for the visuals. We'll also do parallelization in python for it. Will be 
another really big task but I have high hopes for the summer. So lez see. 

Email Khoi. Ask him for help with this. 


/////////////////////////////////////////////////////

Alright we're back again. This time for real 

07/06/2025

Let's see now. 

We've implemented strain tensor but we do not know if it is correct. So go over that once more. 

BUT FIRST
 
Start by separating DV boundary. 
 Implement Strain tensor on DV boundary. 
 Then check strain tensor. 
 Alongside, make sure to set up the local version of the running code 
 Look up how to model thick tissues with mechanical properties. 
 
 
 Also know that you're not getting any simulation results right now which means you've potentially deleted a file that could have been useful or a thing in System.cu that could have been useful. So let's do some debugging. 

*/

#include <stdio.h>
#include "System.h"
#include "SystemStructures.h"
#include "AreaTriangles.h"
//#include "BendingTriangles.h"
//#include "MemRepulsionSprings_universal.h"
//#include "MemRepulsionSprings_local.h"
//#include "MemRepulsionEnergy.h"
#include "LinearSprings.h"
#include "NodeAdvance.h"
#include "Storage.h"
#include "Utilities.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"
#include <bits/stdc++.h>
#include "LineTensionSprings.h"
#include <math.h>
#include <list>
//#include "TurgorForce.h"
#include "LinearSpringsEnergy.h"
#include "StrainTensor.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include "gradientRelax.h"

// somehow the gradient is not being set in my version - Kevin

// Helper function to count elements greater than or equal to zero in a vector.
int count_bigger(const std::vector<int> &elems)
{
    return std::count_if(elems.begin(), elems.end(), [](int c)
                         { return c >= 0; });
}

// Constructor for the System class.
System::System() {};

// Print net force on nodes along a radial line (? ˜ 0) from disc center to boundary
void System::PrintForce() {
    // Copy device forces to host
    thrust::host_vector<double> h_fx = coordInfoVecs.nodeForceX;
    thrust::host_vector<double> h_fy = coordInfoVecs.nodeForceY;
    thrust::host_vector<double> h_fz = coordInfoVecs.nodeForceZ;
    thrust::host_vector<double> h_x  = coordInfoVecs.nodeLocX;
    thrust::host_vector<double> h_y  = coordInfoVecs.nodeLocY;

    const double desiredTheta = 0.0;    // along +x axis
    const double eps = 0.01;            // angular tolerance (rad)
    std::vector<std::pair<double,int>> picks;
    int N = static_cast<int>(h_x.size());
    for (int i = 0; i < N; ++i) {
        double r = std::hypot(h_x[i], h_y[i]);
        double theta = std::atan2(h_y[i], h_x[i]);
        double diff = std::fabs(theta - desiredTheta);
        if (diff > M_PI) diff = 2*M_PI - diff;
        if (diff < eps) picks.emplace_back(r, i);
    }
    std::sort(picks.begin(), picks.end());

    std::printf("   r      Fx       Fy       Fz\n");
    for (auto &pr : picks) {
        int idx = pr.second;
        std::printf("%6.3f  %7.3e  %7.3e  %7.3e\n",
                    pr.first,
                    h_fx[idx], h_fy[idx], h_fz[idx]);
    }
}

// Function to solve the forces in the system.
void System::Solve_Forces()
{

    // Reset all forces to zero.
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

    // Compute forces and energy due to linear springs.
    ComputeLinearSprings(
        generalParams,
        coordInfoVecs,
        linearSpringInfoVecs,
        ljInfoVecs);

    // Compute forces and energy due to area springs. Nav commented out to test Active shape programming mesh type 02/27/2025  . Put back in 03/23/25
//      	ComputeAreaTriangleSprings(
//      		generalParams,
//      		coordInfoVecs,
//      		areaTriangleInfoVecs);

    // Compute forces and energy due to turgor pressure springs. (nav - commenting these out for now for flat surface 5/29/24) nav reintroducing the turgor pressure because the eversion wing does have turgor pressure. 8/17/2024
    // ComputeTurgorSprings(
    // generalParams,
    // coordInfoVecs,
    // areaTriangleInfoVecs
    //);

    // Compute forces and energy due to bending springs. Turn this off 10/10/24
//      	ComputeCosTriangleSprings(
//      		generalParams,
//      		coordInfoVecs,
//      		bendingTriangleInfoVecs);

    // Compute forces and energy due to membrane repulsion springs.// Nav commented out to test Active shape programming mesh type 02/27/2025. PUt back in 03/23/25
//      	ComputeMemRepulsionSprings_local(
//      		coordInfoVecs,
//      		linearSpringInfoVecs,
//      		capsidInfoVecs,
//      		generalParams,
//      		auxVecs);

    // Compute forces and energy due to volume springs. //(nav - commenting these out for now for flat surface 5/29/24) Nav had uncommented but she's bringing the comment back because testing out Active shape mesh 02/27/25
//      	ComputeVolume(
//      		generalParams,
//      		coordInfoVecs,
//      		linearSpringInfoVecs,
//      		ljInfoVecs);
    
    // Now print forces along the radial line
   // PrintForce();
};

// Function to solve the entire system.
void System::solveSystem()
{

    // Nav - I dont want to remove these variables. These may come in handy.
    // coordInfoVecs.k_0 = 20.0;
    // coordInfoVecs.k_1 = 25.0;
    // coordInfoVecs.k_2 = 5.0;
    // coordInfoVecs.k_3 = 5.0;
    // coordInfoVecs.k_4 = 1.0;
    // coordInfoVecs.k_ss = 12;//10.75;
    // coordInfoVecs.beta = 1.0/1.0;///1.45;
    // coordInfoVecs.gamma = 1.0;
    // coordInfoVecs.q1 = 10.0;
    // coordInfoVecs.h = 10.0;
    
    generalParams.dt = 0.01;
    
    
    std::cout<<" Simulating wL3 - 0hAPF "<<std::endl;
    std::cout<<"lambda_iso_center_outDV = -0.12406004 - "<<generalParams.lambda_iso_center_outDV<<std::endl;
    Solve_Forces();
    uint mem_prealloc = 4; // Make sure that this number is the same as set in System::initializeSystem found near the bottom of this script. - Kevin. Q. why is this the case? Why does it need to be the same? - Nav.

    // Initial values for determining the region of material insertion.
    double current_edge_to_tip_height_scale = INT_MAX; // 2.0; // The maximum edge to tip height scale, initially set to INT_MAX.
    std::cout << "current_edge_to_tip_height_scale = " << current_edge_to_tip_height_scale << std::endl;

    // Determines how far away from the tip can new material be inserted during edge swap.
    double current_edge_to_tip_height_scale_ES = INT_MAX; // 4.0;//2.0; // Initially set to INT_MAX
    std::cout << "current_edge_to_tip_height_scale for edge-swap = " << current_edge_to_tip_height_scale_ES << std::endl;

    // Determines how far away from the tip can new material be inserted.
    double bdry_to_tip_height_scale = INT_MAX; // 4.0; // Initially set to INT_MAX
    std::cout << "bdry_to_tip_height_scale = " << bdry_to_tip_height_scale << std::endl;

    // Boolean flag to determine if restiffening is being simulated.
    bool isRestiffening = true; // false; // change this to true (nav) - testing for FvK. was false but I have changed to true - nav 8/5/24
    std::cout << "Are we simulating a case where restiffening or the restoration (even just partially) of mechanical properties? " << isRestiffening << " (bool) " << std::endl;

    // Scaling factors for restiffening regions.
    double scale_linear_restiff;
    double scale_bend_restiff;
    double scale_area_restiff;

    // If restiffening is enabled, calculate the scaling factors.
    if (isRestiffening == true)
    {
        // Scale for linear springs in restiffening regions.
        scale_linear_restiff = linearSpringInfoVecs.spring_constant * 0.1; // 0.25;//25.0/2.5;//75.0/15.0; nav changing it to 0.1 8/5/24
        // Scale for bending springs in restiffening regions.
        scale_bend_restiff = bendingTriangleInfoVecs.spring_constant * 0.1; // 0.05;//10.0/1.0;//75.0/7.5; nav changing it to 0.1 8/5/24
                                                                            //  Scale for area springs in restiffening regions.
        scale_area_restiff = areaTriangleInfoVecs.spring_constant * 0.1; // 0.25;//50.0/5.0;//75.0/15.0; nav changing it to 0.1 8/5/24

        std::cout << "restiff region linear = " << scale_linear_restiff << std::endl;
        std::cout << "restiff region bend = " << scale_bend_restiff << std::endl;
        std::cout << "restiff region area = " << scale_area_restiff << std::endl;
        std::cout << "If restiff value is higher than that of the weakened case, something is wrong" << std::endl;
    }

    // Flag to determine if the system has been triggered.
    bool triggered = false;

    // Initialize the current total simulation step.
    generalParams.current_total_sim_step = 0;

    // Maximum relaxation steps before growth and edge swap.
    int relax_max_steps_before_growth_and_edgeswap = 3e3;
    int relax_max_steps_before_strain_tensor_and_gradient_descent = 3e3;
    // std::cout<<"relax max steps before growth and edgeswap = "<<relax_max_steps_before_growth_and_edgeswap<<"*max_runTime"<<std::endl;
    std::cout << "relax max steps before strain tensor application and gradient descent = " << relax_max_steps_before_strain_tensor_and_gradient_descent << std::endl;

    // Create a shared pointer for Utilities.
    auto utilities_ptr = std::make_shared<Utilities>(coordInfoVecs, generalParams);

    // Create a shared pointer for SystemBuilder.
    auto build_ptr = weak_bld_ptr.lock();
    std::cout << "Declaration of rbc and n_rbc complete." << std::endl;
    std::cout << "Utilities_ptr declaration complete." << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Build the "u" vector representing the external or internal influencer for polarization /////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Part 2

    // Variables to store the maximum and minimum heights of cell triangles.
    double max_cell_triangle_height, min_cell_triangle_height, v1, v2, v3, v4, cell_height;

    // Initialize the maximum and minimum cell triangle heights to extreme values.
    max_cell_triangle_height = -10000.0;
    min_cell_triangle_height = 10000.0;

    // Loop over all triangles to calculate the heights and find the maximum and minimum heights
    for (int i = 0; i < coordInfoVecs.num_triangles; i++)
    {

        // Skip trinagles with invalid node indices.
        if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX - 100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX - 100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX - 100))
        {
            continue;
        }
        if (coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX + 100) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX + 100) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX + 100))
        {
            continue;
        }

        // Get the Z coordinates of the nodes forming the current triangle.
        v1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
        v2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
        v3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];

        // Calculate the average Z coordinate of the triangle (cell triangle height).
        v4 = (v1 + v2 + v3) / 3.0;
        // std::cout<<"Max triangle heights (z-coordinates) = "<<v4<<std::endl;

        // Update the maximum and minimum cell triangle heights.
        if (v4 >= max_cell_triangle_height)
        {
            max_cell_triangle_height = v4;
        }
        if (v4 <= min_cell_triangle_height)
        {
            min_cell_triangle_height = v4;
        }
    }

    // Calculate the cell height as the difference between the maximum and minimum cell triangle heights.
    cell_height = (max_cell_triangle_height - min_cell_triangle_height);

    std::cout << "Max_cell_triangle_height = " << max_cell_triangle_height << ", min_cell_triangle_height = " << min_cell_triangle_height << std::endl;

    std::cout << "Cell height = " << cell_height << std::endl;

    // Message printed when max and min height of triangles has been determined.
    std::cout << "Determination of max triangle and min cell height complete." << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Part 3

    // Setting various simulation parameters and initializing variables.

    // Set the node mass for the simulation.
    generalParams.nodeMass = 1.0;

    // Initialize growth counter to keep track of growth events.
    int GROWTH_COUNTER = 0;

    // Set the minimum number of edge loops for edge-swap during growth events.
    int min_num_edge_loop = 1;
    std::cout << "min_num_edge_loop for edgeswap = " << min_num_edge_loop << std::endl;

    // DEfine max volume and bud area ratios for growth control.
    double MAX_VOLUME_RATIO = 2.0;
    double MAX_BUD_AREA_RATIO = 100.0;
    double MAX_APICAL_AREA_RATIO = 100.0;

    // Set the maximum number of growth events per growth cycle. (INT_MAX means unlimited)
    int MAX_GROWTH_PER_GROWTH_EVENT = 1; // INT_MAX;
    std::cout << "MAX_GROWTH_NUMBER (# of edge to expand) per growth event = " << MAX_GROWTH_PER_GROWTH_EVENT << std::endl;

    // Set the frequency of growth events (how many times Max_Runtime has to be reached to perform growth).
    int GROWTH_FREQUENCY = 25; // 150;//100;//95;//70;//25*3; // E.g., if Max_Runtime = 100, growth will occur every 25 time units.
    std::cout << "GROWTH_FREQ (how many times Max_Runtime has to be reached to perform growth) = " << GROWTH_FREQUENCY << std::endl;

    // Set growth frequency of growth events for variable edge-swap rate cases.
    int GROWTH_FREQUENCY2 = 25; // 150;//100;//95;//70;//25*3;
    std::cout << "GROWTH_FREQ2 (how many times Max_Runtime has to be reached to perform growth, for variable ES rate cases" << GROWTH_FREQUENCY2 << std::endl;

    // Set the point of transition for growth events (describes the ratio of the total simulation time for relaxation (edge-swap) frequency to change).
    double pointOfTransition = 0.10;
    std::cout << "Point of transition describes the ratio of the total time simulation time is reached for relaxation (edge-swap) frequency to change : " << pointOfTransition << std::endl;
    std::cout << "Point of transition is also used to indicate if strain_threshold needs to change or not" << std::endl;

    // Set the energy gradient threshold for growth events.
    double energy_gradient_threshold = 0.02; // 0.01; // Threshold used to trigger growth events based on energy gradients.
    std::cout << "ENERGY_GRADIENT_THRESHOLD = " << energy_gradient_threshold << std::endl;

    // Set kT_growth value for growth events.
    generalParams.kT_growth = 1.0;

    // Set the SCALE_TYPE for weakening during growth events. (0 to 4)
    // 0:= Gaussian-like weakening
    // 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening
    // 2:= pure Gaussian weakening
    // 3:= isotropic
    // 4:= hill equation
    // Note that (3) is used in combination with sigma = INT_MAX;
    generalParams.SCALE_TYPE = 3; // Original scale type was 3. Nav changed it to 0 for flat code. 6/2/24. \\ 0 did not work so nav changed it back to 3 8/18/24
    std::cout << "SCALE TYPE = " << generalParams.SCALE_TYPE << std::endl;
    std::cout << "0:= sigmoidal Gaussian-like weakening, 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening, 2:= pure Gaussian weakening, 3:= isotropic, 4:= hill equation" << std::endl;

    // Check and set additional parameters based on SCALE_TYPE.
    if (generalParams.SCALE_TYPE == 1)
    {
        generalParams.scaling_pow = 2.0;
        std::cout << "scaling_pow (this is for SCALE_TYPE = 1 case) = " << generalParams.scaling_pow << std::endl;
    }
    if (generalParams.SCALE_TYPE == 0)
    {
        generalParams.gausssigma = 0.1;
        std::cout << "gausssigma (this is for the SCALE_TYPE = 0 case) = " << generalParams.gausssigma << std::endl;
    }

    // Set the display_token to true for displaying additional information during growth events.
    bool display_token = true;

    // Declare variables for Hill function_dependent wall stiffness.
    double dtb_scaler, targetHillEqnPow;
    if (generalParams.SCALE_TYPE == 4)
    {
        generalParams.ratio_for_HillFunctionStiffness = 4.0;
        std::cout << "Hill function dependent wall stiffness triggers when the the distance between tip of the bud and the septin ring is " << generalParams.ratio_for_HillFunctionStiffness << std::endl;
        std::cout << "times larger than the equilibrium length Rmin" << std::endl;
        dtb_scaler = 1.0;
        targetHillEqnPow = 16.0;
        std::cout << "The EC50 position is scaled by " << dtb_scaler << " on the distance from tip to boundary, hence the EC50 occurs on dtb*" << dtb_scaler << "/dtb_max" << std::endl;
        std::cout << "Target hill equation power = " << targetHillEqnPow << std::endl;
    }

    // coordInfoVecs.scaling_per_edge.
    // generalParams.hilleqnconst = 0.9;
    // generalParams.hilleqnpow = 40.0;

    // Declare vectors to store nodes, triangles and edges in the growth region. (nav - unsure if I need to take these out since we will not be performing growth)
    std::vector<int> nodes_in_growth;
    std::vector<int> triangles_in_growth;
    std::vector<int> edges_in_growth;

    // Declare variables for distance to boundary and maximum distance used in Hill equation.
    double dtb;     // dtb := distance to boundary
    double dtb_max; // dtb_max := the max distance used to calculate the distance ratio in the Hill equation.

    // Declare variables for sigma (for gradient distribution variance) and sigma_true (for Gaussian-related distribution variance) used in SCALE_TYPE 0.
    double sigma, sigma_true;

    if (generalParams.SCALE_TYPE == 0)
    {
        sigma = 0.0;            // INT_MAX; //if this is set to be INT_MAX then we assume isotropic weakening.
        sigma_true = sqrt(0.5); // This is the variance used to calculate the scaling of the wall weakening.
        std::cout << "initial sigma (for gradient distribution variance), based on initial distribution of Cdc42, if using true gaussian weakening = " << sigma << std::endl;
        std::cout << "If sigma = INT_MAX, then we have isotropic weakening scenario" << std::endl;
        std::cout << "true sigma (for gaussian-related distribution variance) = " << sigma_true << std::endl;
    }

    // Set the insertion energy cost for material insertion during growth events.
    generalParams.insertion_energy_cost = -log(0.0025); // why is material insertion energy cost -log(0.0025) (question) Ask dr. Chen
    std::cout << "GROWTH: material insertion energy cost (dependent on local chemical concentration) = " << generalParams.insertion_energy_cost << std::endl;

    // Set the strain thresholds for the martial (material?) insertion probability calculation during growth events.
    double strain_threshold1 = 0.05; // 0.01; // Set strain_threshold for initial calculation.
    double strain_threshold2 = 0.05; // Set strain_threshold for subsequent changes if needed.
    generalParams.strain_threshold = strain_threshold1;
    std::cout << "GROWTH: critical strain threshold used for insertion probability calculation = " << generalParams.strain_threshold << ", value loaded = " << strain_threshold1 << std::endl;
    std::cout << "GROWTH: critical strain threshold used for insertion probability calculation if changes are needed= " << strain_threshold2 << std::endl;

    // Set the growth energy scaling for material insertion probability durng growth events.
    generalParams.growth_energy_scaling = 1.0; // 0.01375; // Set the scaling factor for growth energy.
    std::cout << "GROWTH ENERGY SCALING FOR MATERIAL INSERTION PROBABILITY = " << generalParams.growth_energy_scaling << std::endl;

    // Set the neighbour safeguardthreshold (the max number of neighboring nodes a node can have).
    generalParams.safeguardthreshold = 9;
    std::cout << "NEIGHBOR SAFE GUARD THRESHOLD = " << generalParams.safeguardthreshold << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// PARAMETER SETTINGS ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    // Part 4

    // Setting various simulation parameters and initializing variables.

    // Initialize total growth counters to keep track of growth events.
    int TOTAL_GROWTH_COUNTER = 0;
    int TOTAL_GROWTH_ATTEMPT = 0;

    // Set the maximum runtime for the simulation (in time units)
    double Max_Runtime = generalParams.dt * 50.0; //*50.0; nav commented this out for double sheet testing. 9/16/2024. //generalParams.dt*10.0;//generalParams.dt*50.0; Nav commented out this last part for flat code 6/5/24. //50.0; time units //time step is here. Nav has once again commented out the __.dt*10 to make it run longer. 8/2/24

    std::cout << "dt = " << generalParams.dt << std::endl;
    double minimal_run_time_ratio = 1.0;                 // Not used in this section.
    double Max_RunStep = Max_Runtime / generalParams.dt; // Calculate the maximum runstep based on the max runtime.
    std::cout << "Max runtime = " << Max_Runtime << std::endl;
    std::cout << "Max runstep = " << Max_RunStep << std::endl;

    // Boolean flag to determine if the simulation should continue running.
    bool runSim = false; // turn off the budding simulation and edgeswap

    // Declare and initialize variables for growth-related calculations.
    int num_edge_loop;
    double initial_kT;
    initial_kT = generalParams.kT; // Stores the initial kT value for the acceptance of changes after looping through every edge within proximity.

    // double SAMPLE_SIZE = 0.05;//0.025;
    // std::cout<<"Sample ratio: "<<SAMPLE_SIZE<<std::endl;
    // std::cout<<"If the Sample raio is 0, it means we have chosen a fixed number of attempt throughout the simulation"<<std::endl;

    // Set the sample size for testing edges during bondflip remeshing.
    double SAMPLE_SIZE = 2;
    std::cout << "Sample size: " << SAMPLE_SIZE << std::endl;

    // Set the record frequency for the dsaving simulation data (time steps).
    int RECORD_TIME = 1; // round(Max_RunStep/2); Save data every time step.
    std::cout << "Record frequency = " << RECORD_TIME << std::endl;

    // translate_frequency determines the frequency for the mesh to re-center and perform dynamical remeshing.
    int translate_frequency = 10;
    std::cout << "recentering of the model cell frequency = " << translate_frequency << std::endl;

    // Set the total number of growth events and targeted growth events for the simulation.
    int NUMBER_OF_GROWTH_EVENT = 1000;                        // changed by nav on 03/04/2025//200;//2000 - nav changed this for flat code. 6/2/24;//1000;//1000*2; // Total number of growth events. // Nav once again changed it to 1000 from 200. 8/26/24
    int NUMBER_OF_TARGETED_GROWTH_EVENT = 1;                  // Number of targeted growth events.
    int NKBT = GROWTH_FREQUENCY * NUMBER_OF_GROWTH_EVENT * 2; // nav changed this 03/04/2025//10; // GROWTH_FREQUENCY*NUMBER_OF_GROWTH_EVENT*2; Nav changed this last one for the flat code. 6/5/24.//GROWTH_FREQUENCY*NUMBER_OF_GROWTH_EVENT;//10000;//7500; Nav is now changing it back from 10 to turn growth on. 8/26/24
    std::cout << "Number of edge-swap per kBT value (or total number of edge-swap if kBT is fixed) = " << NKBT << std::endl;

    int GROWTH_FREQUENCY_SCALE = 4;
    std::cout << "GROWTH FREQ SCALE: decides how many growth event must be checked before recording the result = " << GROWTH_FREQUENCY_SCALE << std::endl;

    double min_kT = -0.1; // 0.21;
    std::cout << "min kT for simulation termination = " << min_kT << std::endl;

    // Initialize WHEN for conditional checks.
    int WHEN = 0;

    // Initialize the following for energy calculations.
    double old_total_energy = 0.0;
    double new_total_energy = 0.0;
    double energy_gradient = 0.0;
    double energy_rep = 0.0;

    // Initialize the number of simulation steps run.
    int Num_of_step_run = 0;

    // Initialize min_energy.
    double min_energy;

    // Count the total number of true edges (edges connected to valid nodes).
    generalParams.true_num_edges = 0;
    for (int i = 0; i < coordInfoVecs.num_edges; i++)
    {
        if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX)
        {
            generalParams.true_num_edges += 1;
        }
    }

    // double COMPRESS = 2.0227;
    //  double COMPRESS2 = -2.0227;

    /////////////////////////////////////////////////////////////////
    /////////////////////// MEMBRANE RELATED ////////////////////////
    /////////////////////////////////////////////////////////////////

    // Part 5

    // Membrane relatde parameters and variables initialization.

    // Initilize the following vectors with zeros.
    std::vector<double> nodenormal_1(generalParams.maxNodeCount, 0.0);
    std::vector<double> nodenormal_2(generalParams.maxNodeCount, 0.0);
    std::vector<double> nodenormal_3(generalParams.maxNodeCount, 0.0);

    // Variable to keep track of how many times the linearSpringsInfoVecs.spring_constant_rep1 is reduced.
    int reduce_counter = 0;

    // Set VOLUME_FACTOR to the maximum volume ratio (target volume = VOLUME_FACTOR * initial_volume).
    double VOLUME_FACTOR = MAX_VOLUME_RATIO; // 1.6;//2.25;

    // double tip_depth = 0.5;
    // tip_depth is currently unused.

    // Line tension threshold for the activation of line tension (currently not used)
    double LINE_TENSION_THRESHOLD = -10000.0;
    std::cout << "LINE TENSION THRESHOLD for activation of line tension = " << LINE_TENSION_THRESHOLD << std::endl;

    // Volume threshold for the activation of weakened membrane (currently not used).
    double VOLUME_THRESHOLD = 0.0;
    std::cout << "VOLUME THRESHOLD for activation of weakened membrane = " << VOLUME_THRESHOLD << std::endl;

    // The minimum height of the z-coordinate of the membrane node to be considered in the area of weakened mechanical properties.
    double weakened = 0.0; // 1.90;//6.0; Nav changed it from 1.90 t0 0.0 to have weakened area increased. 8/26/24

    // double tip_base = 6.0;
    // tip_base currently unused.

    // RULES_OF_EXPAN controls how the EXPAN_THRESHOLD is applied:
    //  // 1:= Both trianglular areas must exceed the threshold value.
    //  // 2:= If one trianglular area exceeds the treshold value while the other exceeds the secondary threshold value.
    //  // 3:= If the combined area of the two triangles exceed 2*EXPAN_THRESHOLD.
    //  // 4:= If a selected edges exceed the threshold value, split the triangles associated with the edge.

    for (int i = 0; i < generalParams.maxNodeCount; i++)
    {
        generalParams.centerX += coordInfoVecs.nodeLocX[i];
        generalParams.centerY += coordInfoVecs.nodeLocY[i];
        generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
    }

    generalParams.centerX = generalParams.centerX / generalParams.maxNodeCount;
    generalParams.centerY = generalParams.centerY / generalParams.maxNodeCount;
    generalParams.centerZ = generalParams.centerZ / generalParams.maxNodeCount;

    // Initialization of newcenterX, newcenterY, newcenterZ for recentering of the mesh.
    double displacementX, displacementY, displacementZ;
    double newcenterX, newcenterY, newcenterZ;

    std::vector<int> VectorShuffleForGrowthLoop;
    std::vector<int> VectorShuffleForFilamentLoop;
    std::vector<int> VectorShuffleForEdgeswapLoop;

    // Find the min and max height of the membrane nodes and their indices.
    double min_height = coordInfoVecs.nodeLocZ[0];
    double max_height = -10000.0;
    int max_height_index = -1;
    for (int k = 0; k < generalParams.maxNodeCount; k++)
    {
        if (coordInfoVecs.nodeLocZ[k] >= max_height)
        {
            max_height = coordInfoVecs.nodeLocZ[k];
            max_height_index = k;
        }
    }

    // Max and min height of the membrane nodes, these have to be changed if the mesh used is changed.

    // Set the equilibrium length of an edge of the triangle.
    generalParams.Rmin = 0.001; // 0.75;//0.5; nav changed this once again. Made it larger 11/7/24 //0.0001; nav changed this on 10/10/24 for the double layer code. This value of 0.0001 works for the circular sheet. //0.3012; changed by nav on 6/5/24 for flat code. //0.15; //equilibrium length (Nav) changed by nav again to 0.5 from 0.15. 8/5/24 11/8/24 5 worked for small number of nodes. < Nav

    generalParams.abs_Rmin = generalParams.Rmin; // 0.15;
    std::cout << "abs_Rmin = " << generalParams.abs_Rmin << std::endl;

    // Equilibrium distance between membrane node for volume exclusion.
    //  Initialize the following which represents the equilibrium triangular area.
    areaTriangleInfoVecs.initial_area = 1.0; // 0.039;nav changed this to make it larger 11/7/24 //2835;//0.009808;//0.039;//0.03927344;//0.009817; 11/8/24 25 worked for small number of nodes. < Nav
    std::cout << "equilibrium triangular area = " << areaTriangleInfoVecs.initial_area << std::endl;

    // Set ljInfoVecs parameters (currently all set to 0.0, indicating no interactions).

    // Equilibrium triangular area.
    ljInfoVecs.Rmin_M = 0.0;

    // Equilibrium distance between the nucleus particle and membrane.
    ljInfoVecs.Rcutoff_M = 0.0;

    // Maximal interaction range between the nucleus and membrane.
    ljInfoVecs.Rmin_LJ = 0.0; // 3.0//1.0;

    // Equilibrium distance between nuclei.
    ljInfoVecs.Rcutoff_LJ = 0.0; // 3.0;//1.0;

    // Maximal interaction range between the nuclei.
    ljInfoVecs.epsilon_M_att1 = 0.0; // 6.0;//16.0;
    ljInfoVecs.epsilon_M_att2 = 0.0; // 1.0;//1.0;
    std::cout << "Morse_NM_D_att = " << ljInfoVecs.epsilon_M_att1 << std::endl;
    std::cout << "Morse_NM_a_att = " << ljInfoVecs.epsilon_M_att2 << std::endl;

    // Coefficient for the attractive interaction between nuclei and membrane.
    ljInfoVecs.epsilon_M_rep1 = 0.0; // 12.5;//16.0;
    ljInfoVecs.epsilon_M_rep2 = 0.0; // 0.5;//1.0;
    std::cout << "Morse_NM_D_rep = " << ljInfoVecs.epsilon_M_rep1 << std::endl;
    std::cout << "Morse_NM_a_rep = " << ljInfoVecs.epsilon_M_rep2 << std::endl;

    // Coefficient for the repulsive interaction between nuclei and membrane.
    ljInfoVecs.epsilon_LJ_rep1 = 0.0; // 10.0;//0.5;// 0.06;//7.5;
    ljInfoVecs.epsilon_LJ_rep2 = 0.0; // 0.5;//1.0;//1.0;//1.0;
    std::cout << "Morse_NN_D = " << ljInfoVecs.epsilon_LJ_rep1 << std::endl;
    std::cout << "Morse_NN_a = " << ljInfoVecs.epsilon_LJ_rep2 << std::endl;
    // Coefficient of the interaction between nuclei.

    linearSpringInfoVecs.spring_constant_rep1 = 0.01; // 0.023;
    linearSpringInfoVecs.spring_constant_rep2 = 9.0;  // 5.0;
    std::cout << "Membrane volume exclusion Morse D = " << linearSpringInfoVecs.spring_constant_rep1 << std::endl;
    std::cout << "Membrane volume exclusion Morse a = " << linearSpringInfoVecs.spring_constant_rep2 << std::endl;
    // The coefficient used for non-neighboring membrane node volume exclusion.
    // rep1 is the "D" and rep2 is the "alpha" in the standard form of Morse potential.

    generalParams.volume_spring_constant = 5.0; //(1.0/3.0)*areaTriangleInfoVecs.initial_area*1.0;
    std::cout << "spring constant for surface normal expansion (pressure within the cell) = " << generalParams.volume_spring_constant << std::endl;

    generalParams.line_tension_constant = 5.0; // 250.0; // Value that generated flat sheet is 0.0. 8/14/24
    std::cout << "spring constant for the septin ring (before budding) = " << generalParams.line_tension_constant << std::endl;

    // Equilibrium length of each segment of the septin ring.
    generalParams.length_scale = 0.0; // 1.0*generalParams.Rmin;//nav changed this from 0 to the current value to test the boundary nodes. 03/06/2025 //0.85;//0.1577;//1.0*generalParams.Rmin;// 0.8333; //nav changed this to be 0.0 from 1.0. 8/5/24; Flat sheet generated when septin ring was 0.0. 8/14/24


    // Set spring constants for apical basal and vertical layers. 
   // generalParams.vertical_spring_constant = 250.0;
    
    // Set weakened region scaling factors.
    generalParams.maxSpringScaler_linear = 1.0;
    generalParams.maxSpringScaler_area = 1.0;
    generalParams.maxSpringScaler_bend = 1.0;
    double scale_linear = linearSpringInfoVecs.spring_constant * 1.0; // 0.25;//25.0/2.5;//75.0/15.0; flat sheet generated when multiplied by 1; 8/15/24; Changing it to 0.25 makes it wrinkle up from before. Same with all three below 8/15/24
    double scale_bend = bendingTriangleInfoVecs.spring_constant * 1.0;   // 0.05;//10.0/1.0;//75.0/7.5;  flat sheet generated when multiplied by 1; 8/15/24;
    double scale_area = areaTriangleInfoVecs.spring_constant * 0.25;   // 0.25;//50.0/5.0;//75.0/15.0;  flat sheet generated when multiplied by 1; 8/15/24;
    // nav changed all of the above to their original values to see how it affects budding. 8/26/24
    std::cout << "weakened region linear (before budding) = " << scale_linear << std::endl;
    std::cout << "weakened region bend (before budding) = " << scale_bend << std::endl;
    std::cout << "weakened region area (before budding) = " << scale_area << std::endl;

    // Scaling factor of the weakend mechanical properties.
    linearSpringInfoVecs.spring_constant_weak = scale_linear;
    bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
    areaTriangleInfoVecs.spring_constant_weak = scale_area;

    // Set the following for bending angle equilibrium.
   // bendingTriangleInfoVecs.initial_angle = 0.0087;                 // 0.087165870975460;//0.087249;//0.04335; // This is also an angle we need to change to make the code flat. Nav. Make this into some random drastic value.
    bendingTriangleInfoVecs.initial_angle_raft = 0.087165870975460; // 0.087249;//0.04335; //The raft and coat versions are from some weird legacy thing Kevin was testing - Nav.
    bendingTriangleInfoVecs.initial_angle_coat = 0.087165870975460; // 0.087249;//0.04335;
  //  std::cout << "equilibrium bending angle of the membrane = " << bendingTriangleInfoVecs.initial_angle << std::endl;
    // raft and coat are current unused due to the assumption of uniform preferred curvature.

    bendingTriangleInfoVecs.initial_angle_bud = bendingTriangleInfoVecs.initial_angle[1];
    std::cout << "equilibrium bending angle of the bud = " << bendingTriangleInfoVecs.initial_angle_bud << std::endl;

    // following vectors currently empty.
    std::vector<int> pull_nodes_up;   // = {35,    76,    79,   111,   113,   151,   153,   360,   361,   362,   363,   364,   365,   505,   506,   515,   516,   593,   632};//{35, 360,   361,   362,   363,   364,   365};
    std::vector<int> pull_nodes_down; // = {38,    86,    89,   121,   123,   144,   146,   378,   379,   380,   381,   382,   383,   535,   536,   545,   546,   602,   626};//{38, 378,   379,   380,   381,   382,   383};
    std::vector<int> push_nodes_down;
    std::vector<int> push_nodes_up;

    /////////////////////////////////////////////////////////////////
    ////////////////// END OF MEMBRANE RELATED //////////////////////
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    //////////////////////// NULCEUS RELATED ////////////////////////
    /////////////////////////////////////////////////////////////////

    // Part 6
    // Nucleus related parameters and variables initialization. // Nav - nucleus related stuff should be commented out for flat code. 6/2/24

    // Set beta1 and beta2 to manually push the nucleus tip and the remainder of the nucleus vertically.
    double beta1 = 0.0;
    double beta2 = 0.0;
    std::cout << "manual push speed for the nucleus tip = " << beta1 << std::endl;
    std::cout << "manual push speed for the remainder of the nucleus = " << beta2 << std::endl;
    // beta1 is the vertical speed (0, 0, beta1) applied to the nucleus tip.
    // beta2 is the vertical speed (0, 0, beta2) applied to the remainder of the nucleus.

    // V1, V2 and V3 are vectors representing the (x,y,z)-coordinates of the nucleus particles.
    // Note: These vectors are currently initialized with single values for demonstration purposes.
    std::vector<double> V1 = {-0.0}; /*, 0.0  ,  0.1966  ,  0.5547 ,  -0.4689 ,   0.2422 ,  -0.2229,
                                    -0.4312 ,  -0.0185 ,   0.2887 ,   0.3187 ,   0.7140 ,
                                     0.2231 ,  -0.1921 ,	  -0.5541 ,   -0.1542 ,   -0.1689 ,    0.4391 ,
                                    -0.6661 ,  -0.6381 ,   0.6256 ,   0.0466 ,  -0.0610 ,   0.5134};
                                     */
    std::vector<double> V2 = {0.0};  /*, 0.0 ,  -0.4595 ,  -0.4129 ,   0.0954 ,   0.1764 ,   0.4186 ,
                                    -0.5602 ,  -0.6082 ,  -0.5318 ,   0.3561 ,   0.0753 ,
                                    -0.0917 ,  -0.2596 , 0.2871 ,  -0.3918 ,   0.5195 ,   0.5579 ,
                                    -0.2805 ,   0.0133  , -0.0073 ,   0.7426 ,   0.0614 ,  -0.1506};
                                      */
    std::vector<double> V3 = {0.0};  // initailly 0.6390. changing it to 0.0 (nav)
                                     /*, 0.0 ,  -0.5511 ,   0.0267 ,  -0.5240  , -0.4004 ,   0.2850 ,
                                                                   0.2032 ,  -0.1771 ,   0.4048 ,   0.3461 ,  -0.2034 ,
                                                                   0.5041 ,  -0.4535 ,	-0.1241 ,   0.5722 ,  -0.3748 ,  -0.1335 ,
                                                                   -0.0851 ,   0.3213 ,   0.2389 ,   0.0044 ,  -0.7424 ,  -0.7450};
                                                                   */

    // Push the (x,y,z)-coordinates of the nuclues particles into the ljInfoVecs vectors.
    // These vectors will be used for interactions between the nuclues and other particles.
    for (int i = 0; i < V1.size(); i++)
    {
        ljInfoVecs.LJ_PosX_all.push_back(V1[i]);
        ljInfoVecs.LJ_PosY_all.push_back(V2[i]);
        ljInfoVecs.LJ_PosZ_all.push_back(V3[i]);
    }

    // Set NUCLEUS_UpperHEM_BASE and NUCLEUS_LOWERHEM_BASE, which define the z-coordinate requirement for nucleus particles.
    // to be considered tip-region or base-region. This is used to determine where to apply spring or constant force.
    double NUCLEUS_UPPERHEM_BASE = 0.0; // initially 0.5. Changing it to 0.0 (nav)
    double NUCLEUS_LOWERHEM_BASE = 0.0; // initially -0.6. Changing it to 0.0 (nav)

    //////////////////////////////////////////////////////////////////
    ///////////////// END OF NUCLEUS RELATED /////////////////////////
    //////////////////////////////////////////////////////////////////

    // Part 7

    //////////////////////////////////////////////////////////////////
    /////////// IDENTIFYING REGIONS WITH DIFFERENT MECH PROP /////////
    //////////////////////////////////////////////////////////////////


    /*
    
                    DV BOUNDARY REGION 
    
    */
    
    
    // calculate center of the apical layer of disc
    int max_upperhem_nodes = 0;
    for (int i = 0; i<generalParams.maxNodeCount; i++){
    
      if (generalParams.nodes_in_upperhem[i]==1){
          //std::cout<<"node in upperhem = " << generalParams.nodes_in_upperhem[i] <<std::endl;
          max_upperhem_nodes++;
      }
    
    }
    
    std::cout<< "\n total number of nodes in upperhem = " << max_upperhem_nodes <<std::endl;

    
    double cx_b=0, cy_b=0, cz_b=0;  int nLower=0;
    double cx_a=0, cy_a=0, cz_a=0;  int nUpper=0;
    
    for(int i=0;i<generalParams.maxNodeCount;++i){
        if(generalParams.nodes_in_upperhem[i]==1){
            cx_a += coordInfoVecs.nodeLocX[i];
            cy_a += coordInfoVecs.nodeLocY[i];
            cz_a += coordInfoVecs.nodeLocZ[i];
            ++nUpper;
        }
        else if(generalParams.nodes_in_upperhem[i] == -1){
            cx_b += coordInfoVecs.nodeLocX[i];
            cy_b += coordInfoVecs.nodeLocY[i];
            cz_b += coordInfoVecs.nodeLocZ[i];
            ++nLower;
        }
    }
    cx_a/=nUpper;  cy_a/=nUpper;  cz_a/=nUpper; // Apical center   
    cx_b/=nLower;  cy_b/=nLower;  cz_b/=nLower; // Basal center   
    
    
    // Boundary node along the x-direction -- this is arbitrary. You can choose any axis to go along. 
    
    double R = 0.0, DR = 0.0;
    
    for (int i = 0; i<generalParams.maxNodeCount; i++){
        if (generalParams.nodes_in_upperhem[i]!=1) continue;
        double dx = coordInfoVecs.nodeLocX[i] - cx_a;
        double dy = sqrt((coordInfoVecs.nodeLocY[i] - cy_a)*(coordInfoVecs.nodeLocY[i] - cy_a));
        if (dy<1e-3 && dx>R) R=dx; DR = dy;
    }
    
    
    // for convenience 
    
    generalParams.centerX = cx_a; 
    generalParams.centerY = cy_a;
    generalParams.centerZ = cz_a;
    
    std::cout<< "boundary node in the x direction (apical layer) =  (" <<R<<", "<<DR<< ") \n"<<std::endl;
    
    
    // Identifying regions with different mechanical properties and finding the coundary nodes and edges of the upper hemisphere.

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// Define vectors to store information about boundary edges and edge indices.
    // std::vector<int> out;
    ////int ALPHA;
    // std::vector<bool> boundary_edges;
    // boundary_edges.reserve(coordInfoVecs.num_edges);
    //
    //// Populate the boundary_edges vector to identify boundary edges in the mesh.
    // for (int i = 0; i < coordInfoVecs.num_edges; i++){
    //		if (coordInfoVecs.edges2Triangles_1[i] == coordInfoVecs.edges2Triangles_2[i]){
    //			  boundary_edges.push_back(true); // If the edge connects to only one triangle, it's a boundary edge.
    //		}
    //		else {
    //			  boundary_edges.push_back(false);
    //		}
    // }
    //
    //// Create a vector to store edge indices that are not boundary edges.
    // std::vector<int> edgeIndices;
    // edgeIndices.reserve(coordInfoVecs.num_edges);
    // for (int i = 0; i < coordInfoVecs.num_edges; ++i){
    //		if (boundary_edges[i] == false){
    //			  edgeIndices.push_back(i); // Store the indices of non-boundary edges in the edgeIndices vector.
    //		}
    //		else {
    //			  edgeIndices.push_back(-1);
    //		}
    // }
    //
    //
    //// Remove invalid (negative) indices from edgeIndices.
    // auto it = remove_if(edgeIndices.begin(), edgeIndices.end(),  [](const int i) {return i < 0; });
    // edgeIndices.erase(it, edgeIndices.end());
    //// Make sure boundaries_in_upperhem is resized appropriately.
    // generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
    //
    //
    //// Moved boundary part 8/26/24 nav
    //  //Find the boundary of the nodes_in_upperhem region
    //	//generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
    //	std::vector<int> boundary_node_list;
    //	std::vector<int> boundary_edge_list;
    //	for (int i = 0; i < coordInfoVecs.num_edges; i++){
    //		double T1 = coordInfoVecs.edges2Triangles_1[i];
    //		double T2 = coordInfoVecs.edges2Triangles_2[i];
    //		if (T1 >= (INT_MAX - 1000) || T1 < 0 || T2 >= (INT_MAX-1000) || T2 < 0){
    //			continue;
    //		}
    //		if (generalParams.triangles_in_upperhem[T1] == 1 && generalParams.triangles_in_upperhem[T2] == 1){
    //			generalParams.boundaries_in_upperhem[i] = 1;
    //			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
    //		  generalParams.triangles_in_upperhem[T1] = 0;
    //			generalParams.triangles_in_upperhem[T2] = 0;
    //			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
    //			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
    //			//std::cout<<"septin ring nodes - bdrynode1 = "<<bdry_node1<<std::endl;
    //      boundary_node_list.push_back(bdry_node1);
    //			boundary_node_list.push_back(bdry_node2);
    //			boundary_edge_list.push_back(i);
    //
    //
    //			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
    //			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
    //			//coordInfoVecs.isNodeFixed[bdry_node1] = true;
    //			//coordInfoVecs.isNodeFixed[bdry_node2] = true;
    //		}
    //	/*	else if (generalParams.triangles_in_upperhem[T1] != 1 && generalParams.triangles_in_upperhem[T2] == 1){
    //			generalParams.boundaries_in_upperhem[i] = 1;
    //			std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
    //			generalParams.triangles_in_upperhem[T1] = 0;
    //			generalParams.triangles_in_upperhem[T2] = 0;
    //			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
    //			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
    //			boundary_node_list.push_back(bdry_node1);
    //			boundary_node_list.push_back(bdry_node2);
    //			boundary_edge_list.push_back(i);
    //			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
    //			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
    //		 coordInfoVecs.isNodeFixed[bdry_node1] = true;
    //		 coordInfoVecs.isNodeFixed[bdry_node2] = true;
    //		}*/
    //		else {
    //			generalParams.boundaries_in_upperhem[i] = -1;
    //		  //std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
    //		}
    //	}
    //
    // std::cout<<"size of boundary_node_list (this is double-counted) = "<<boundary_node_list.size()<<std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Ensure boundaries_in_upperhem is resized.
    // generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
    // Ensure boundaries_in_upperhem is resized.
    generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);

    std::cout << "boundaries in upperhem = " << generalParams.boundaries_in_upperhem.size() << std::endl;

    std::vector<int> boundary_edge_list;
    std::vector<int> boundary_node_list;

    std::cout << "edges2Triangles_1 = " << coordInfoVecs.edges2Triangles_1.size() << std::endl;
    std::cout << "edges2Triangles_2 = " << coordInfoVecs.edges2Triangles_2.size() << std::endl;

    //std::cout << "generalParams.nodes in upperhem size = " << generalParams.nodes_in_upperhem.size() << std::endl;

    for (int i = 0; i < coordInfoVecs.num_edges; i++)
    {
        int T1 = static_cast<int>(coordInfoVecs.edges2Triangles_1[i]);
        int T2 = static_cast<int>(coordInfoVecs.edges2Triangles_2[i]);

        // std::cout<<"it got till here - nav "<< std::endl;// it got till here.

        // Optionally check if the triangle indices are valid.
        if (T1 < 0 || T2 < 0 || T1 >= (INT_MAX - 1000) || T2 >= (INT_MAX - 1000))
        {
            continue;
        }

        // std::cout<<"it got till here - nav 1 "<< std::endl; // it got till here.

        // Mark edge as boundary if the two adjacent triangle IDs are identical.
        // std::cout<<"T1 = "<<T1<<std::endl;
        // std::cout<<"T2 = "<<T2<<std::endl;
        // std::cout<<"generalParams.edges_in_upperhem["<<i<<"] = "<< generalParams.edges_in_upperhem[i]<<std::endl;

        if (T1 == T2 && generalParams.edges_in_upperhem[i] == -1)
        {
            generalParams.boundaries_in_lowerhem[i] = 1;
            boundary_edge_list.push_back(i); // This is to store the total number of boundary edges.
            int bdry_node1 = static_cast<int>(coordInfoVecs.edges2Nodes_1[i]);
            int bdry_node2 = static_cast<int>(coordInfoVecs.edges2Nodes_2[i]);
            boundary_node_list.push_back(bdry_node1);
            boundary_node_list.push_back(bdry_node2);

            // mark these nodes as boundary or (fixed).
            generalParams.nodes_in_upperhem[bdry_node1] = -1;
            generalParams.nodes_in_upperhem[bdry_node2] = -1;
            coordInfoVecs.isNodeFixed[bdry_node1] = false;
            coordInfoVecs.isNodeFixed[bdry_node2] = false;
            
        }
        // This is for apical boundary nodes.
        if (T1 == T2 && generalParams.edges_in_upperhem[i] == 1)
        { // nav added the second conditional && generalParams.nodes_in_upperhem[i]==1 so that the new apical model would work.
            generalParams.boundaries_in_upperhem[i] = 1;
            boundary_edge_list.push_back(i);

            // std::cout<<"it got till here - nav 3 "<< std::endl;
            int bdry_node1 = static_cast<int>(coordInfoVecs.edges2Nodes_1[i]);
            int bdry_node2 = static_cast<int>(coordInfoVecs.edges2Nodes_2[i]);
            boundary_node_list.push_back(bdry_node1);
            boundary_node_list.push_back(bdry_node2);

            // Optionally mark these nodes as boundary (or fixed).
            generalParams.nodes_in_upperhem[bdry_node1] = 1; // nav changed this from 0 to 10.
            generalParams.nodes_in_upperhem[bdry_node2] = 1; // nav changed this from 0 to 10.
            coordInfoVecs.isNodeFixed[bdry_node1] = false;
            coordInfoVecs.isNodeFixed[bdry_node2] = false;

            // std::cout<<"it got till here - nav 4 "<< std::endl;
        }
        else
        {
            generalParams.boundaries_in_upperhem[i] = -1;
        }
    }

    std::cout << "size of boundary_edge_list = " << boundary_edge_list.size() << std::endl;
    std::cout << "size of boundary_node_list (double-counted) = " << boundary_node_list.size() << std::endl;

    // Count the true number of edges in the upper hemisphere.
    int true_num_edges_in_upperhem = 0;
    int edges_in_upperhem_COUNT = 0;
    for (int i = 0; i < coordInfoVecs.num_edges; i++)
    {
        if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0)
        {
            true_num_edges_in_upperhem += 1;
            edges_in_upperhem_COUNT += 1;
        }
    }

    // Define a row2 vector to store specific node indices. These are the specific nodes in the septin ring. Nav replaced the following hard coded row with a boundary node list from the segment of code she moved. 8/26/24
    // std::vector<int> row2 = boundary_node_list;//nav commenting out on 03/09/2025. //{35 ,   76 ,   79 ,  111 ,  113 ,  151 ,  153 ,  360 ,  361 ,  362 ,  363 ,  364 ,  365 ,  505 ,  506 ,  515 ,  516 ,  593 ,  632};
    // nav commenting the above out to see if row2 can be defined later for the septin ring. 8/19/24 nav putting these back in.

    // Identify nodes in the upper hemisphere based on their z-coordinates.
    // for (int i = 0; i < generalParams.maxNodeCount; i++){
    //   // generalParams.nodes_in_upperhem[i] = -1; (nav commented this out for flat virus code 5/29/24)
    //  // generalParams.nodes_in_upperhem[i] = 1; // (nav uncommented this for flat virus code 5/29/24)
    //  std::cout<<"nodes in upperhem = "<<generalParams.nodes_in_upperhem[i]<<std::endl;
    //}

    // Nav commented out the following for flat virus code 5/29/24. nav put it back in 8/19/24. nav commented out again 03/09/2025
    // for (int i = 0; i < row2.size(); i++){
    //		generalParams.nodes_in_upperhem[row2[i]] = 1;
    // }

    // nav commented out 03/09/2025. Calculate the minimum z-coordinate of the nodes in row2.
    // double min_septin_z = 1000.0;
    // for (int i = 0; i < row2.size(); i++){
    //		if (coordInfoVecs.nodeLocZ[row2[i]] < min_septin_z){
    //			min_septin_z = coordInfoVecs.nodeLocZ[row2[i]];
    //		}
    //}

    // nav commented out again 03/09/2025. Identify additional nodes in the upper hemisphere based on their z-coordinates.
    // for (int i = 0; i < generalParams.maxNodeCount; i++){
    //		if (coordInfoVecs.nodeLocZ[i] > (min_septin_z)){
    //			generalParams.nodes_in_upperhem[i] = 1;
    //		}
    //}

    //// Identify triangles in the upper hemisaphere based on their nodes. nav uncommented till line 745 for tests 8/19/24. nav commented out 03/09/2025
    for (int i = 0; i < coordInfoVecs.num_triangles; i++)
    {
        if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_1[i] < 0)
        {
            generalParams.triangles_in_upperhem[i] = -1;
            continue;
        }
        else if (coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_2[i] < 0)
        {
            generalParams.triangles_in_upperhem[i] = -1;
            continue;
        }
        else if (coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_3[i] < 0)
        {
            generalParams.triangles_in_upperhem[i] = -1;
            continue;
        }

        int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
        int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
        int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];

        if ((aaa + bbb + ccc) == 3)
        {
            generalParams.triangles_in_upperhem[i] = 1;
        }
        else
        {
            generalParams.triangles_in_upperhem[i] = -1;
        }
    }

    // Identify edges in the upper hemisphere based on their triangles.
    // Store the indices of edges in the upper hemisphere in the edges_in_upperhem_list vector.
    // Count the number of edges in the upper hemisphere.

    // int edges_in_upperhem_COUNT = 0;
    //
    // for (int i = 0; i < coordInfoVecs.num_edges; i++){
    //// NEW: Compute the edge’s midpoint z coordinate
    // double avg_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[i]] +
    //                 coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[i]]) / 2.0;
    // if (avg_z <= 0.5) {
    //     // If the midpoint is not above z=0.5, mark the edge as not in upper hemisphere.
    //     generalParams.edges_in_upperhem[i] = -1;
    //     generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //     continue;
    // }
    //
    // else{
    //     int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];
    //     int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];
    //     if (aaa == 1 && bbb == 1){
    //         generalParams.edges_in_upperhem[i] = 1;
    //         generalParams.edges_in_upperhem_list[i] = i;
    //         edges_in_upperhem_COUNT += 1;
    //     }
    //     else if (aaa == 1 || bbb == 1){
    //         generalParams.edges_in_upperhem[i] = 1;
    //         generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //         edges_in_upperhem_COUNT += 1;
    //     }
    //     else{
    //         generalParams.edges_in_upperhem[i] = -1;
    //         generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //     }
    // }
    //
    // }
    //// nav commented out the following and added the part above. 03/10/2025
    // for (int i = 0; i < coordInfoVecs.num_edges; i++){
    //		if (coordInfoVecs.edges2Triangles_1[i] >= (INT_MAX-1000) || coordInfoVecs.edges2Triangles_1[i] < 0){
    //			generalParams.edges_in_upperhem[i] = -1;
    //			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //			continue;
    //		}
    //		else if (coordInfoVecs.edges2Triangles_2[i] >= (INT_MAX-1000) || coordInfoVecs.edges2Triangles_2[i] < 0){
    //   			generalParams.edges_in_upperhem[i] = -1;
    //   			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //   			continue;
    //		}
    //		else{
    //   			int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];
    //   			int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];
    //   			if (aaa == 1 && bbb == 1){
    //   				generalParams.edges_in_upperhem[i] = 1;
    //   				generalParams.edges_in_upperhem_list[i] = i;
    //   				edges_in_upperhem_COUNT += 1;
    //   			}
    //   			else if (aaa == 1 || bbb == 1){
    //   				generalParams.edges_in_upperhem[i] = 1;
    //   				generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //   				edges_in_upperhem_COUNT += 1;
    //   			}
    //   			else{
    //   				generalParams.edges_in_upperhem[i] = -1;
    //   				generalParams.edges_in_upperhem_list[i] = -INT_MAX;
    //   			}
    //        // std::cout<< "Edges in upperhem = "<<generalParams.edges_in_upperhem[i]<<std::endl;
    //		}
    // }

    std::cout << "INITIAL EDGES IN UPPERHEM = " << edges_in_upperhem_COUNT << std::endl;

    int COUNTING_EDGE = 0;
    for (int y = 0; y < coordInfoVecs.num_edges; y++)
    {
        if (generalParams.edges_in_upperhem_list[y] >= 0)
        {
            COUNTING_EDGE += 1;
        }
        generalParams.edges_in_upperhem_list_length = COUNTING_EDGE;
    }

    /*

    //Find the boundary of the nodes_in_upperhem region
    // Define vectors to store the indices of boundary nodes and edges of the upper hemisphere.
    std::vector<int> boundary_node_list;
    std::vector<int> boundary_edge_list;

    // Find the boundary nodes and edges of the upper hemisphere.
    for (int i = 0; i < coordInfoVecs.num_edges; i++){
            double T1 = coordInfoVecs.edges2Triangles_1[i];
            double T2 = coordInfoVecs.edges2Triangles_2[i];
            if (T1 >= (INT_MAX - 1000) || T1 < 0 || T2 >= (INT_MAX-1000) || T2 < 0){
                continue;
            }
       // Have to change this to reflect what the boundary condition is.
       // This condition is checking to if you take a sphere and circle out a region then every triangle in the boundary of that region will have the condition that one triangle is within the region and the other is not.

       //The following is a nav added part 6/4/24:

       double node1 = coordInfoVecs.edges2Nodes_1[i];
        double node2 = coordInfoVecs.edges2Nodes_2[i];

        // Check if the two nodes associated with the edge are the same
        if (node1 == node2) {
            generalParams.boundaries_in_upperhem[i] = 1;

            boundary_node_list.push_back(node1);
            boundary_node_list.push_back(node2);
            boundary_edge_list.push_back(i);
        } else {
            generalParams.boundaries_in_upperhem[i] = -1;
        }
    }*/

    /*		if (generalParams.triangles_in_upperhem[T1] == generalParams.triangles_in_upperhem[T2]){ //1 && generalParams.triangles_in_upperhem[T2] ==1){ //!= 1){ this is what it was before. Nav is changing it for the flat code.
                generalParams.boundaries_in_upperhem[i] = 1;

          double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
                double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];

          boundary_node_list.push_back(bdry_node1);
                boundary_node_list.push_back(bdry_node2);
                boundary_edge_list.push_back(i);
            }
            /*else if (generalParams.triangles_in_upperhem[T1] ==1 /*!= 1 This is what it was before nav changed it  && /*generalParams.triangles_in_upperhem[T2] == 1){
                generalParams.boundaries_in_upperhem[i] = 1;
                double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
                double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
                boundary_node_list.push_back(bdry_node1);
                boundary_node_list.push_back(bdry_node2);
                boundary_edge_list.push_back(i);

            }*/
    // else {
    // generalParams.boundaries_in_upperhem[i] = -1;
    //	}
    //}*///Nav commented the above out to test her own version 6/4/24
    // This is where the boundary part labelled (moved boundary part) was originally. 8/26/24
    // Nav commented this out to restore to original version to make changes once more. Let's see! 8/5/24

    // If conditions on 815 and 830 in the original code need to be changed to reflect boundary condition ie T1==T2 - nav

    // Initialize the generalParams.edge_to_ljparticle vector to store the connection between an edge and LJ particle (nucleus particle).
    // for (int i = 0; i < coordInfoVecs.num_edges; i++){
    //		generalParams.edge_to_ljparticle.push_back(-1);
    //};

    /////////////////////////////////////////////////////////////////////
    ////////////// END OF IDENTIFYING REG. WITH DIFF. MECH PROP /////////
    /////////////////////////////////////////////////////////////////////

    // Part 9

    // Compute the initial volume of the system.
    ComputeVolume(
        generalParams,
        coordInfoVecs,
        linearSpringInfoVecs,
        ljInfoVecs);
    double initial_volume;

    std::cout << "Initial volume = " << initial_volume << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// START OF ACTUAL SIMULATION /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    /* Build the initial gradient weakend scale */

    // Initialize variables for gradient weakening.
    dtb = 0.0;                         // dtb := distance to boundary
    generalParams.septin_ring_z = 0.0; // was 0.0, nav changed it to test 8/5/24
    generalParams.boundary_z = 0.0;

    // Loop through all boundary nodes to calculate the distance to the cell tip node.
    for (int k = 0; k < boundary_node_list.size(); k++)
    {
        double n1 = boundary_node_list[k];
        double dist_x = coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1];
        double dist_y = coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1];
        double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1];
        double temp_dist = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
        generalParams.septin_ring_z += coordInfoVecs.nodeLocZ[n1];
        if (temp_dist >= dtb)
        {
            dtb = temp_dist;
            /* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
        }
    }
    std::cout << "dtb = " << dtb << std::endl;

    // dtb will be only calculated once so we can effectively keep the Hill eqn curve consistent with only horizontal shift
    dtb_max = dtb + (generalParams.Rmin); // Calculate dtb_max, which is dtb plus the equilibrium length of an edge (Rmin).

    std::cout << "initial distance between cell tip and the boundary of weakened area = " << dtb << std::endl;
    std::cout << "Notice that here, the distance from the tip to the boundary is slightly extended by half of the equilibrium length of an edge" << std::endl;

    // Calculate the hill equation constant K using dtb and dtb_max, and set the hill (equation) coefficient.
    generalParams.hilleqnconst = (dtb * dtb_scaler) / dtb_max;
    generalParams.hilleqnpow = targetHillEqnPow;

    if (generalParams.SCALE_TYPE == 4)
    {
        std::cout << "hill equation constant K = " << generalParams.hilleqnconst << std::endl;
        std::cout << "hill (equation) coefficient = " << generalParams.hilleqnpow << std::endl;
    }
    // NOTE: IN THIS SIMULATION, THE LOCATION WHERE 50% WEAKENING IS EXPERIENCED IS LOCATED SLIGHTLY AWAY FROM THE SEPTIN RING. THIS IS DUE TO THE FACT THAT IN ISOTROPIC CASE, SEPTIN RING LOCATION MUST BE SUFFICIENTLY WEAKENED TO INDUCE BUDDING.

    // Transfer host vectors to device memory and perform gradient weakening update.
    utilities_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
    utilities_ptr->gradient_weakening_update_host_vecs(sigma,
                                                       coordInfoVecs.nodeLocX[max_height_index],
                                                       coordInfoVecs.nodeLocY[max_height_index],
                                                       coordInfoVecs.nodeLocZ[max_height_index],
                                                       dtb,
                                                       dtb_max,
                                                       generalParams,
                                                       coordInfoVecs,
                                                       build_ptr->hostSetInfoVecs);

    // Calculate the boundary elements for each node.
    for (int u = 0; u < generalParams.maxNodeCount; u++)
    {
        int BETA = utilities_ptr->nodes2Triangles_host_vecs(
            u,
            build_ptr->hostSetInfoVecs,
            coordInfoVecs,
            generalParams,
            auxVecs);
    }

    // Transfer updated host vectors back to device memory.
    utilities_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// DSP EVERSION SIMULATION LOOP ///////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "STARTING THE EVERSION ACTUAL SIMULATION" << std::endl;

    // Print VTK file for visualization.

    bool runSim_Eversion = true;

    // EVERSION SIMULATION LOOP
    while (runSim_Eversion == true)
    {
        double current_time = 0.0;
        int translate_counter = 0;

        // Simulate until the specified relaxation time is reached.
        while (current_time < relax_max_steps_before_strain_tensor_and_gradient_descent * (Max_Runtime))
        {

            translate_counter += 1;

            // Solve forces and update positions of nodes
            Solve_Forces(); // check if the force modules take in individual rest lengths or just a single one (Rmin) Also remove Rmin from the code.
            AdvancePositions(coordInfoVecs,
                             generalParams,
                             domainParams);

            // Calculate the new energy of the system.
            new_total_energy = linearSpringInfoVecs.linear_spring_energy +
                               areaTriangleInfoVecs.area_triangle_energy +
                               bendingTriangleInfoVecs.bending_triangle_energy;

            old_total_energy = new_total_energy;
            current_time += generalParams.dt;

        } // end of time loop

        // Print simulation results for "steady state" initial condition before strain application and gradient relaxation.
        std::cout << "Time used for 'steady state' initial condition before strain application and gradient relaxation = " << current_time << std::endl;
        std::cout << "Current total energy (before strain application and gradient relaxation) = " << new_total_energy << std::endl;
        std::cout << "LINEAR ENERGY = " << linearSpringInfoVecs.linear_spring_energy << std::endl;
        std::cout << "BEND ENERGY = " << bendingTriangleInfoVecs.bending_triangle_energy << std::endl;
        std::cout << "AREA ENERGY = " << areaTriangleInfoVecs.area_triangle_energy << std::endl;
        std::cout << "True current total volume (before strain tensor and gradient relaxation) = " << generalParams.true_current_total_volume << std::endl;
        std::cout << "Current KBT = " << generalParams.kT << std::endl;

        if (std::isnan(new_total_energy))
        {
            std::cout << "Total energy is NaN. Exit code (-1) in System.cu main function." << std::endl;
            exit(-1);
        }

        // Calculate area of apical and basal layers and print them out.

        double current_apical_area = 0.0;
        // double current_basal_area = 0.0; basal layer can be done after apical layer.

        // check total number of triangles and then calculate area.

        std::cout << "INT_MAX = " << INT_MAX << std::endl;

        for (int k = 0; k < coordInfoVecs.num_triangles; k++)
        {

            if (coordInfoVecs.triangles2Nodes_1[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_1[k] <= (-INT_MAX + 1000.0) ||
                coordInfoVecs.triangles2Nodes_2[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_2[k] <= (-INT_MAX + 1000.0) ||
                coordInfoVecs.triangles2Nodes_3[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_3[k] <= (-INT_MAX + 1000.0))
            {
                continue;
            }
            else
            {
                if (generalParams.triangles_in_upperhem[k] == 1)
                {
                    double r1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
                    double r1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
                    double r1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
                    double r2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
                    double r2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
                    double r2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
                    double r3x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
                    double r3y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
                    double r3z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
                    double norm_r1r2 = sqrt((r2x - r1x) * (r2x - r1x) + (r2y - r1y) * (r2y - r1y) + (r2z - r1z) * (r2z - r1z));
                    double norm_r2r3 = sqrt((r3x - r2x) * (r3x - r2x) + (r3y - r2y) * (r3y - r2y) + (r3z - r2z) * (r3z - r2z));
                    double norm_r3r1 = sqrt((r3x - r1x) * (r3x - r1x) + (r3y - r1y) * (r3y - r1y) + (r3z - r1z) * (r3z - r1z));
                    //std::cout << "norm_r1r2 = " << norm_r1r2 << std::endl;
                    //std::cout << "norm_r2r3 = " << norm_r2r3 << std::endl;
                    //std::cout << "norm_r3r1 = " << norm_r3r1 << std::endl;
                    double s = (norm_r1r2 + norm_r2r3 + norm_r3r1) / 2.0;
                    double area = sqrt(s * (s - norm_r1r2) * (s - norm_r2r3) * (s - norm_r3r1));
                    current_apical_area += area;
                }
            }
        }

        // Below is your initial apical surface area.
        double Initial_Apical_Area = current_apical_area;
        std::cout << "Initial apical surface area (before growth and edge swaps) = " << Initial_Apical_Area << std::endl;

        // Spring constants for volume interactions. At the moment I dont think these are needed so I'm commenting it out. We shall see if they're used in other computations in different places.
        // generalParams.volume_spring_constant = 0.2;
        // std::cout<<"Spring constant for surface normal expansion (volume spring constant) = "<<generalParams.volume_spring_constant<<std::endl;

        if (generalParams.SCALE_TYPE == 4)
        {
            // Set scaling factors for different mechanical properties in the hill function. It's 1 if you want the top and bottom to have the same mechanical properties. Otherwise change them. Originally
            // used for the bud vs mother cell.
            generalParams.maxSpringScaler_linear = 1.0;
            generalParams.maxSpringScaler_area = 1.0;
            generalParams.maxSpringScaler_bend = 1.0;
        }

        std::cout << "maxSpringScaler_linear (not 1.0 if we want max linear spring in the hill function scaling not equal to mother cell) = " << generalParams.maxSpringScaler_linear << std::endl;
        std::cout << "maxSpringScaler_area (not 1.0 if we want max area spring in the hill function scaling not equal to mother cell) = " << generalParams.maxSpringScaler_area << std::endl;
        std::cout << "maxSpringScaler_bend (not 1.0 if we want max bend spring in the hill function scaling not equal to mother cell) = " << generalParams.maxSpringScaler_bend << std::endl;

        double scale_linear = linearSpringInfoVecs.spring_constant * 0.75;
        double scale_bend = bendingTriangleInfoVecs.spring_constant * 0.135;
        double scale_area = areaTriangleInfoVecs.spring_constant * 0.75;

        std::cout << "weakened region linear = " << scale_linear << std::endl;
        std::cout << "weakened region bend = " << scale_bend << std::endl;
        std::cout << "weakened region area = " << scale_area << std::endl;

        // Update the weakened spring constants.

        linearSpringInfoVecs.spring_constant_weak = scale_linear;
        bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
        areaTriangleInfoVecs.spring_constant_weak = scale_area;

        // Initial volume of the system.
        initial_volume = generalParams.true_current_total_volume;
        generalParams.eq_total_volume = generalParams.true_current_total_volume * VOLUME_FACTOR; // This is for setting different equilibrium volume to mimic growth or shirnkage.
        std::cout << "true current total volume = " << generalParams.true_current_total_volume << std::endl;
        std::cout << "eq total volume = " << generalParams.eq_total_volume << std::endl;

        // Print VTK file for visualization.
        // storage->print_VTK_File();

        // Weakened start for simulating perhaps some chemical properties? Undetermined for now.
        int gradient_descent_iteration = 2; // this number needs to be changed.
        bool WEAKENED_STRAT = false;
        bool GRADIENT_DESCENT_ALGORITHM_TRIGGERED;
        bool NeedToRebuildDifStructAfterGradientDescent = false;
        int number_of_simulation_step = 0;

        // Main simulation loop.
        while (initial_kT > 0)
        {
            if (gradient_descent_iteration >= NKBT)
            {
                runSim_Eversion = false;
                initial_kT = -1;
                gradient_descent_iteration+=1;
                break;
            }

            // Now relax mesh (initial relaxation before strain tensor application)

            current_time = 0.0;
            translate_counter = 0;
            double VOLUME_RATIO = generalParams.true_current_total_volume / generalParams.eq_total_volume; // this is for volume measurements.

            GRADIENT_DESCENT_ALGORITHM_TRIGGERED = true;
            bool end_of_relaxation = false;

            // The following loop is set to calclulate the total number of simulation steps.
            while (current_time < Max_Runtime)
            {
                number_of_simulation_step += 1;
                if (Max_Runtime <= 0.0)
                {
                    std::cout << "Max_Runtime is set to be 0 or negative! " << std::endl;
                    break;
                }

                Solve_Forces();

                AdvancePositions(
                    coordInfoVecs,
                    generalParams,
                    domainParams);

                // Calculate the new total energy of the system.
                new_total_energy = linearSpringInfoVecs.linear_spring_energy +
                                   areaTriangleInfoVecs.area_triangle_energy +
                                   bendingTriangleInfoVecs.bending_triangle_energy; // +0.5*energy_rep;

                energy_gradient = sqrt((new_total_energy - old_total_energy) * (new_total_energy - old_total_energy)) / old_total_energy;
                old_total_energy = new_total_energy;
                current_time += generalParams.dt;

// nav commented out below because we are no longer recentering the mesh. 
//                // Computing center of just the apical layer.
//                if (translate_counter % translate_frequency == 0)
//                {
//
//                    // Compute the new center of the system.
//                    newcenterX = 0.0;
//                    newcenterY = 0.0;
//                    newcenterZ = 0.0;
//
//                    for (int i = 0; i < generalParams.maxNodeCount; i++)
//                    {
//                        if (generalParams.nodes_in_upperhem[i] == 1)
//                        {
//                            newcenterX += coordInfoVecs.nodeLocX[i];
//                            newcenterY += coordInfoVecs.nodeLocY[i];
//                            newcenterZ += coordInfoVecs.nodeLocZ[i];
//                        }
//                    }
//                    newcenterX = newcenterX / generalParams.maxNodeCount;
//                    newcenterY = newcenterY / generalParams.maxNodeCount;
//                    newcenterZ = newcenterZ / generalParams.maxNodeCount;
//
//                    // Compute the displacement vector.
//                    displacementX = newcenterX - generalParams.centerX;
//                    displacementY = newcenterY - generalParams.centerY;
//                    displacementZ = newcenterZ - generalParams.centerZ;
//
//                    // Update the positions of all nodes and LJ particles.
//                    for (int i = 0; i < generalParams.maxNodeCount; i++)
//                    {
//                        coordInfoVecs.nodeLocX[i] += -displacementX;
//                        coordInfoVecs.nodeLocY[i] += -displacementY;
//                        coordInfoVecs.nodeLocZ[i] += -displacementZ;
//                    }
//                    // you can use this part to compute the positions of the LJ_particles too if you choose to include that energy.
//
//                    // Recompute the volume of the system after the translation
//                    ComputeVolume(
//                        generalParams,
//                        coordInfoVecs,
//                        linearSpringInfoVecs,
//                        ljInfoVecs);
//                }

            } // end of main simulation loop eversion

            // Print VTK file for visualization.
            storage->print_VTK_File();

            end_of_relaxation = true;
            double current_center_x = 0.0;
            double current_center_y = 0.0;
            double bdry_to_tip_height = 0.0; // this here is now for the total height of the apical layer of the tissue.

            if (generalParams.SCALE_TYPE == 4)
            {
                max_height = -10000.0;
                for (int k = 0; k < generalParams.maxNodeCount; k++)
                {
                    if (generalParams.nodes_in_upperhem[k] == 1)
                    {
                        current_center_x += coordInfoVecs.nodeLocX[k];
                        current_center_y += coordInfoVecs.nodeLocX[k];
                    }

                    if (coordInfoVecs.nodeLocZ[k] >= max_height)
                    {
                        max_height = coordInfoVecs.nodeLocZ[k];
                        max_height_index = k;
                    }
                }
                current_center_x = current_center_x / generalParams.maxNodeCount;
                current_center_y = current_center_y / generalParams.maxNodeCount;

//                if (generalParams.nonuniform_wall_weakening_bend == false && generalParams.nonuniform_wall_weakening_linear == false && generalParams.nonuniform_wall_weakening_area == false)
//                {
//                    bdry_to_tip_height = 0.0;
//
//                    for (int y = 0; y < boundary_edge_list.size(); y++)
//                    {
//                        double edge_mdpt_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[boundary_edge_list[y]]] +
//                                              coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[boundary_edge_list[y]]]) /
//                                             2.0;
//                        bdry_to_tip_height += sqrt(pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z, 2.0));
//                    }
//                    bdry_to_tip_height = bdry_to_tip_height / boundary_edge_list.size();
//
//                    for (int y = 0; y < coordInfoVecs.num_edges; y++)
//                    {
//                        if (generalParams.edges_in_upperhem_list[y] >= 0 &&
//                            generalParams.edges_in_upperhem_list[y] != INT_MAX &&
//                            generalParams.edges_in_upperhem_list[y] <= (INT_MAX - 1000) &&
//                            generalParams.edges_in_upperhem_list[y] >= (-INT_MAX + 1000) &&
//                            generalParams.boundaries_in_upperhem[y] != 1)
//                        {
//                            if (coordInfoVecs.edges2Nodes_1[y] < 0 || coordInfoVecs.edges2Nodes_1[y] >= (INT_MAX - 1000))
//                            {
//                                continue;
//                            }
//                            else if (coordInfoVecs.edges2Nodes_2[y] < 0 || coordInfoVecs.edges2Nodes_2[y] >= (INT_MAX - 1000))
//                            {
//                                continue;
//                            }
//                            double edge_mdpt_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[y]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[y]]) / 2.0;
//                            double current_edge_to_tip_height = sqrt(pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z, 2.0));
//                            // generalParams.Rmin = linearSpringInfoVecs.edge_initial_length[y];
//
//                            if (bdry_to_tip_height >= (generalParams.Rmin * generalParams.ratio_for_HillFunctionStiffness))
//                            {
//                                if (generalParams.nonuniform_wall_weakening_bend == false && generalParams.nonuniform_wall_weakening_area == false && generalParams.nonuniform_wall_weakening_linear == false && display_token == true)
//                                {
//                                    std::cout << "generalParams.nonuniform_wall_weakening_XXXX is set to be true from this point" << std::endl;
//                                    display_token = false;
//                                }
//                                generalParams.nonuniform_wall_weakening_bend = true;
//                                generalParams.nonuniform_wall_weakening_area = true;
//                                generalParams.nonuniform_wall_weakening_linear = true;
//                            }
//                            else if (bdry_to_tip_height < (generalParams.Rmin * generalParams.ratio_for_HillFunctionStiffness))
//                            {
//                                generalParams.nonuniform_wall_weakening_bend = false;
//                                generalParams.nonuniform_wall_weakening_linear = false;
//                                generalParams.nonuniform_wall_weakening_area = false;
//                            }
//                        }
//                    }
//                }
            }

            //// Here is where you would insert chemical diffusion. We dont have that yet so no need to worry about it.

            // Check if relaxation loop has ended.
            if (end_of_relaxation == true)
            {

                // Compute volume of system
                ComputeVolume(generalParams,
                              coordInfoVecs,
                              linearSpringInfoVecs,
                              ljInfoVecs);

                // Check if mesh overcompression is happening for volume related termination.
                if ((generalParams.true_current_total_volume / initial_volume) < 0.6 || (generalParams.true_current_total_volume / initial_volume) >= MAX_VOLUME_RATIO)
                {

                    // Update true_num_edges and store data before terminating the simulation
                    generalParams.true_num_edges = 0;
                    for (int i = 0; i < coordInfoVecs.num_edges; i++)
                    {
                        if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX)
                        {
                            generalParams.true_num_edges += 1;
                        }
                    }

                    storage->print_VTK_File();

                    // Print appropriate message based on the termination reason.
                    if (generalParams.true_current_total_volume / initial_volume < 0.6)
                    {
                        std::cout << "Cell over compression 60%" << std::endl;
                    }
                    else if (generalParams.true_current_total_volume / initial_volume >= MAX_VOLUME_RATIO)
                    {
                        std::cout << "Target volume ratio exceeded. Current volume ratio = " << generalParams.true_current_total_volume / initial_volume << std::endl;
                    }

                    // Print relevant simulation statistics
                    std::cout << "Current number of gradient descent iteration performed at volume-related termination = " << gradient_descent_iteration << std::endl;
                    std::cout << "Current number of simulation step at volume-related termination = " << number_of_simulation_step << std::endl;

                    // Termination simulation. Nav commented out for testing without volume thresholds 8/26/24
                    Max_Runtime = 0.0;
                    runSim_Eversion = false;
                    initial_kT = -1;
                    break;
                }
                else {// This is just added in for the sake of it. Remove it later - nav 05-13-25 the yeast budding model chose a volume related termination, you should have a curvature related termination
                    Max_Runtime = 0.0;
                    runSim_Eversion = false;
                    initial_kT = -1;
                    break;
                }
            }

            // Calculate current apical surface area
            double current_apical_area = 0.0;

            for (int k = 0; k < coordInfoVecs.num_triangles; k++)
            {

                // Check if triangle data is valid
                if (coordInfoVecs.triangles2Nodes_1[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_1[k] <= (-INT_MAX + 1000.0) ||
                    coordInfoVecs.triangles2Nodes_2[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_2[k] <= (-INT_MAX + 1000.0) ||
                    coordInfoVecs.triangles2Nodes_3[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_3[k] <= (-INT_MAX + 1000.0))
                {
                    continue;
                }
                else
                {
                    // Check if triangle is in the upper hemisphere
                    if (generalParams.triangles_in_upperhem[k] == 1)
                    {

                        // Calculate triangle area and update current_bud_area
                        double r1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
                        double r1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
                        double r1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
                        double r2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
                        double r2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
                        double r2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
                        double r3x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
                        double r3y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
                        double r3z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
                        double norm_r1r2 = sqrt((r2x - r1x) * (r2x - r1x) + (r2y - r1y) * (r2y - r1y) + (r2z - r1z) * (r2z - r1z));
                        double norm_r2r3 = sqrt((r3x - r2x) * (r3x - r2x) + (r3y - r2y) * (r3y - r2y) + (r3z - r2z) * (r3z - r2z));
                        double norm_r3r1 = sqrt((r3x - r1x) * (r3x - r1x) + (r3y - r1y) * (r3y - r1y) + (r3z - r1z) * (r3z - r1z));
                        double s = (norm_r1r2 + norm_r2r3 + norm_r3r1) / 2.0;
                        double area = sqrt(s * (s - norm_r1r2) * (s - norm_r2r3) * (s - norm_r3r1));
                        current_apical_area += area;
                    }
                }
            }

            // Check if bud surface area ratio conditions are met for termination
            if (current_apical_area / Initial_Apical_Area >= MAX_APICAL_AREA_RATIO)
            {

                // Print relevant message and statistics
                std::cout << "Target apical surface area ratio exceeded. Current apical surface area ratio = " << current_apical_area / Initial_Apical_Area << std::endl;
                std::cout << "Current number of relaxation iterations performed at area-related termination = " << gradient_descent_iteration << std::endl;
                std::cout << "Current number of simulation step at area-related termination = " << number_of_simulation_step << std::endl;

                // Update true_num_edges and store data before terminating the simulation
                generalParams.true_num_edges = 0;
                for (int i = 0; i < coordInfoVecs.num_edges; i++)
                {
                    if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX)
                    {
                        generalParams.true_num_edges += 1;
                    }
                }
                storage->print_VTK_File();

                // Terminate the simulation
                Max_Runtime = 0.0;
                runSim = false;
                initial_kT = -1;
                break;
            }
        }
    
    // here recenter the mesh first.


// nav commented out because we are no longer recentering the mesh. 
//    // Compute the new center of the system.
    newcenterX = 0.0;
    newcenterY = 0.0;
    newcenterZ = 0.0;

    for (int i = 0; i < generalParams.maxNodeCount; i++)
    {
        if (generalParams.nodes_in_upperhem[i] == 1)
        {
            newcenterX += coordInfoVecs.nodeLocX[i];
            newcenterY += coordInfoVecs.nodeLocY[i];
            newcenterZ += coordInfoVecs.nodeLocZ[i];
        }
    }
    newcenterX = newcenterX / generalParams.maxNodeCount;
    newcenterY = newcenterY / generalParams.maxNodeCount;
    newcenterZ = newcenterZ / generalParams.maxNodeCount;
    
    generalParams.centerX = newcenterX; //sumX / nNodes;
        generalParams.centerY =newcenterY;// sumY / nNodes;
        generalParams.centerZ = newcenterZ;

//    // Compute the displacement vector.
//    displacementX = newcenterX - generalParams.centerX;
//    displacementY = newcenterY - generalParams.centerY;
//    displacementZ = newcenterZ - generalParams.centerZ;
//
//    // Update the positions of all nodes and LJ particles.
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//        coordInfoVecs.nodeLocX[i] += -displacementX;
//        coordInfoVecs.nodeLocY[i] += -displacementY;
//        coordInfoVecs.nodeLocZ[i] += -displacementZ;
//    }
    storage->print_VTK_File();

    //// commenting out the following for some tests.  Nav 05/04/2025
    std::cout << "number of simulation steps = " << generalParams.Tf << std::endl;


    
        double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
        int nNodes = coordInfoVecs.nodeLocX.size();
        //int layerflag = 0; // this is the layer that will be excluded. 
        for (int i = 0; i < nNodes; i++)
        {  
            //if (layerflag == 1 && i >= nNodes/2) break;
            //if (layerflag == -1 && i < nNodes/2) continue;
                sumX += coordInfoVecs.nodeLocX[i];
                sumY += coordInfoVecs.nodeLocY[i];
                sumZ += coordInfoVecs.nodeLocZ[i];
        }
        generalParams.centerX = sumX / nNodes;
        generalParams.centerY = sumY / nNodes;
        generalParams.centerZ = sumZ / nNodes;
        
        //std::cout<<"center at "<<layerflag<<" = ("<<generalParams.centerX<<", "<<generalParams.centerY<<", "<<generalParams.centerZ<<") "<< std::endl;

    
    std::cout<<"Tf = "<< generalParams.Tf<<std::endl;
//    
//    for (int i = 0; i<coordInfoVecs.num_edges; i++){
//        if (i >= 10) break;
//        std::cout<< "initial rest length before strain tensor, at edge = "<< i << " = " <<linearSpringInfoVecs.edge_rest_length[i]<<std::endl;
//        
//    }

    //double tol = 1e-8;
    double iter =100;// 1/(generalParams.dt*generalParams.tol);
    
    // if you leave ^ this to be 100 it'll take wayy longer for your simulations. Make it at least 10000.
    
    ///////////////////////////////////////////////////////////////////
    // Okay so let's start the simulation loop for the strain tensor //
    ///////////////////////////////////////////////////////////////////
    
    int stages = generalParams.Tf;
    //generalParams.tol = 1e-4;
    
 //   int layerflag = 0;
    
    LambdaField lambda;
    
//    for (int stage = 0; stage<=stages; stage++){
    
        double frac = 1.0;///generalParams.Tf; // this fraction should remain what it is. Dont change it. 
        
      sumX = 0.0, sumY = 0.0, sumZ = 0.0;
        
        
         nNodes = coordInfoVecs.nodeLocX.size();
        int layerflag = 0; // this is the layer that will be excluded. 
        for (int i = 0; i < nNodes; i++)
        {  
            //if (layerflag == 1 && i >= nNodes/2) break;
            //if (layerflag == -1 && i < nNodes/2) continue;
                sumX += coordInfoVecs.nodeLocX[i];
                sumY += coordInfoVecs.nodeLocY[i];
                sumZ += coordInfoVecs.nodeLocZ[i];
        }
        generalParams.centerX = sumX / nNodes;
        generalParams.centerY = sumY / nNodes;
        generalParams.centerZ = sumZ / nNodes;
    
        StrainTensorGPU::buildVertexLambda(generalParams, coordInfoVecs, lambda, frac);
        
          //std::cout << "\n[Stage " << stage << "] spring lengths BEFORE strain:\n";
//    for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
//        if (generalParams.edges_in_upperhem[e] == 0) continue;  // skip vertical
//        double L0 = linearSpringInfoVecs.edge_initial_length[e];
//        std::cout
//            << "  edge[" << e << "]: initial_length = "
//            << L0 << "\n";
//    }
        StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs,generalParams, lambda, linearSpringInfoVecs, layerflag);
        
        //Solve_Forces();
        //AdvancePositions(coordInfoVecs, generalParams, domainParams);
        
        
        // === PRINT AFTER ===
   // std::cout << "[Stage " << stage << "] spring lengths AFTER strain:\n";
//    for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
//        if (generalParams.edges_in_upperhem[e] == 0) continue;  // skip vertical
//        double Lstar = linearSpringInfoVecs.edge_final_length[e];
//        std::cout
//            << "  edge[" << e << "]: strained_length = "
//            << Lstar << "\n";
//    }
generalParams.tol = 1e-4;
              
              ///////////////////////////////////////////////////////////////////////  
        // resize edge_rest_length to match the size of the steps between edge_initial and edge_final 
        // loop for changing the rest lengths between stages. Here the spring length is broken down into steps and passed into the force kernels. 
        for (int step = 0; step<= iter; step++){
        
            for(int i = 0; i<coordInfoVecs.num_edges; i++){ // this can be done from inside the strain tensor file. 
                
                double dl = (linearSpringInfoVecs.edge_final_length[i] - linearSpringInfoVecs.edge_initial_length[i])/iter;
                linearSpringInfoVecs.edge_rest_length[i] += dl;
            }
            
           // double avg_mvt = 999999.0;
            // generalParams.dx = 0.0;
            
           // std::cout<< "Starting relaxation loop."<<std::endl;
            //int k = 0;
            
            int k = relaxUntilConverged(*this);
            
            new_total_energy = linearSpringInfoVecs.linear_spring_energy;
                                   //areaTriangleInfoVecs.area_triangle_energy +
                                   //bendingTriangleInfoVecs.bending_triangle_energy; 
                                   
            std::cout<<"Energy = "<<new_total_energy<<", iterations = "<<k<<", = Average Movement = "<<generalParams.dx<<std::endl;
                        
           // std::cout<<"Exiting relaxation loop. Average movement = "<<generalParams.dx<<", and k = "<< k<<" at step = "<<step<<std::endl;//" and stage = "<<stage<<std::endl;     
            if (step > 0 && step % 2 == 0) {
                
                storage->print_VTK_File(); // this is just for now. To get better simulation results we're gonna view them more frequently later. 
            }

          }
//////////////////////////////////////
//            
           // do{
             //   generalParams.dx = 0.0; 
                
//                double x1[coordInfoVecs.nodeLocX.size()];
//                double y1[coordInfoVecs.nodeLocX.size()];
//                double z1[coordInfoVecs.nodeLocX.size()];
//                for (int n = 0; n<coordInfoVecs.nodeLocX.size(); n++){ // this can be handled from inside the nodeAdvance.cu
//                    x1[n] = coordInfoVecs.nodeLocX[n];
//                    y1[n] = coordInfoVecs.nodeLocY[n];
//                    z1[n] = coordInfoVecs.nodeLocZ[n];
//                } 
                
                //avg_mvt = 0.0;
               // Solve_Forces();
                // AdvancePositions(coordInfoVecs, generalParams, domainParams);
                
                
//                double x2[coordInfoVecs.nodeLocX.size()];
//                double y2[coordInfoVecs.nodeLocX.size()];
//                double z2[coordInfoVecs.nodeLocX.size()];
//                
//                for (int n = 0; n<coordInfoVecs.nodeLocX.size(); n++){ // same as above. 
//                    x2[n] = coordInfoVecs.nodeLocX[n];
//                    y2[n] = coordInfoVecs.nodeLocY[n];
//                    z2[n] = coordInfoVecs.nodeLocZ[n];
//                    
//                    generalParams.dx+= sqrt((x1[n]-x2[n])*(x1[n]-x2[n])+(y1[n]-y2[n])*(y1[n]-y2[n])+(z1[n]-z2[n])*(z1[n]-z2[n]));
//                }
//                
//                
                //avg_mvt=generalParams.dx;   
              //  k++;
                
               // if (k == 4500) break;
            //}while (generalParams.dx>generalParams.tol);
            
            
            //if (step == 100 || step == 200 || step == 300 || step == 400 || step == 500 || step == 600 || step == 700 || step == 800 || step == 900 || step == 1000 || step == 2000 || step == 5000 || step == 7000 || step == 9000 ){ // badly done ik but I'll fix it later. 
            
      //  }
        
        storage->print_VTK_File();
        
  //  }
    
    // strain tensor acts on edge_inital to give edge_final 

}
};




// Function to assign the shared pointer to storage.
void System::assignStorage(std::shared_ptr<Storage> _storage)
{
    storage = _storage;
};

// Function to set the weak pointer to the SystemBuilder.
void System::set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr)
{
    weak_bld_ptr = _weak_bld_ptr;
};

// Function to initialize memory for thrust vectors and set coordInfoVecs values from input.
void System::initializeSystem(HostSetInfoVecs & hostSetInfoVecs)
{
    std::cout << "Initializing" << std::endl;

    // Set the max node count, edge count and triangle count.
    generalParams.maxNodeCount = hostSetInfoVecs.nodeLocX.size();
    coordInfoVecs.num_edges = hostSetInfoVecs.edges2Nodes_1.size();
    coordInfoVecs.num_triangles = hostSetInfoVecs.triangles2Nodes_1.size();

    std::cout << "num nodes: " << generalParams.maxNodeCount << std::endl;
    std::cout << "num edges: " << coordInfoVecs.num_edges << std::endl;
    std::cout << "num elems: " << coordInfoVecs.num_triangles << std::endl;
    // Allocate memory for various vectors using preallocated memory size.
    int mem_prealloc = 4;

    // Resize and initialize the following various vectors.
    coordInfoVecs.scaling_per_edge.resize(mem_prealloc * coordInfoVecs.num_edges, 0.0);
    hostSetInfoVecs.scaling_per_edge.resize(coordInfoVecs.scaling_per_edge.size(), 0.0);

    coordInfoVecs.soln_per_triangle.resize(mem_prealloc * coordInfoVecs.num_triangles, INT_MAX);
    coordInfoVecs.b_per_triangle.resize(mem_prealloc * coordInfoVecs.num_triangles, INT_MAX);

    coordInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), false);
    coordInfoVecs.prevNodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    coordInfoVecs.prevNodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.prevNodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    coordInfoVecs.nodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.nodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
    coordInfoVecs.nodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    // coordInfoVecs.nodeVelX.resize(mem_prealloc*hostSetInfoVecs.nodeVelX.size(), 0.0);
    // coordInfoVecs.nodeVelY.resize(mem_prealloc*hostSetInfoVecs.nodeVelY.size(), 0.0);
    // coordInfoVecs.nodeVelZ.resize(mem_prealloc*hostSetInfoVecs.nodeVelZ.size(), 0.0);

    coordInfoVecs.nodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
    coordInfoVecs.nodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
    coordInfoVecs.nodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);

    coordInfoVecs.triangles2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Nodes_3.resize(mem_prealloc * coordInfoVecs.num_triangles);

    coordInfoVecs.triangles2Edges_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Edges_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
    coordInfoVecs.triangles2Edges_3.resize(mem_prealloc * coordInfoVecs.num_triangles);

    coordInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    coordInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    coordInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);

    hostSetInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    hostSetInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
    hostSetInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);

    coordInfoVecs.edges2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_edges);
    coordInfoVecs.edges2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_edges);

    coordInfoVecs.edges2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_edges);
    coordInfoVecs.edges2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_edges);

    coordInfoVecs.nndata1.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata2.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata3.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata4.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata5.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata6.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata7.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata8.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.nndata9.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.SurfaceNormalX.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.SurfaceNormalY.resize(mem_prealloc * generalParams.maxNodeCount);
    coordInfoVecs.SurfaceNormalZ.resize(mem_prealloc * generalParams.maxNodeCount);

    generalParams.nodes_in_upperhem.resize(mem_prealloc * generalParams.maxNodeCount);
    generalParams.triangles_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_triangles);
    generalParams.edges_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges);
    generalParams.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
    generalParams.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
    generalParams.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);

    hostSetInfoVecs.nodes_in_upperhem.resize(generalParams.nodes_in_upperhem.size());
    generalParams.nodes_in_upperhem = hostSetInfoVecs.nodes_in_upperhem;
    hostSetInfoVecs.triangles_in_upperhem.resize(generalParams.triangles_in_upperhem.size());
    hostSetInfoVecs.edges_in_upperhem.resize(generalParams.edges_in_upperhem.size());
    generalParams.edges_in_upperhem = hostSetInfoVecs.edges_in_upperhem;
    hostSetInfoVecs.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
    hostSetInfoVecs.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
    hostSetInfoVecs.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);

    hostSetInfoVecs.nodes2Triangles_1.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_2.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_3.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_4.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_5.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_6.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_7.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_8.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    hostSetInfoVecs.nodes2Triangles_9.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);

    // Part 17

    // Resize vectors to allocate memory for nodes-to-triangles mapping.
    coordInfoVecs.nodes2Triangles_1.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_2.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_3.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_4.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_5.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_6.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_7.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_8.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
    coordInfoVecs.nodes2Triangles_9.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);

    // Copy nodes-to-triangles mapping information from hostSetInfoVecs to coodInfoVecs and others.
    thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(), coordInfoVecs.nodes2Triangles_1.end(), hostSetInfoVecs.nodes2Triangles_1.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(), coordInfoVecs.nodes2Triangles_2.end(), hostSetInfoVecs.nodes2Triangles_2.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(), coordInfoVecs.nodes2Triangles_3.end(), hostSetInfoVecs.nodes2Triangles_3.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(), coordInfoVecs.nodes2Triangles_4.end(), hostSetInfoVecs.nodes2Triangles_4.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(), coordInfoVecs.nodes2Triangles_5.end(), hostSetInfoVecs.nodes2Triangles_5.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(), coordInfoVecs.nodes2Triangles_6.end(), hostSetInfoVecs.nodes2Triangles_6.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(), coordInfoVecs.nodes2Triangles_7.end(), hostSetInfoVecs.nodes2Triangles_7.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(), coordInfoVecs.nodes2Triangles_8.end(), hostSetInfoVecs.nodes2Triangles_8.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(), coordInfoVecs.nodes2Triangles_9.end(), hostSetInfoVecs.nodes2Triangles_9.begin());

    // copy info to GPU
    std::cout << "Copying" << std::endl;
    thrust::copy(hostSetInfoVecs.isNodeFixed.begin(), hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());

    // Print information about fixed nodes in hostSetInfoVecs and coordInfoVecs.
    std::cout << "fixed_node_in_host: " << std::endl;
    for (int k = 0; k < hostSetInfoVecs.isNodeFixed.size(); k++)
    {
    }
    std::cout << "end_of_fixed_node_host_printout" << std::endl;
    std::cout << "fixed_node_in_device: " << std::endl;
    for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++)
    {
    }
    std::cout << "end_of_fixed_node_device_printout" << std::endl;
    std::cout << "size of host fixed " << hostSetInfoVecs.isNodeFixed.size() << std::endl;
    std::cout << "size of device fixed " << coordInfoVecs.isNodeFixed.size() << std::endl;

    // initialize various vectors with zeros or values from hostSetInfoVecs.
    //  Fill operations for other nodeForce and prevNodeForce vectors.
    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

    thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
    thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
    thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);

    // Copy node locations and other related information from hostSetInfoVecs to coordInfoVecs and other copy operations for triangles, edges and other related vectors.
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.prevNodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.prevNodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.nodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(), hostSetInfoVecs.triangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(), hostSetInfoVecs.triangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(), hostSetInfoVecs.triangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin());

    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(), hostSetInfoVecs.triangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(), hostSetInfoVecs.triangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(), hostSetInfoVecs.triangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin());

    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(), hostSetInfoVecs.edges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(), hostSetInfoVecs.edges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin());

    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(), hostSetInfoVecs.edges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(), hostSetInfoVecs.edges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin());

    thrust::copy(hostSetInfoVecs.nndata1.begin(), hostSetInfoVecs.nndata1.end(), coordInfoVecs.nndata1.begin());
    thrust::copy(hostSetInfoVecs.nndata2.begin(), hostSetInfoVecs.nndata2.end(), coordInfoVecs.nndata2.begin());
    thrust::copy(hostSetInfoVecs.nndata3.begin(), hostSetInfoVecs.nndata3.end(), coordInfoVecs.nndata3.begin());
    thrust::copy(hostSetInfoVecs.nndata4.begin(), hostSetInfoVecs.nndata4.end(), coordInfoVecs.nndata4.begin());
    thrust::copy(hostSetInfoVecs.nndata5.begin(), hostSetInfoVecs.nndata5.end(), coordInfoVecs.nndata5.begin());
    thrust::copy(hostSetInfoVecs.nndata6.begin(), hostSetInfoVecs.nndata6.end(), coordInfoVecs.nndata6.begin());
    thrust::copy(hostSetInfoVecs.nndata7.begin(), hostSetInfoVecs.nndata7.end(), coordInfoVecs.nndata7.begin());
    thrust::copy(hostSetInfoVecs.nndata8.begin(), hostSetInfoVecs.nndata8.end(), coordInfoVecs.nndata8.begin());
    thrust::copy(hostSetInfoVecs.nndata9.begin(), hostSetInfoVecs.nndata9.end(), coordInfoVecs.nndata9.begin());

    // Resize and initialize the 'u' vector.
    coordInfoVecs.u.resize(mem_prealloc * coordInfoVecs.num_triangles);

    // Part 18

    // Allocate memory for additiional data structures.

    // Area triangle info vec.
    // Number of area springs is the number of triangles
    std::cout << "Mem" << std::endl;
    // Allocate memory for temporary node information in unreduced form for area springs
    areaTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

    // Allocate memory for temporary node information in reduced form for area springs.
    areaTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
    areaTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

    // beinding triangle info vec
    // num bending springs is the number of times each edge is between two triangles.
    bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;

    // Allocate memory for temporary node information in unreduced form for bending springs.
    bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

    // Allocate memory for temporary node information in reduced form for bending springs.
    bendingTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    bendingTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
    
    bendingTriangleInfoVecs.initial_angle.resize(coordInfoVecs.num_triangles);
    thrust::fill(bendingTriangleInfoVecs.initial_angle.begin(), bendingTriangleInfoVecs.initial_angle.end(), 0.0087);
        
    // linear springs info vectors.
    //  Allocate memory for temporary node information in unreduced form for linear springs.
    linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);

    // Allocate memory for temporary node information in reduced form for bending springs.
    linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
    linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);

    // Clear edge_initial_length vector for linear springs.
    // linearSpringInfoVecs.edge_initial_length.clear();
    // linearSpringInfoVecs.edge_rest_length.clear();

    // linearSpringInfoVecs.edge_rest_length.resize(hostSetInfoVecs.edge_rest_length.size());
    linearSpringInfoVecs.edge_final_length.resize(hostSetInfoVecs.edge_initial_length.size());
    linearSpringInfoVecs.edge_initial_length = hostSetInfoVecs.edge_initial_length;
    linearSpringInfoVecs.edge_final_length = hostSetInfoVecs.edge_initial_length;
    
    
    std::cout << "host edge_initial_length size = " << hostSetInfoVecs.edge_initial_length.size() << std::endl;
    std::cout << "device edge_initial_length size = " << linearSpringInfoVecs.edge_initial_length.size() << std::endl;

    //  for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
    //    int i = coordInfoVecs.edges2Nodes_1[e];
    //    int j = coordInfoVecs.edges2Nodes_2[e];
    //    double dx = hostSetInfoVecs.nodeLocX[j] - hostSetInfoVecs.nodeLocX[i];
    //    double dy = hostSetInfoVecs.nodeLocY[j] - hostSetInfoVecs.nodeLocY[i];
    //    double dz = hostSetInfoVecs.nodeLocZ[j] - hostSetInfoVecs.nodeLocZ[i];
    //    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    //    //hostSetInfoVecs.edge_initial_length.push_back(dist);    // already done for initial
    //    hostSetInfoVecs.edge_rest_length.push_back(dist);
    //    //std::cout<< "edge_rest_length = " << hostSetInfoVecs.edge_initial_length[i]<<std::endl; in the current data structure it gave me 1021 edges. That's good. Now that they have been initialized I should start changing the rest lengths.
    //    }

    // thrust::copy(linearSpringInfoVecs.edge_rest_length.begin(),
    //              linearSpringInfoVecs.edge_rest_length.end(),
    //              hostSetInfoVecs.edge_rest_length.begin());

    linearSpringInfoVecs.edge_rest_length = hostSetInfoVecs.edge_initial_length;

    //linearSpringInfoVecs.edge_rest_length.resize(1/(generalParams.dt*generalParams.tol)) //= hostSetInfoVecs.edge_rest_length;

    std::cout << "host edge_rest_length size = " << hostSetInfoVecs.edge_rest_length.size() << std::endl;
    std::cout << "device edge_rest_length size = " << linearSpringInfoVecs.edge_rest_length.size() << std::endl;

    //linearSpringInfoVecs.edge_final_length.resize(coordInfoVecs.num_edges);

    //linearSpringInfoVecs.edge_final_length = linearSpringInfoVecs.edge_initial_length;

    // This loop is to test out whether the edge_rest_length is being initialized properly.
    // for (int i = 0; i < hostSetInfoVecs.edge_rest_length.size(); i++) {
    //    std::cout<< "edge_rest_length # "<< i <<" = "<< hostSetInfoVecs.edge_rest_length[i]<<std::endl;
    //    std::cout<< "edge_initial_length # "<< i <<" = "<< hostSetInfoVecs.edge_initial_length[i]<<std::endl;
    //}

    // Resize the hostSetInfoVecs for data transfer between host and device.
    hostSetInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());

    hostSetInfoVecs.nodeLocX.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeLocY.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeLocZ.resize(coordInfoVecs.nodeLocX.size());
    std::cout << "Host_nodeLocX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;

    // hostSetInfoVecs.nodeVelX.resize(coordInfoVecs.nodeVelX.size());
    // hostSetInfoVecs.nodeVelY.resize(coordInfoVecs.nodeVelY.size());
    // hostSetInfoVecs.nodeVelZ.resize(coordInfoVecs.nodeVelZ.size());

    hostSetInfoVecs.nodeForceX.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeForceY.resize(coordInfoVecs.nodeLocX.size());
    hostSetInfoVecs.nodeForceZ.resize(coordInfoVecs.nodeLocX.size());
    std::cout << "Host_nodeForceX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;

    hostSetInfoVecs.triangles2Nodes_1.resize(coordInfoVecs.triangles2Nodes_1.size());
    hostSetInfoVecs.triangles2Nodes_2.resize(coordInfoVecs.triangles2Nodes_2.size());
    hostSetInfoVecs.triangles2Nodes_3.resize(coordInfoVecs.triangles2Nodes_3.size());
    std::cout << "Host_triangles2Nodes size = " << hostSetInfoVecs.triangles2Nodes_1.size() << std::endl;

    hostSetInfoVecs.triangles2Edges_1.resize(coordInfoVecs.triangles2Edges_1.size());
    hostSetInfoVecs.triangles2Edges_2.resize(coordInfoVecs.triangles2Edges_2.size());
    hostSetInfoVecs.triangles2Edges_3.resize(coordInfoVecs.triangles2Edges_3.size());
    std::cout << "Host_triangles2Edges size = " << hostSetInfoVecs.triangles2Edges_1.size() << std::endl;

    hostSetInfoVecs.edges2Nodes_1.resize(coordInfoVecs.edges2Nodes_1.size());
    hostSetInfoVecs.edges2Nodes_2.resize(coordInfoVecs.edges2Nodes_2.size());
    std::cout << "Host_edges2Nodes size = " << hostSetInfoVecs.edges2Nodes_1.size() << std::endl;

    hostSetInfoVecs.edges2Triangles_1.resize(coordInfoVecs.edges2Triangles_1.size());
    hostSetInfoVecs.edges2Triangles_2.resize(coordInfoVecs.edges2Triangles_2.size());
    std::cout << "Host_edges2Triangles size = " << hostSetInfoVecs.edges2Triangles_1.size() << std::endl;

    hostSetInfoVecs.nndata1.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata2.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata3.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata4.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata5.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata6.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata7.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata8.resize(mem_prealloc * generalParams.maxNodeCount);
    hostSetInfoVecs.nndata9.resize(mem_prealloc * generalParams.maxNodeCount);

    // Print message indicating the system is ready.
    std::cout << "System Ready" << std::endl;

    // Allocate memory for buckets.
    auxVecs.id_bucket.resize(generalParams.maxNodeCount);
    auxVecs.id_value.resize(generalParams.maxNodeCount);
    auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
    auxVecs.id_value_expanded.resize(27 * (generalParams.maxNodeCount));
};

//};
//void System::relaxUntilConverged(int k){
//    
//    generalParams.dx = 0.0;
//    k = 0;
//    while (generalParams.dx>generalParams.tol){        
//        Solve_Forces();
//        AdvancePositions(coordInfoVecs, generalParams, domainParams);
//        
//        k++;
//    }
//    //return k;
////    
//};//

//        /* unary op ?F?2 ------------------------------------------------------------- */
//        struct ForceNorm
//        {
//            __host__ __device__ double operator()(const thrust::tuple<double, double, double> &t) const
//            {
//                double fx = thrust::get<0>(t);
//                double fy = thrust::get<1>(t);
//                double fz = thrust::get<2>(t);
//                return sqrt(fx * fx + fy * fy + fz * fz);
//            }
//        };
//        /* binary op = max ----------------------------------------------------------- */
//        struct MaxOp
//        {
//            __host__ __device__ double operator()(double a, double b) const { return a > b ? a : b; }
//        };
//
//        /* ========================================================================== *
//         *  relaxUntilConverged – returns true on success                             *
//         * ========================================================================== */
//        bool System::relaxUntilConverged(double tol /*=1e-6*/, int maxIter /*=3000*/)
//        {
//            for (int iter = 0; iter < maxIter; ++iter)
//            {
//                /* --------- assemble forces on GPU ---------------------------------- */
//                Solve_Forces(); // <— every physics module writes into nodeForceX/Y/Z
//
//                /* --------- compute max |F| using transform_reduce ------------------ */
//                auto first = thrust::make_zip_iterator(
//                    thrust::make_tuple(coordInfoVecs.nodeForceX.begin(),
//                                       coordInfoVecs.nodeForceY.begin(),
//                                       coordInfoVecs.nodeForceZ.begin()));
//                auto last = thrust::make_zip_iterator(
//                    thrust::make_tuple(coordInfoVecs.nodeForceX.end(),
//                                       coordInfoVecs.nodeForceY.end(),
//                                       coordInfoVecs.nodeForceZ.end()));
//
//                double maxF = thrust::transform_reduce(/* O(N) GPU reduction   */
//                                                       first, last,
//                                                       ForceNorm(), /* unary op  |F|        */
//                                                       0.0,         /* init value           */
//                                                       MaxOp());    /* binary op  max       */
//                /* ------------------------------------------------------------------- */
//                AdvancePositions(coordInfoVecs, generalParams, domainParams);
//
//                if (maxF < tol)
//                { // lambda converged
//                    return true;
//                }
//            }
//            return false; // lambda hit iteration limit
//        };

// turn off line 1821 - 1831 closes towards the very end.







//#include <stdio.h>
//#include "System.h"
//#include "SystemStructures.h"
//#include "AreaTriangles.h"
//#include "LinearSprings.h"
//#include "NodeAdvance.h"
//#include "Storage.h"
//#include "Utilities.h"
//#include "SystemBuilder.h"
//#include <vector>
//#include "VolumeComp.h"
//#include "VolumeSprings.h"
//#include <bits/stdc++.h>
//#include "LineTensionSprings.h"
//#include <math.h>
//#include <list>
//#include "LinearSpringsEnergy.h"
//#include "StrainTensor.h"
//#include <thrust/iterator/zip_iterator.h>
//#include <thrust/functional.h>
//#include <thrust/transform_reduce.h>
//#include <thrust/tuple.h>
//#include "gradientRelax.h"
//
//// somehow the gradient is not being set in my version - Kevin
//
//// Helper function to count elements greater than or equal to zero in a vector.
//int count_bigger(const std::vector<int> &elems)
//{
//    return std::count_if(elems.begin(), elems.end(), [](int c)
//                         { return c >= 0; });
//}
//
//// Constructor for the System class.
//System::System() {};
//
//// Print net force on nodes along a radial line (? ˜ 0) from disc center to boundary
//void System::PrintForce() {
//    // Copy device forces to host
//    thrust::host_vector<double> h_fx = coordInfoVecs.nodeForceX;
//    thrust::host_vector<double> h_fy = coordInfoVecs.nodeForceY;
//    thrust::host_vector<double> h_fz = coordInfoVecs.nodeForceZ;
//    thrust::host_vector<double> h_x  = coordInfoVecs.nodeLocX;
//    thrust::host_vector<double> h_y  = coordInfoVecs.nodeLocY;
//
//    const double desiredTheta = 0.0;    // along +x axis
//    const double eps = 0.01;            // angular tolerance (rad)
//    std::vector<std::pair<double,int>> picks;
//    int N = static_cast<int>(h_x.size());
//    for (int i = 0; i < N; ++i) {
//        double r = std::hypot(h_x[i], h_y[i]);
//        double theta = std::atan2(h_y[i], h_x[i]);
//        double diff = std::fabs(theta - desiredTheta);
//        if (diff > M_PI) diff = 2*M_PI - diff;
//        if (diff < eps) picks.emplace_back(r, i);
//    }
//    std::sort(picks.begin(), picks.end());
//
//    std::printf("   r      Fx       Fy       Fz\n");
//    for (auto &pr : picks) {
//        int idx = pr.second;
//        std::printf("%6.3f  %7.3e  %7.3e  %7.3e\n",
//                    pr.first,
//                    h_fx[idx], h_fy[idx], h_fz[idx]);
//    }
//}
//
//// Function to solve the forces in the system.
//void System::Solve_Forces()
//{
//
//    // Reset all forces to zero.
//    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
//    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
//    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
//
//    // Compute forces and energy due to linear springs.
//    ComputeLinearSprings(
//        generalParams,
//        coordInfoVecs,
//        linearSpringInfoVecs,
//        ljInfoVecs);
//
//    // Compute forces and energy due to area springs. Nav commented out to test Active shape programming mesh type 02/27/2025  . Put back in 03/23/25
////      	ComputeAreaTriangleSprings(
////      		generalParams,
////      		coordInfoVecs,
////      		areaTriangleInfoVecs);
//
//    // Compute forces and energy due to turgor pressure springs. (nav - commenting these out for now for flat surface 5/29/24) nav reintroducing the turgor pressure because the eversion wing does have turgor pressure. 8/17/2024
//    // ComputeTurgorSprings(
//    // generalParams,
//    // coordInfoVecs,
//    // areaTriangleInfoVecs
//    //);
//
//    // Compute forces and energy due to bending springs. Turn this off 10/10/24
////      	ComputeCosTriangleSprings(
////      		generalParams,
////      		coordInfoVecs,
////      		bendingTriangleInfoVecs);
//
//    // Compute forces and energy due to membrane repulsion springs.// Nav commented out to test Active shape programming mesh type 02/27/2025. PUt back in 03/23/25
////      	ComputeMemRepulsionSprings_local(
////      		coordInfoVecs,
////      		linearSpringInfoVecs,
////      		capsidInfoVecs,
////      		generalParams,
////      		auxVecs);
//
//    // Compute forces and energy due to volume springs. //(nav - commenting these out for now for flat surface 5/29/24) Nav had uncommented but she's bringing the comment back because testing out Active shape mesh 02/27/25
////      	ComputeVolume(
////      		generalParams,
////      		coordInfoVecs,
////      		linearSpringInfoVecs,
////      		ljInfoVecs);
//    
//    // Now print forces along the radial line
//   // PrintForce();
//};
//
//// Function to solve the entire system.
//void System::solveSystem()
//{
//
//    // Nav - I dont want to remove these variables. These may come in handy.
//    // coordInfoVecs.k_0 = 20.0;
//    // coordInfoVecs.k_1 = 25.0;
//    // coordInfoVecs.k_2 = 5.0;
//    // coordInfoVecs.k_3 = 5.0;
//    // coordInfoVecs.k_4 = 1.0;
//    // coordInfoVecs.k_ss = 12;//10.75;
//    // coordInfoVecs.beta = 1.0/1.0;///1.45;
//    // coordInfoVecs.gamma = 1.0;
//    // coordInfoVecs.q1 = 10.0;
//    // coordInfoVecs.h = 10.0;
//    
//    generalParams.dt = 0.01;
//    
//    std::cout<<" Simulating wL3 - 0hAPF "<<std::endl;
//    std::cout<<"lambda_iso_center_outDV = -0.12406004 - "<<generalParams.lambda_iso_center_outDV<<std::endl;
//
//    uint mem_prealloc = 4; // Make sure that this number is the same as set in System::initializeSystem found near the bottom of this script. - Kevin. Q. why is this the case? Why does it need to be the same? - Nav.
//   
//   
//    // Create a shared pointer for Utilities.
//    auto utilities_ptr = std::make_shared<Utilities>(coordInfoVecs, generalParams);
//
//    // Create a shared pointer for SystemBuilder.
//    auto build_ptr = weak_bld_ptr.lock();
//    std::cout << "Declaration of rbc and n_rbc complete." << std::endl;
//    std::cout << "Utilities_ptr declaration complete." << std::endl;
//
//
//   
//    // Set the node mass for the simulation.
//    generalParams.nodeMass = 1.0;
//
//
//    // Set kT_growth value for growth events.
//    generalParams.kT_growth = 1.0;
//
//    // Set the SCALE_TYPE for weakening during growth events. (0 to 4)
//    // 0:= Gaussian-like weakening
//    // 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening
//    // 2:= pure Gaussian weakening
//    // 3:= isotropic
//    // 4:= hill equation
//    // Note that (3) is used in combination with sigma = INT_MAX;
//    generalParams.SCALE_TYPE = 3; // Original scale type was 3. Nav changed it to 0 for flat code. 6/2/24. \\ 0 did not work so nav changed it back to 3 8/18/24
//    std::cout << "SCALE TYPE = " << generalParams.SCALE_TYPE << std::endl;
//    std::cout << "0:= sigmoidal Gaussian-like weakening, 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening, 2:= pure Gaussian weakening, 3:= isotropic, 4:= hill equation" << std::endl;
//
//    // Check and set additional parameters based on SCALE_TYPE.
//    if (generalParams.SCALE_TYPE == 1)
//    {
//        generalParams.scaling_pow = 2.0;
//        std::cout << "scaling_pow (this is for SCALE_TYPE = 1 case) = " << generalParams.scaling_pow << std::endl;
//    }
//    if (generalParams.SCALE_TYPE == 0)
//    {
//        generalParams.gausssigma = 0.1;
//        std::cout << "gausssigma (this is for the SCALE_TYPE = 0 case) = " << generalParams.gausssigma << std::endl;
//    }
//
//    // Declare variables for Hill function_dependent wall stiffness.
//    double dtb_scaler, targetHillEqnPow;
//    if (generalParams.SCALE_TYPE == 4)
//    {
//        generalParams.ratio_for_HillFunctionStiffness = 4.0;
//        std::cout << "Hill function dependent wall stiffness triggers when the the distance between tip of the bud and the septin ring is " << generalParams.ratio_for_HillFunctionStiffness << std::endl;
//        std::cout << "times larger than the equilibrium length Rmin" << std::endl;
//        dtb_scaler = 1.0;
//        targetHillEqnPow = 16.0;
//        std::cout << "The EC50 position is scaled by " << dtb_scaler << " on the distance from tip to boundary, hence the EC50 occurs on dtb*" << dtb_scaler << "/dtb_max" << std::endl;
//        std::cout << "Target hill equation power = " << targetHillEqnPow << std::endl;
//    }
//
//
//    if (generalParams.SCALE_TYPE == 0)
//    {
//       double sigma = 0.0;            // INT_MAX; //if this is set to be INT_MAX then we assume isotropic weakening.
//       double sigma_true = sqrt(0.5); // This is the variance used to calculate the scaling of the wall weakening.
//        std::cout << "initial sigma (for gradient distribution variance), based on initial distribution of Cdc42, if using true gaussian weakening = " << sigma << std::endl;
//        std::cout << "If sigma = INT_MAX, then we have isotropic weakening scenario" << std::endl;
//        std::cout << "true sigma (for gaussian-related distribution variance) = " << sigma_true << std::endl;
//    }
//
//    //////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////
//    ////////////////////////// PARAMETER SETTINGS ////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////
//
//    // Part 4
//
//    // Setting various simulation parameters and initializing variables.
//
//    // Declare and initialize variables for growth-related calculations.
//    double initial_kT;
//    initial_kT = generalParams.kT; // Stores the initial kT value for the acceptance of changes after looping through every edge within proximity.
//
//    // Initialize min_energy.
//    double min_energy;
//
//    /////////////////////////////////////////////////////////////////
//    /////////////////////// MEMBRANE RELATED ////////////////////////
//    /////////////////////////////////////////////////////////////////
//
//
//    // Initilize the following vectors with zeros.
//    std::vector<double> nodenormal_1(generalParams.maxNodeCount, 0.0);
//    std::vector<double> nodenormal_2(generalParams.maxNodeCount, 0.0);
//    std::vector<double> nodenormal_3(generalParams.maxNodeCount, 0.0);
//
//
//    // Line tension threshold for the activation of line tension (currently not used)
//    double LINE_TENSION_THRESHOLD = -10000.0;
//    std::cout << "LINE TENSION THRESHOLD for activation of line tension = " << LINE_TENSION_THRESHOLD << std::endl;
//
//    // Volume threshold for the activation of weakened membrane (currently not used).
//    double VOLUME_THRESHOLD = 0.0;
//    std::cout << "VOLUME THRESHOLD for activation of weakened membrane = " << VOLUME_THRESHOLD << std::endl;
//
//    // RULES_OF_EXPAN controls how the EXPAN_THRESHOLD is applied:
//    //  // 1:= Both trianglular areas must exceed the threshold value.
//    //  // 2:= If one trianglular area exceeds the treshold value while the other exceeds the secondary threshold value.
//    //  // 3:= If the combined area of the two triangles exceed 2*EXPAN_THRESHOLD.
//    //  // 4:= If a selected edges exceed the threshold value, split the triangles associated with the edge.
//
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//        generalParams.centerX += coordInfoVecs.nodeLocX[i];
//        generalParams.centerY += coordInfoVecs.nodeLocY[i];
//        generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
//    }
//
//    generalParams.centerX = generalParams.centerX / generalParams.maxNodeCount;
//    generalParams.centerY = generalParams.centerY / generalParams.maxNodeCount;
//    generalParams.centerZ = generalParams.centerZ / generalParams.maxNodeCount;
//
//    // Initialization of newcenterX, newcenterY, newcenterZ for recentering of the mesh.
//    double displacementX, displacementY, displacementZ;
//    double newcenterX, newcenterY, newcenterZ;
//
//    std::vector<int> VectorShuffleForGrowthLoop;
//    std::vector<int> VectorShuffleForFilamentLoop;
//    std::vector<int> VectorShuffleForEdgeswapLoop;
//
//    // Find the min and max height of the membrane nodes and their indices.
//    double min_height = coordInfoVecs.nodeLocZ[0];
//    double max_height = -10000.0;
//    int max_height_index = -1;
//    for (int k = 0; k < generalParams.maxNodeCount; k++)
//    {
//        if (coordInfoVecs.nodeLocZ[k] >= max_height)
//        {
//            max_height = coordInfoVecs.nodeLocZ[k];
//            max_height_index = k;
//        }
//    }
//
//
//    // Equilibrium distance between membrane node for volume exclusion.
//    //  Initialize the following which represents the equilibrium triangular area.
//    areaTriangleInfoVecs.initial_area = 1.0; // 0.039;nav changed this to make it larger 11/7/24 //2835;//0.009808;//0.039;//0.03927344;//0.009817; 11/8/24 25 worked for small number of nodes. < Nav
//    std::cout << "equilibrium triangular area = " << areaTriangleInfoVecs.initial_area << std::endl;
//
//    generalParams.volume_spring_constant = 0.2; //(1.0/3.0)*areaTriangleInfoVecs.initial_area*1.0;
//    std::cout << "spring constant for surface normal expansion (pressure within the cell) = " << generalParams.volume_spring_constant << std::endl;
//
//    generalParams.line_tension_constant = 0.0; // 250.0; // Value that generated flat sheet is 0.0. 8/14/24
//    std::cout << "spring constant for the septin ring (before budding) = " << generalParams.line_tension_constant << std::endl;
//
//    // Equilibrium length of each segment of the septin ring.
//    generalParams.length_scale = 0.0; // 1.0*generalParams.Rmin;//nav changed this from 0 to the current value to test the boundary nodes. 03/06/2025 //0.85;//0.1577;//1.0*generalParams.Rmin;// 0.8333; //nav changed this to be 0.0 from 1.0. 8/5/24; Flat sheet generated when septin ring was 0.0. 8/14/24
//
//
//    // Set weakened region scaling factors.
//    generalParams.maxSpringScaler_linear = 1.0;
//    generalParams.maxSpringScaler_area = 1.0;
//    generalParams.maxSpringScaler_bend = 1.0;
//    double scale_linear = linearSpringInfoVecs.spring_constant * 0.25; // 0.25;//25.0/2.5;//75.0/15.0; flat sheet generated when multiplied by 1; 8/15/24; Changing it to 0.25 makes it wrinkle up from before. Same with all three below 8/15/24
//    double scale_bend = bendingTriangleInfoVecs.spring_constant * 1;   // 0.05;//10.0/1.0;//75.0/7.5;  flat sheet generated when multiplied by 1; 8/15/24;
//    double scale_area = areaTriangleInfoVecs.spring_constant * 0.25;   // 0.25;//50.0/5.0;//75.0/15.0;  flat sheet generated when multiplied by 1; 8/15/24;
//    // nav changed all of the above to their original values to see how it affects budding. 8/26/24
//    std::cout << "weakened region linear (before budding) = " << scale_linear << std::endl;
//    std::cout << "weakened region bend (before budding) = " << scale_bend << std::endl;
//    std::cout << "weakened region area (before budding) = " << scale_area << std::endl;
//
//    // Scaling factor of the weakend mechanical properties.
//    linearSpringInfoVecs.spring_constant_weak = scale_linear;
//    bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
//    areaTriangleInfoVecs.spring_constant_weak = scale_area;
//
//    bendingTriangleInfoVecs.initial_angle_bud = bendingTriangleInfoVecs.initial_angle[1];
//    std::cout << "equilibrium bending angle of the bud = " << bendingTriangleInfoVecs.initial_angle_bud << std::endl;
//
//
//    /////////////////////////////////////////////////////////////////
//    ////////////////// END OF MEMBRANE RELATED //////////////////////
//    /////////////////////////////////////////////////////////////////
//
//
//
//    //////////////////////////////////////////////////////////////////
//    /////////// IDENTIFYING REGIONS WITH DIFFERENT MECH PROP /////////
//    //////////////////////////////////////////////////////////////////
// /*
//
//                    DV BOUNDARY REGION
//
//    */
//
//    // calculate center of the apical layer of disc
//    int max_upperhem_nodes = 0;
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//
//        if (generalParams.nodes_in_upperhem[i] == 1)
//        {
//            // std::cout<<"node in upperhem = " << generalParams.nodes_in_upperhem[i] <<std::endl;
//            max_upperhem_nodes++;
//        }
//    }
//
//    std::cout << "\n total number of nodes in upperhem = " << max_upperhem_nodes << std::endl;
//
//    double cx_b = 0, cy_b = 0, cz_b = 0;
//    int nLower = 0;
//    double cx_a = 0, cy_a = 0, cz_a = 0;
//    int nUpper = 0;
//
//    for (int i = 0; i < generalParams.maxNodeCount; ++i)
//    {
//        if (generalParams.nodes_in_upperhem[i] == 1)
//        {
//            cx_a += coordInfoVecs.nodeLocX[i];
//            cy_a += coordInfoVecs.nodeLocY[i];
//            cz_a += coordInfoVecs.nodeLocZ[i];
//            ++nUpper;
//        }
//        else if (generalParams.nodes_in_upperhem[i] == -1)
//        {
//            cx_b += coordInfoVecs.nodeLocX[i];
//            cy_b += coordInfoVecs.nodeLocY[i];
//            cz_b += coordInfoVecs.nodeLocZ[i];
//            ++nLower;
//        }
//    }
//    cx_a /= nUpper;
//    cy_a /= nUpper;
//    cz_a /= nUpper; // Apical center
//    cx_b /= nLower;
//    cy_b /= nLower;
//    cz_b /= nLower; // Basal center
//
//    // Boundary node along the x-direction -- this is arbitrary. You can choose any axis to go along.
//
//    double R = 0.0, DR = 0.0;
//
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//        if (generalParams.nodes_in_upperhem[i] != 1)
//            continue;
//        double dx = coordInfoVecs.nodeLocX[i] - cx_a;
//        double dy = sqrt((coordInfoVecs.nodeLocY[i] - cy_a) * (coordInfoVecs.nodeLocY[i] - cy_a));
//        if (dy < 1e-3 && dx > R)
//            R = dx;
//        DR = dy;
//    }
//
//    // for convenience
//
//    generalParams.centerX = cx_a;
//    generalParams.centerY = cy_a;
//    generalParams.centerZ = cz_a;
//
//    std::cout << "boundary node in the x direction (apical layer) =  (" << R << ", " << DR << ") \n"
//              << std::endl;
//
////
//
//    // Identifying regions with different mechanical properties and finding the coundary nodes and edges of the upper hemisphere.
//    generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
//
//    std::cout << "boundaries in upperhem = " << generalParams.boundaries_in_upperhem.size() << std::endl;
//
//    std::vector<int> boundary_edge_list;
//    std::vector<int> boundary_node_list;
//
//    std::cout << "edges2Triangles_1 = " << coordInfoVecs.edges2Triangles_1.size() << std::endl;
//    std::cout << "edges2Triangles_2 = " << coordInfoVecs.edges2Triangles_2.size() << std::endl;
//
//    std::cout << "generalParams.nodes in upperhem size = " << generalParams.nodes_in_upperhem.size() << std::endl;
//
//    for (int i = 0; i < coordInfoVecs.num_edges; i++)
//    {
//        int T1 = static_cast<int>(coordInfoVecs.edges2Triangles_1[i]);
//        int T2 = static_cast<int>(coordInfoVecs.edges2Triangles_2[i]);
//
//        // std::cout<<"it got till here - nav "<< std::endl;// it got till here.
//
//        // Optionally check if the triangle indices are valid.
//        if (T1 < 0 || T2 < 0 || T1 >= (INT_MAX - 1000) || T2 >= (INT_MAX - 1000))
//        {
//            continue;
//        }
//
//        if (T1 == T2 && generalParams.edges_in_upperhem[i] == -1)
//        {
//            generalParams.boundaries_in_lowerhem[i] = 1;
//            boundary_edge_list.push_back(i); // This is to store the total number of boundary edges.
//            int bdry_node1 = static_cast<int>(coordInfoVecs.edges2Nodes_1[i]);
//            int bdry_node2 = static_cast<int>(coordInfoVecs.edges2Nodes_2[i]);
//            boundary_node_list.push_back(bdry_node1);
//            boundary_node_list.push_back(bdry_node2);
//
//            // mark these nodes as boundary or (fixed).
//            generalParams.nodes_in_upperhem[bdry_node1] = -1;
//            generalParams.nodes_in_upperhem[bdry_node2] = -1;
//            coordInfoVecs.isNodeFixed[bdry_node1] = false;
//            coordInfoVecs.isNodeFixed[bdry_node2] = false;
//            
//        }
//        // This is for apical boundary nodes.
//        if (T1 == T2 && generalParams.edges_in_upperhem[i] == 1)
//        { // nav added the second conditional && generalParams.nodes_in_upperhem[i]==1 so that the new apical model would work.
//            generalParams.boundaries_in_upperhem[i] = 1;
//            boundary_edge_list.push_back(i);
//
//            // std::cout<<"it got till here - nav 3 "<< std::endl;
//            int bdry_node1 = static_cast<int>(coordInfoVecs.edges2Nodes_1[i]);
//            int bdry_node2 = static_cast<int>(coordInfoVecs.edges2Nodes_2[i]);
//            boundary_node_list.push_back(bdry_node1);
//            boundary_node_list.push_back(bdry_node2);
//
//            // Optionally mark these nodes as boundary (or fixed).
//            generalParams.nodes_in_upperhem[bdry_node1] = 1; // nav changed this from 0 to 10.
//            generalParams.nodes_in_upperhem[bdry_node2] = 1; // nav changed this from 0 to 10.
//            coordInfoVecs.isNodeFixed[bdry_node1] = false;
//            coordInfoVecs.isNodeFixed[bdry_node2] = false;
//
//            // std::cout<<"it got till here - nav 4 "<< std::endl;
//        }
//        else
//        {
//            generalParams.boundaries_in_upperhem[i] = -1;
//        }
//    }
//
//    std::cout << "size of boundary_edge_list = " << boundary_edge_list.size() << std::endl;
//    std::cout << "size of boundary_node_list (double-counted) = " << boundary_node_list.size() << std::endl;
//
//    // Count the true number of edges in the upper hemisphere.
//    int true_num_edges_in_upperhem = 0;
//    int edges_in_upperhem_COUNT = 0;
//    for (int i = 0; i < coordInfoVecs.num_edges; i++)
//    {
//        if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0)
//        {
//            true_num_edges_in_upperhem += 1;
//            edges_in_upperhem_COUNT += 1;
//        }
//    }
//
//    for (int i = 0; i < coordInfoVecs.num_triangles; i++)
//    {
//        if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_1[i] < 0)
//        {
//            generalParams.triangles_in_upperhem[i] = -1;
//            continue;
//        }
//        else if (coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_2[i] < 0)
//        {
//            generalParams.triangles_in_upperhem[i] = -1;
//            continue;
//        }
//        else if (coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_3[i] < 0)
//        {
//            generalParams.triangles_in_upperhem[i] = -1;
//            continue;
//        }
//
//        int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
//        int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
//        int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];
//
//        if ((aaa + bbb + ccc) == 3)
//        {
//            generalParams.triangles_in_upperhem[i] = 1;
//        }
//        else
//        {
//            generalParams.triangles_in_upperhem[i] = -1;
//        }
//    }
//
//
//
//    std::cout << "INITIAL EDGES IN UPPERHEM = " << edges_in_upperhem_COUNT << std::endl;
//
//    int COUNTING_EDGE = 0;
//    for (int y = 0; y < coordInfoVecs.num_edges; y++)
//    {
//        if (generalParams.edges_in_upperhem_list[y] >= 0)
//        {
//            COUNTING_EDGE += 1;
//        }
//        generalParams.edges_in_upperhem_list_length = COUNTING_EDGE;
//    }
//
//
//    /*
//    
//        SEPARATION OF DV BOUNDARY REGION. THIS CAN ALSO BE PREPROGRAMMED INTO THE STARTING MESH AND SEPARATED AT THE XML PARSER LEVEL.
//    
//    1. Find center of the disc (both apical and basal separately)
//    2. Find boundary of the disc (preferably the distance on the curve)
//    3. You'll go 5% of the distance from the center to the boundary in whatever direction you choose. So start there. Pick an axis to go along. 
//    4. Along that axis until the edge of the disc have 5% on either side of the line drawn from the center be part of the DV boundary. 
//    4. Isolate the nodes in that part and now these are your DV nodes. All edges connecting this region will be a part of the DV boundary region. 
//        Here is a problem I am anticipating with this. There will be asymmetry in connection so is that something that would affect the simulation? I dont think so, considering how 
//        actual tissues will have this king of asymmetry as well because of the way the cells are structured. 
//    5. Now to separate the Origins. Draw a line from the center along the axis orthogonal to the one that the DV boundary lies parallel to. Where it crosses that line is where the 
//       origin of the Dorsal/Ventral regions lies. 
//    6. All the nodes connecting the Dorsal origin and Ventral origin, are included in the origin point of the DV boundary. 
//    
//    */
//    
//    // 1. Compute the new center of the system.
//    newcenterX = 0.0;
//    newcenterY = 0.0;
//    newcenterZ = 0.0;
//
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//        if (generalParams.nodes_in_upperhem[i] == 1)
//        {
//            newcenterX += coordInfoVecs.nodeLocX[i];
//            newcenterY += coordInfoVecs.nodeLocY[i];
//            newcenterZ += coordInfoVecs.nodeLocZ[i];
//        }
//    }
//    newcenterX = newcenterX / generalParams.maxNodeCount;
//    newcenterY = newcenterY / generalParams.maxNodeCount;
//    newcenterZ = newcenterZ / generalParams.maxNodeCount;
//    
//    generalParams.centerX = newcenterX; 
//    generalParams.centerY =newcenterY;
//    generalParams.centerZ = newcenterZ;
//
// 
//    
//
//    /////////////////////////////////////////////////////////////////////
//    ////////////// END OF IDENTIFYING REG. WITH DIFF. MECH PROP /////////
//    /////////////////////////////////////////////////////////////////////
//
//
//   
//    //////////////////////////////////////////////////////////////////////////////////////
//    /////////////////////////////////////////////////////////////////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ////////////////////////////// START OF ACTUAL SIMULATION /////////////////////////////////////////////////////////////////
//    /////////////////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////////
//
//    
//    std::cout<<"\n Starting the actual simulation. \n"<<std::endl;
//    // here recenter the mesh first.
//  
//
//    // nav commented out because we are no longer recentering the mesh. 
//    //    // Compute the new center of the system.
//    newcenterX = 0.0;
//    newcenterY = 0.0;
//    newcenterZ = 0.0;
//
//    for (int i = 0; i < generalParams.maxNodeCount; i++)
//    {
//        if (generalParams.nodes_in_upperhem[i] == 1)
//        {
//            newcenterX += coordInfoVecs.nodeLocX[i];
//            newcenterY += coordInfoVecs.nodeLocY[i];
//            newcenterZ += coordInfoVecs.nodeLocZ[i];
//        }
//    }
//    newcenterX = newcenterX / generalParams.maxNodeCount;
//    newcenterY = newcenterY / generalParams.maxNodeCount;
//    newcenterZ = newcenterZ / generalParams.maxNodeCount;
//    
//    generalParams.centerX = newcenterX; //sumX / nNodes;
//        generalParams.centerY =newcenterY;// sumY / nNodes;
//        generalParams.centerZ = newcenterZ;
//
//    storage->print_VTK_File();
//
//    //// commenting out the following for some tests.  Nav 05/04/2025
//    std::cout << "number of simulation steps = " << generalParams.Tf << std::endl;
//
//
//    
//        double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
//        int nNodes = coordInfoVecs.nodeLocX.size();
//        //int layerflag = 0; // this is the layer that will be excluded. 
//        for (int i = 0; i < nNodes; i++)
//        {  
//            //if (layerflag == 1 && i >= nNodes/2) break;
//            //if (layerflag == -1 && i < nNodes/2) continue;
//                sumX += coordInfoVecs.nodeLocX[i];
//                sumY += coordInfoVecs.nodeLocY[i];
//                sumZ += coordInfoVecs.nodeLocZ[i];
//        }
//        generalParams.centerX = sumX / nNodes;
//        generalParams.centerY = sumY / nNodes;
//        generalParams.centerZ = sumZ / nNodes;
//        
//        //std::cout<<"center at "<<layerflag<<" = ("<<generalParams.centerX<<", "<<generalParams.centerY<<", "<<generalParams.centerZ<<") "<< std::endl;
//
//    
//    std::cout<<"Tf = "<< generalParams.Tf<<std::endl;
//    //    
//    //    for (int i = 0; i<coordInfoVecs.num_edges; i++){
//    //        if (i >= 10) break;
//    //        std::cout<< "initial rest length before strain tensor, at edge = "<< i << " = " <<linearSpringInfoVecs.edge_rest_length[i]<<std::endl;
//    //        
//    //    }
//
//    //double tol = 1e-8;
//    double iter =100;// 1/(generalParams.dt*generalParams.tol);
//    
//    // if you leave ^ this to be 100 it'll take wayy longer for your simulations. Make it at least 10000.
//    
//    ///////////////////////////////////////////////////////////////////
//    // Okay so let's start the simulation loop for the strain tensor //
//    ///////////////////////////////////////////////////////////////////
//    
//    //generalParams.tol = 1e-4;
//    
// //   int layerflag = 0;
//    
//    LambdaField lambda;
//    
//    //    for (int stage = 0; stage<=stages; stage++){
//    
//        double frac = 1.0;///generalParams.Tf; // this fraction should remain what it is. Dont change it. 
//        
//      sumX = 0.0, sumY = 0.0, sumZ = 0.0;
//        
//        
//         nNodes = coordInfoVecs.nodeLocX.size();
//        int layerflag = 0; // this is the layer that will be excluded. 
//        for (int i = 0; i < nNodes; i++)
//        {  
//            //if (layerflag == 1 && i >= nNodes/2) break;
//            //if (layerflag == -1 && i < nNodes/2) continue;
//                sumX += coordInfoVecs.nodeLocX[i];
//                sumY += coordInfoVecs.nodeLocY[i];
//                sumZ += coordInfoVecs.nodeLocZ[i];
//        }
//        generalParams.centerX = sumX / nNodes;
//        generalParams.centerY = sumY / nNodes;
//        generalParams.centerZ = sumZ / nNodes;
//    
//        StrainTensorGPU::buildVertexLambda(generalParams, coordInfoVecs, lambda, frac);
//        
//          //std::cout << "\n[Stage " << stage << "] spring lengths BEFORE strain:\n";
//    //    for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
//    //        if (generalParams.edges_in_upperhem[e] == 0) continue;  // skip vertical
//    //        double L0 = linearSpringInfoVecs.edge_initial_length[e];
//    //        std::cout
//    //            << "  edge[" << e << "]: initial_length = "
//    //            << L0 << "\n";
//    //    }
//        StrainTensorGPU::updateEdgeRestLengths(coordInfoVecs,generalParams, lambda, linearSpringInfoVecs, layerflag);
//        
//        //Solve_Forces();
//        //AdvancePositions(coordInfoVecs, generalParams, domainParams);
//        
//        
//        // === PRINT AFTER ===
//    // std::cout << "[Stage " << stage << "] spring lengths AFTER strain:\n";
//    //    for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
//    //        if (generalParams.edges_in_upperhem[e] == 0) continue;  // skip vertical
//    //        double Lstar = linearSpringInfoVecs.edge_final_length[e];
//    //        std::cout
//    //            << "  edge[" << e << "]: strained_length = "
//    //            << Lstar << "\n";
//    //    }
//    generalParams.tol = 1e-4;
//    
//    // Things to include:
//    /*
//    
//    1. Volume check
//    2. Area Check
//    3. Energy Check
//    4. Connectivity Check
//    
//    */
//              
//              ///////////////////////////////////////////////////////////////////////  
//        // resize edge_rest_length to match the size of the steps between edge_initial and edge_final 
//        // loop for changing the rest lengths between stages. Here the spring length is broken down into steps and passed into the force kernels. 
//        for (int step = 0; step<= iter; step++){
//        
//            for(int i = 0; i<coordInfoVecs.num_edges; i++){ // this can be done from inside the strain tensor file. 
//                
//                double dl = (linearSpringInfoVecs.edge_final_length[i] - linearSpringInfoVecs.edge_initial_length[i])/iter;
//                linearSpringInfoVecs.edge_rest_length[i] += dl;
//            }
//            
//           // double avg_mvt = 999999.0;
//            // generalParams.dx = 0.0;
//            
//           // std::cout<< "Starting relaxation loop."<<std::endl;
//            //int k = 0;
//            
//            int k = relaxUntilConverged(*this);
//            
//           double new_total_energy = linearSpringInfoVecs.linear_spring_energy;
//                                   //areaTriangleInfoVecs.area_triangle_energy +
//                                   //bendingTriangleInfoVecs.bending_triangle_energy; 
//                                   
//            std::cout<<"Energy = "<<new_total_energy<<", iterations = "<<k<<", = Average Movement = "<<generalParams.dx<<std::endl;
//                        
//           // std::cout<<"Exiting relaxation loop. Average movement = "<<generalParams.dx<<", and k = "<< k<<" at step = "<<step<<std::endl;//" and stage = "<<stage<<std::endl;     
//            if (step > 0 && step % 2 == 0) {
//                
//                storage->print_VTK_File(); // this is just for now. To get better simulation results we're gonna view them more frequently later. 
//            }
//
//          }
//    //////////////////////////////////////
//    //            
//           // do{
//             //   generalParams.dx = 0.0; 
//                
//    //                double x1[coordInfoVecs.nodeLocX.size()];
//    //                double y1[coordInfoVecs.nodeLocX.size()];
//    //                double z1[coordInfoVecs.nodeLocX.size()];
//    //                for (int n = 0; n<coordInfoVecs.nodeLocX.size(); n++){ // this can be handled from inside the nodeAdvance.cu
//    //                    x1[n] = coordInfoVecs.nodeLocX[n];
//    //                    y1[n] = coordInfoVecs.nodeLocY[n];
//    //                    z1[n] = coordInfoVecs.nodeLocZ[n];
//    //                } 
//                
//                //avg_mvt = 0.0;
//               // Solve_Forces();
//                // AdvancePositions(coordInfoVecs, generalParams, domainParams);
//                
//                
//    //                double x2[coordInfoVecs.nodeLocX.size()];
//    //                double y2[coordInfoVecs.nodeLocX.size()];
//    //                double z2[coordInfoVecs.nodeLocX.size()];
//    //                
//    //                for (int n = 0; n<coordInfoVecs.nodeLocX.size(); n++){ // same as above. 
//    //                    x2[n] = coordInfoVecs.nodeLocX[n];
//    //                    y2[n] = coordInfoVecs.nodeLocY[n];
//    //                    z2[n] = coordInfoVecs.nodeLocZ[n];
//    //                    
//    //                    generalParams.dx+= sqrt((x1[n]-x2[n])*(x1[n]-x2[n])+(y1[n]-y2[n])*(y1[n]-y2[n])+(z1[n]-z2[n])*(z1[n]-z2[n]));
//    //                }
//    //                
//    //                
//                //avg_mvt=generalParams.dx;   
//              //  k++;
//                
//               // if (k == 4500) break;
//            //}while (generalParams.dx>generalParams.tol);
//            
//            
//            //if (step == 100 || step == 200 || step == 300 || step == 400 || step == 500 || step == 600 || step == 700 || step == 800 || step == 900 || step == 1000 || step == 2000 || step == 5000 || step == 7000 || step == 9000 ){ // badly done ik but I'll fix it later. 
//            
//      //  }
//        
//        storage->print_VTK_File();
//        
//  //  }
//    
//    // strain tensor acts on edge_inital to give edge_final 
//};
//
//
//
//
//
//// Function to assign the shared pointer to storage.
//void System::assignStorage(std::shared_ptr<Storage> _storage)
//{
//    storage = _storage;
//};
//
//// Function to set the weak pointer to the SystemBuilder.
//void System::set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr)
//{
//    weak_bld_ptr = _weak_bld_ptr;
//};
//
//// Function to initialize memory for thrust vectors and set coordInfoVecs values from input.
//void System::initializeSystem(HostSetInfoVecs & hostSetInfoVecs)
//{
//    std::cout << "Initializing" << std::endl;
//
//    // Set the max node count, edge count and triangle count.
//    generalParams.maxNodeCount = hostSetInfoVecs.nodeLocX.size();
//    coordInfoVecs.num_edges = hostSetInfoVecs.edges2Nodes_1.size();
//    coordInfoVecs.num_triangles = hostSetInfoVecs.triangles2Nodes_1.size();
//
//    std::cout << "num nodes: " << generalParams.maxNodeCount << std::endl;
//    std::cout << "num edges: " << coordInfoVecs.num_edges << std::endl;
//    std::cout << "num elems: " << coordInfoVecs.num_triangles << std::endl;
//    // Allocate memory for various vectors using preallocated memory size.
//    int mem_prealloc = 4;
//
//    // Resize and initialize the following various vectors.
//    coordInfoVecs.scaling_per_edge.resize(mem_prealloc * coordInfoVecs.num_edges, 0.0);
//    hostSetInfoVecs.scaling_per_edge.resize(coordInfoVecs.scaling_per_edge.size(), 0.0);
//
//    coordInfoVecs.soln_per_triangle.resize(mem_prealloc * coordInfoVecs.num_triangles, INT_MAX);
//    coordInfoVecs.b_per_triangle.resize(mem_prealloc * coordInfoVecs.num_triangles, INT_MAX);
//
//    coordInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), false);
//    coordInfoVecs.prevNodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.prevNodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.prevNodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//
//    coordInfoVecs.prevNodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.prevNodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.prevNodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//
//    coordInfoVecs.nodeLocX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.nodeLocY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//    coordInfoVecs.nodeLocZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//
//    // coordInfoVecs.nodeVelX.resize(mem_prealloc*hostSetInfoVecs.nodeVelX.size(), 0.0);
//    // coordInfoVecs.nodeVelY.resize(mem_prealloc*hostSetInfoVecs.nodeVelY.size(), 0.0);
//    // coordInfoVecs.nodeVelZ.resize(mem_prealloc*hostSetInfoVecs.nodeVelZ.size(), 0.0);
//
//    coordInfoVecs.nodeForceX.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
//    coordInfoVecs.nodeForceY.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
//    coordInfoVecs.nodeForceZ.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size(), 0.0);
//
//    coordInfoVecs.triangles2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
//    coordInfoVecs.triangles2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
//    coordInfoVecs.triangles2Nodes_3.resize(mem_prealloc * coordInfoVecs.num_triangles);
//
//    coordInfoVecs.triangles2Edges_1.resize(mem_prealloc * coordInfoVecs.num_triangles);
//    coordInfoVecs.triangles2Edges_2.resize(mem_prealloc * coordInfoVecs.num_triangles);
//    coordInfoVecs.triangles2Edges_3.resize(mem_prealloc * coordInfoVecs.num_triangles);
//
//    coordInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    coordInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    coordInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//
//    hostSetInfoVecs.triangles2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    hostSetInfoVecs.triangles2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//    hostSetInfoVecs.triangles2Triangles_3.resize(mem_prealloc * coordInfoVecs.num_triangles, -INT_MAX);
//
//    coordInfoVecs.edges2Nodes_1.resize(mem_prealloc * coordInfoVecs.num_edges);
//    coordInfoVecs.edges2Nodes_2.resize(mem_prealloc * coordInfoVecs.num_edges);
//
//    coordInfoVecs.edges2Triangles_1.resize(mem_prealloc * coordInfoVecs.num_edges);
//    coordInfoVecs.edges2Triangles_2.resize(mem_prealloc * coordInfoVecs.num_edges);
//
//    coordInfoVecs.nndata1.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata2.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata3.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata4.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata5.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata6.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata7.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata8.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.nndata9.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.SurfaceNormalX.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.SurfaceNormalY.resize(mem_prealloc * generalParams.maxNodeCount);
//    coordInfoVecs.SurfaceNormalZ.resize(mem_prealloc * generalParams.maxNodeCount);
//
//    generalParams.nodes_in_upperhem.resize(mem_prealloc * generalParams.maxNodeCount);
//    generalParams.triangles_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_triangles);
//    generalParams.edges_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges);
//    generalParams.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
//    generalParams.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
//    generalParams.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
//
//    hostSetInfoVecs.nodes_in_upperhem.resize(generalParams.nodes_in_upperhem.size());
//    generalParams.nodes_in_upperhem = hostSetInfoVecs.nodes_in_upperhem;
//    hostSetInfoVecs.triangles_in_upperhem.resize(generalParams.triangles_in_upperhem.size());
//    hostSetInfoVecs.edges_in_upperhem.resize(generalParams.edges_in_upperhem.size());
//    generalParams.edges_in_upperhem = hostSetInfoVecs.edges_in_upperhem;
//    hostSetInfoVecs.edges_in_upperhem_list.resize(mem_prealloc * coordInfoVecs.num_edges);
//    hostSetInfoVecs.boundaries_in_upperhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
//    hostSetInfoVecs.boundaries_in_lowerhem.resize(mem_prealloc * coordInfoVecs.num_edges, -1);
//
//    hostSetInfoVecs.nodes2Triangles_1.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_2.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_3.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_4.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_5.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_6.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_7.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_8.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    hostSetInfoVecs.nodes2Triangles_9.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//
//    // Part 17
//
//    // Resize vectors to allocate memory for nodes-to-triangles mapping.
//    coordInfoVecs.nodes2Triangles_1.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_2.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_3.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_4.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_5.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_6.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_7.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_8.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//    coordInfoVecs.nodes2Triangles_9.resize(mem_prealloc * generalParams.maxNodeCount, -INT_MAX);
//
//    // Copy nodes-to-triangles mapping information from hostSetInfoVecs to coodInfoVecs and others.
//    thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(), coordInfoVecs.nodes2Triangles_1.end(), hostSetInfoVecs.nodes2Triangles_1.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(), coordInfoVecs.nodes2Triangles_2.end(), hostSetInfoVecs.nodes2Triangles_2.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(), coordInfoVecs.nodes2Triangles_3.end(), hostSetInfoVecs.nodes2Triangles_3.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(), coordInfoVecs.nodes2Triangles_4.end(), hostSetInfoVecs.nodes2Triangles_4.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(), coordInfoVecs.nodes2Triangles_5.end(), hostSetInfoVecs.nodes2Triangles_5.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(), coordInfoVecs.nodes2Triangles_6.end(), hostSetInfoVecs.nodes2Triangles_6.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(), coordInfoVecs.nodes2Triangles_7.end(), hostSetInfoVecs.nodes2Triangles_7.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(), coordInfoVecs.nodes2Triangles_8.end(), hostSetInfoVecs.nodes2Triangles_8.begin());
//    thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(), coordInfoVecs.nodes2Triangles_9.end(), hostSetInfoVecs.nodes2Triangles_9.begin());
//
//    // copy info to GPU
//    std::cout << "Copying" << std::endl;
//    thrust::copy(hostSetInfoVecs.isNodeFixed.begin(), hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());
//
//    // Print information about fixed nodes in hostSetInfoVecs and coordInfoVecs.
//    std::cout << "fixed_node_in_host: " << std::endl;
//    for (int k = 0; k < hostSetInfoVecs.isNodeFixed.size(); k++)
//    {
//    }
//    std::cout << "end_of_fixed_node_host_printout" << std::endl;
//    std::cout << "fixed_node_in_device: " << std::endl;
//    for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++)
//    {
//    }
//    std::cout << "end_of_fixed_node_device_printout" << std::endl;
//    std::cout << "size of host fixed " << hostSetInfoVecs.isNodeFixed.size() << std::endl;
//    std::cout << "size of device fixed " << coordInfoVecs.isNodeFixed.size() << std::endl;
//
//    // initialize various vectors with zeros or values from hostSetInfoVecs.
//    //  Fill operations for other nodeForce and prevNodeForce vectors.
//    thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
//    thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
//    thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
//
//    thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
//    thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
//    thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);
//
//    // Copy node locations and other related information from hostSetInfoVecs to coordInfoVecs and other copy operations for triangles, edges and other related vectors.
//    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.prevNodeLocX.begin());
//    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.prevNodeLocY.begin());
//    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin());
//
//    thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocX.begin());
//    thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocY.begin());
//    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.nodeLocZ.begin());
//
//    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(), hostSetInfoVecs.triangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin());
//    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(), hostSetInfoVecs.triangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin());
//    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(), hostSetInfoVecs.triangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin());
//
//    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(), hostSetInfoVecs.triangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin());
//    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(), hostSetInfoVecs.triangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin());
//    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(), hostSetInfoVecs.triangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin());
//
//    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(), hostSetInfoVecs.edges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin());
//    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(), hostSetInfoVecs.edges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin());
//
//    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(), hostSetInfoVecs.edges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin());
//    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(), hostSetInfoVecs.edges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin());
//
//    thrust::copy(hostSetInfoVecs.nndata1.begin(), hostSetInfoVecs.nndata1.end(), coordInfoVecs.nndata1.begin());
//    thrust::copy(hostSetInfoVecs.nndata2.begin(), hostSetInfoVecs.nndata2.end(), coordInfoVecs.nndata2.begin());
//    thrust::copy(hostSetInfoVecs.nndata3.begin(), hostSetInfoVecs.nndata3.end(), coordInfoVecs.nndata3.begin());
//    thrust::copy(hostSetInfoVecs.nndata4.begin(), hostSetInfoVecs.nndata4.end(), coordInfoVecs.nndata4.begin());
//    thrust::copy(hostSetInfoVecs.nndata5.begin(), hostSetInfoVecs.nndata5.end(), coordInfoVecs.nndata5.begin());
//    thrust::copy(hostSetInfoVecs.nndata6.begin(), hostSetInfoVecs.nndata6.end(), coordInfoVecs.nndata6.begin());
//    thrust::copy(hostSetInfoVecs.nndata7.begin(), hostSetInfoVecs.nndata7.end(), coordInfoVecs.nndata7.begin());
//    thrust::copy(hostSetInfoVecs.nndata8.begin(), hostSetInfoVecs.nndata8.end(), coordInfoVecs.nndata8.begin());
//    thrust::copy(hostSetInfoVecs.nndata9.begin(), hostSetInfoVecs.nndata9.end(), coordInfoVecs.nndata9.begin());
//
//    // Resize and initialize the 'u' vector.
//    coordInfoVecs.u.resize(mem_prealloc * coordInfoVecs.num_triangles);
//
//    // Part 18
//
//    // Allocate memory for additiional data structures.
//
//    // Area triangle info vec.
//    // Number of area springs is the number of triangles
//    std::cout << "Mem" << std::endl;
//    // Allocate memory for temporary node information in unreduced form for area springs
//    areaTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//
//    // Allocate memory for temporary node information in reduced form for area springs.
//    areaTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//    areaTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
//
//    // beinding triangle info vec
//    // num bending springs is the number of times each edge is between two triangles.
//    bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;
//
//    // Allocate memory for temporary node information in unreduced form for bending springs.
//    bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//
//    // Allocate memory for temporary node information in reduced form for bending springs.
//    bendingTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    bendingTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
//    
//    bendingTriangleInfoVecs.initial_angle.resize(coordInfoVecs.num_triangles);
//    thrust::fill(bendingTriangleInfoVecs.initial_angle.begin(), bendingTriangleInfoVecs.initial_angle.end(), 0.0087);
//        
//    // linear springs info vectors.
//    //  Allocate memory for temporary node information in unreduced form for linear springs.
//    linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//
//    // Allocate memory for temporary node information in reduced form for bending springs.
//    linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//    linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc * linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
//
//    // Clear edge_initial_length vector for linear springs.
//    // linearSpringInfoVecs.edge_initial_length.clear();
//    // linearSpringInfoVecs.edge_rest_length.clear();
//
//    // linearSpringInfoVecs.edge_rest_length.resize(hostSetInfoVecs.edge_rest_length.size());
//    linearSpringInfoVecs.edge_final_length.resize(hostSetInfoVecs.edge_initial_length.size());
//    linearSpringInfoVecs.edge_initial_length = hostSetInfoVecs.edge_initial_length;
//    linearSpringInfoVecs.edge_final_length = hostSetInfoVecs.edge_initial_length;
//    
//    
//    std::cout << "host edge_initial_length size = " << hostSetInfoVecs.edge_initial_length.size() << std::endl;
//    std::cout << "device edge_initial_length size = " << linearSpringInfoVecs.edge_initial_length.size() << std::endl;
//
//    //  for (int e = 0; e < coordInfoVecs.num_edges; ++e) {
//    //    int i = coordInfoVecs.edges2Nodes_1[e];
//    //    int j = coordInfoVecs.edges2Nodes_2[e];
//    //    double dx = hostSetInfoVecs.nodeLocX[j] - hostSetInfoVecs.nodeLocX[i];
//    //    double dy = hostSetInfoVecs.nodeLocY[j] - hostSetInfoVecs.nodeLocY[i];
//    //    double dz = hostSetInfoVecs.nodeLocZ[j] - hostSetInfoVecs.nodeLocZ[i];
//    //    double dist = sqrt(dx*dx + dy*dy + dz*dz);
//    //    //hostSetInfoVecs.edge_initial_length.push_back(dist);    // already done for initial
//    //    hostSetInfoVecs.edge_rest_length.push_back(dist);
//    //    //std::cout<< "edge_rest_length = " << hostSetInfoVecs.edge_initial_length[i]<<std::endl; in the current data structure it gave me 1021 edges. That's good. Now that they have been initialized I should start changing the rest lengths.
//    //    }
//
//    // thrust::copy(linearSpringInfoVecs.edge_rest_length.begin(),
//    //              linearSpringInfoVecs.edge_rest_length.end(),
//    //              hostSetInfoVecs.edge_rest_length.begin());
//
//    linearSpringInfoVecs.edge_rest_length = hostSetInfoVecs.edge_initial_length;
//
//    //linearSpringInfoVecs.edge_rest_length.resize(1/(generalParams.dt*generalParams.tol)) //= hostSetInfoVecs.edge_rest_length;
//
//    std::cout << "host edge_rest_length size = " << hostSetInfoVecs.edge_rest_length.size() << std::endl;
//    std::cout << "device edge_rest_length size = " << linearSpringInfoVecs.edge_rest_length.size() << std::endl;
//
//    //linearSpringInfoVecs.edge_final_length.resize(coordInfoVecs.num_edges);
//
//    //linearSpringInfoVecs.edge_final_length = linearSpringInfoVecs.edge_initial_length;
//
//    // This loop is to test out whether the edge_rest_length is being initialized properly.
//    // for (int i = 0; i < hostSetInfoVecs.edge_rest_length.size(); i++) {
//    //    std::cout<< "edge_rest_length # "<< i <<" = "<< hostSetInfoVecs.edge_rest_length[i]<<std::endl;
//    //    std::cout<< "edge_initial_length # "<< i <<" = "<< hostSetInfoVecs.edge_initial_length[i]<<std::endl;
//    //}
//
//    // Resize the hostSetInfoVecs for data transfer between host and device.
//    hostSetInfoVecs.isNodeFixed.resize(mem_prealloc * hostSetInfoVecs.nodeLocX.size());
//
//    hostSetInfoVecs.nodeLocX.resize(coordInfoVecs.nodeLocX.size());
//    hostSetInfoVecs.nodeLocY.resize(coordInfoVecs.nodeLocX.size());
//    hostSetInfoVecs.nodeLocZ.resize(coordInfoVecs.nodeLocX.size());
//    std::cout << "Host_nodeLocX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;
//
//    // hostSetInfoVecs.nodeVelX.resize(coordInfoVecs.nodeVelX.size());
//    // hostSetInfoVecs.nodeVelY.resize(coordInfoVecs.nodeVelY.size());
//    // hostSetInfoVecs.nodeVelZ.resize(coordInfoVecs.nodeVelZ.size());
//
//    hostSetInfoVecs.nodeForceX.resize(coordInfoVecs.nodeLocX.size());
//    hostSetInfoVecs.nodeForceY.resize(coordInfoVecs.nodeLocX.size());
//    hostSetInfoVecs.nodeForceZ.resize(coordInfoVecs.nodeLocX.size());
//    std::cout << "Host_nodeForceX size = " << hostSetInfoVecs.nodeLocX.size() << std::endl;
//
//    hostSetInfoVecs.triangles2Nodes_1.resize(coordInfoVecs.triangles2Nodes_1.size());
//    hostSetInfoVecs.triangles2Nodes_2.resize(coordInfoVecs.triangles2Nodes_2.size());
//    hostSetInfoVecs.triangles2Nodes_3.resize(coordInfoVecs.triangles2Nodes_3.size());
//    std::cout << "Host_triangles2Nodes size = " << hostSetInfoVecs.triangles2Nodes_1.size() << std::endl;
//
//    hostSetInfoVecs.triangles2Edges_1.resize(coordInfoVecs.triangles2Edges_1.size());
//    hostSetInfoVecs.triangles2Edges_2.resize(coordInfoVecs.triangles2Edges_2.size());
//    hostSetInfoVecs.triangles2Edges_3.resize(coordInfoVecs.triangles2Edges_3.size());
//    std::cout << "Host_triangles2Edges size = " << hostSetInfoVecs.triangles2Edges_1.size() << std::endl;
//
//    hostSetInfoVecs.edges2Nodes_1.resize(coordInfoVecs.edges2Nodes_1.size());
//    hostSetInfoVecs.edges2Nodes_2.resize(coordInfoVecs.edges2Nodes_2.size());
//    std::cout << "Host_edges2Nodes size = " << hostSetInfoVecs.edges2Nodes_1.size() << std::endl;
//
//    hostSetInfoVecs.edges2Triangles_1.resize(coordInfoVecs.edges2Triangles_1.size());
//    hostSetInfoVecs.edges2Triangles_2.resize(coordInfoVecs.edges2Triangles_2.size());
//    std::cout << "Host_edges2Triangles size = " << hostSetInfoVecs.edges2Triangles_1.size() << std::endl;
//
//    hostSetInfoVecs.nndata1.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata2.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata3.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata4.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata5.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata6.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata7.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata8.resize(mem_prealloc * generalParams.maxNodeCount);
//    hostSetInfoVecs.nndata9.resize(mem_prealloc * generalParams.maxNodeCount);
//
//    // Print message indicating the system is ready.
//    std::cout << "System Ready" << std::endl;
//
//    // Allocate memory for buckets.
//    auxVecs.id_bucket.resize(generalParams.maxNodeCount);
//    auxVecs.id_value.resize(generalParams.maxNodeCount);
//    auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
//    auxVecs.id_value_expanded.resize(27 * (generalParams.maxNodeCount));
//};
//
