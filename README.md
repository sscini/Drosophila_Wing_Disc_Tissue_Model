# BuddingCode_UCR
Budding model. The initial condition can be a closed system or an open system.

Case 1. yeast budding

      Starting from an initial spherical (or non-spherical) triangulated model cell, a selected region can undergo changes in mechanical properties leading to deformation. 
      Additional triangles can be introduced into the selected region to simulate local growth.

Case 2. Viral budding

      Starting from an initial flat triangulated sheet representing a piece of membrane, a selected region can undergo changes in mechanical properties leading to deformation.
      Interaction with a single viral particle or a cluster of viral particles can also be accommodated.


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

The general flow of simulation steps can be found at "void System::solveSystem()" in System.cu.

For particular functions such as linear spring, please refer to LinearSpring.cu and LinearSpring.h files. The same applies to bending spring and area springs.

For edge-swap algorithm and general data structure manipulation functions, please refer to Edgeswap_test.cpp and Edgeswap_test.h.

To change the name of saved animation and data output, please refer to Storage.cpp and Storage.h.

To change the simulation job title and to some extend simulation time step size, please refer to SBATCH.sh.

Initial data structure (built via MATLAB functions) is located in Data_Structure.xml.

Overall flow of the simulation steps:

I. Initialization of global parameters and data structures.

II. Run a predetermined number of relaxation steps of the model system to attain quasi-steady state.

III. Start the actual simulation:

      1. Update parameters (if necessary).
      2. Run a predetermined number of relaxation steps or dynamical number of relaxation steps depending on simulation types (molecular dynamics vs energy minimization).
      3. Run edge-swap algorithm (if applicable).
      4. Repeat (1-3) for a number of times, and then test for growth (if applicable).
      5. Repeat (1-4) until simulation terminates.

To run the simulation on UCR HPCC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load extra; module load GCCcore/4.9.3; module load cuda/9.1
4. type: make (Or make -j N, N can be 2,3,4,....,or 12. But this should only be done if you are using an interactive gpu session. See UCR HPCC website for detail)
5. type: sbatch -p gpu --gres=gpu:1 --time=144:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh 

###################################################

If the above steps do not work, try the following (ported from Episcale code): (note: These instructions do not work anymore. 
To run simulation on UCR HPCC cluster: 
After uploading all the files under this repository, follow the command below.
1. cd [The folder where all your files are placed in]
2.  module load singluarity
3.  module load centos
4.  centos.sh
5.  module load extra; module load GCC; module load cuda;
6.  (IF it is the first time compiling) module load cmake;
7.  (IF it is the first time compiling) cmake .
8.  make
9.  After compilation, enter: exit
10. sbatch -p gpu --gres=gpu:1 --time=X:00:00 SBATCH.sh (X here is the number of hours you want to keep the simulation running)

###################################################

Current working simulation steps (modified from above):
To run simulation on UCR HPCC cluster:
1. *either clone or upload file onto HPCC user storage.*
2. *change directory to where the "SBATCH_try_this_one_if_the_original_does_not_work.sh" is located using*:
      a. cd [*folder name with "SBATCH_try_this_one_if_the_original_does_not_work.sh" file*]
3. module load singularity
4. module load centos
5. centos.sh
6. module load extra; module load GCCcore/4.9.3; module load cuda/9.1;
7. module load cmake
8. make
9. *After compilation completes, enter*: exit;
10. sbatch -p gpu --gres=gpu:1 --time=X:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh (*X here is the number of hours you want to keep the simulation running*)

Current working module list:






Current steps taken:

1. Submitting job in case 2 with flat viral model. Initial conditions are the only things that have been changed so far. 

1. Need to change condition on line 815 and 830 to reflect T1==T2 for boundaries. 
Current submission does not produce animation results. 


For the initial angle for bending traingles - cosine beinding potential defined by the normal direction of the triangle constant - change this in the code. He has defined it as being 0.08.. You can change the initial angle. Make it relatively drastic like 0.5 radians. The idea is that now if you let a code run it should relax to some local energy minmum based on the angle. So if youre able to run it with this then the membrane is ready to go. The paramater to change is bendingTriangleInfoVecs.initial_angle.


6/5/24
As per Kevin's email I will be modifying the initial parameters of the simulation. 
1. The Tau coefficient in Data_structure_membr_freebdry.xml has been changed from 452.0 to 0.0. This is to test if any of the mechanical force calculation and/or associated functions is giving an error. 

2. The error "Too few edges are available for the chosen edge swap sample size" in the printout can be fixed by going to System.cu around line 1632 (on github) and replace the if condition "if ((current_edge_to_tip_height) <= generalParams..." to "if ((current_edge_to_tip_height) > 0.0)".

3. Comment out line 1893 to 1901 (on github) in Utilities.cpp. Those lines are the one producing the "error: triangles2Triangles_host_vec() ..." message. You don't need those at the moment.

4. in the data structures file rmin has been changed to 0.15 because that was what was found commented out in the code. If that value does not work I will calculate Rmin myself and input it. 

6/5/24

Since it is still not working I am going to turn the growth algorithm on again. I feel like the problems are coming about because of me truning it off. 


6/10/24

Lin constant in the initial parameters file changes the line constant. 

Bend constant being made 1.0 makes it a completely flat sheet. 

8/17/24
NOTE: I have saved the flat sheet (circular) that has been generated. 8/17/2024

Turning on the turgor pressure makes one tiny corner curl up and around. That needs to be fixed. 

Changed the septin ring stiffness to 50.0 from 0.0 (parameter name = generalParams.line_tension_constant). 

Changed the septin ring length to generalParams.length_scale = 1.0; from 0.0. 

Look into what cell height is in the simulation output. 

The output received after turning on septin ring is very weird and the name is FINALLY_FLAT_turgor_on_septin_ring_50_stiffness_00000.vtk. The septin ring is causing it to be very butt-like and needs to be changed. 
