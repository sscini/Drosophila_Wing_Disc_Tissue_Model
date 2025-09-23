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

