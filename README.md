# Drosophila_Wing_Disc_Eversion_Tissue_Model

Case 1. Drosophila Wing Disc Eversion

      Starting from an initial spherical triangulated dome representing a piece of tissue on the wing disc, a selected region can undergo changes in mechanical properties leading to deformation.
      Chemical diffusion in the model can also be accomodated.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

The general flow of simulation steps can be found at "void System::solveSystem()" in System.cu.

For particular functions such as linear spring, please refer to LinearSpring.cu and LinearSpring.h files. The same applies to bending spring and area springs.

To change the name of saved animation and data output, please refer to Storage.cpp and Storage.h.

To change the simulation job title and to some extend simulation time step size, please refer to SBATCH_try_this_one_if_the_other_does_not_work.sh.

Initial data structure (built via MATLAB functions) is located in Data_Structures/Data_Structure.xml.

Overall flow of the simulation steps:

I. Initialization of global parameters and data structures.

II. Run a predetermined number of relaxation steps of the model system to attain quasi-steady state.

III. Start the actual simulation:

      1. Update parameters (if necessary).
      2. Run a predetermined number of relaxation steps or dynamical number of relaxation steps depending on simulation types (molecular dynamics vs energy minimization).
      3. Run strain field for a certain number of sets on set layers.
      4. Repeat (1-3) for a number of times, and then test for deformation (if applicable).
      5. Repeat (1-4) until simulation terminates.

To run the simulation on UCR HPCC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load cuda; module load cmake
4. type: make (Or make -j N, N can be 2,3,4,....,or 12. But this should only be done if you are using an interactive gpu session. See UCR HPCC website for detail)
5. type: sbatch -p gpu --gres=gpu:1 --time=01:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh 

###################################################

To run simulation on CRC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load cuda cmake
4. type: make
5. Once compliation has completed run: qsub run.sh
