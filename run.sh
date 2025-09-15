#!/bin/csh

#to run without script you must type this line into the terminal:
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/"

#$ -M  nsher012@ucr.edu# Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#s -pe smp 6         #specifies threads??? maybe
#$ -N  "dsp_test_2_circular" # Specify job name
#$ -t 1       #specify number of data input files


set data = ( DataStructures/Data_Structure.xml )

module purge
module load gcc/12.2.0
module load cuda/12.8


echo -n "It is currently: ";date
echo -n "I am logged on as nsher012";whoami
echo -n "This computer is called Lernie";hostname
echo -n "I am currently in the directory /rhome/nsher012/Dsp_wing_eversion/07-03-25-codes_to_change/";pwd


./tissue-model  -solve_time=5 -dt=0.01 $data[${SGE_TASK_ID}]

