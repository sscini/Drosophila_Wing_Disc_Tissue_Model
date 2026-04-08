#!/bin/csh

#to run without script you must type this line into the terminal:
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/afs/crc.nd.edu/x86_64_linux/m/matlab/R2017b/extern/bin/glnxa64/"

#$ -M  scini@nd.edu# Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#s -pe smp 6         #specifies threads??? maybe
#$ -N  "build_test_040826" # Specify job name
#$ -t 1       #specify number of data input files


set data = ( Data_Structures/Data_Structure.xml )

module purge
module load cuda

echo -n "It is currently: ";date
echo -n "I am logged on as scini";whoami
echo -n "This computer is called Lernie";hostname
# echo -n "I am currently in the directory /rhome/scini/Dsp_wing_eversion/07-03-25-codes_to_change/";pwd


./tissue-model  -solve_time=5 -dt=0.01 $data[${SGE_TASK_ID}]

