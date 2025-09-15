#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --output=_Does_this_work_test_wl3-0hapf_# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=nsher012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="dsp_test_218_nodes_"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load singularity
export SINGULARITY_NV=1
module load centos

centos.sh "module load extra; module load GCC/8.3.0; module load cuda/11.2; ./tissue-model -dt=0.01 Data_Structures/Data_Structure.xml"
