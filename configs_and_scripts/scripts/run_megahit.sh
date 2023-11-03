#!/bin/bash -l 
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J megahit

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 07:00:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=200

# Run the executable named my exe
# Assembling reads into scaffolds 
megahit -1 ASE_67_05M_forward_unmatched.fastq -2 ASE_67_05M_forward_unmatched.fastq --min-contig-len 1000 -o megahit_output_folder --num-cpu-threads 80
