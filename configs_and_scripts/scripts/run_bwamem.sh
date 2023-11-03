#!/bin/bash -l 
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J bwamem

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 01:30:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=64

# Run the executable named my exe
# Mapping reads onto the reference
bwa mem -a ASE_67_05M_final_contigs.min1000.fa ASE_67_05M_R1_phix_unmatched.fastq ASE_67_05M_R2_phix_unmatched.fastq -t 10 > ASE_67_05M_mapped.sam
