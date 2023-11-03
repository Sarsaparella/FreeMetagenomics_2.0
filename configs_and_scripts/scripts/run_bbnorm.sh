#!/bin/bash -l 
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J bbnorm_norm_2

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 00:25:00

# Number of Nodes 
#SBATCH --nodes=2

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=200

# Run the executable named my exe
# Normalise reads
bbnorm.sh in=ASE_67_05M_R2_phix_unmatched.fastq out=ASE_67_05M_reverse_unmatched.fastq target=40 mindepth=0
