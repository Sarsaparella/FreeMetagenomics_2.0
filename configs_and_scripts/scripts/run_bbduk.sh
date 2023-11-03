#!/bin/bash -l 
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J bbduk_filter

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 00:05:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=250

# Run the executable named my exe
# Removing low-quality reads
bbduk.sh in1=ASE_67_05M-QUALITY_PASSED_R1.fastq in2=ASE_67_05M-QUALITY_PASSED_R2.fastq out1=ASE_67_05M_R1_filtered.fastq out2=ASE_67_05M_R2_filtered.fastq qtrim=rl trimq=14 maq=20 maxns=0 minlength=45
