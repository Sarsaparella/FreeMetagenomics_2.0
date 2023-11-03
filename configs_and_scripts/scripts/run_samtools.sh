#!/bin/bash -l 
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J samtools

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 00:20:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1

# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=64

# Run the executable named my exe
# Removing reads with insufficent mapping
samtools view -h ASE_67_05M_mapped.sam | awk -F'\t' '{ if ($1 ~ /^@/ || ($6 ~ /^[0-9]+M$/ && $5 >= 45 && (($4 + $5) / length($10)) >= 0.8)) print }' | samtools view -bS -F 4 - > ASE_67_05M_mapped_filtered.bam
