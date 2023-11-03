#!/bin/bash -l
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A naiss2023-22-585

# The name of the script
#SBATCH -J magpipe

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 00:07:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1
# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=60


# Run the executable named my exe
snakemake -s magpipe_snake.smk --configfile /cfs/klemming/scratch/e/erikzhi/project/magpipe-master/configs/test_config.yaml --use-conda --profile /cfs/klemming/scratch/e/erikzhi/project/magpipe-master/snakes/cluster_configs --latency-wait 20 
