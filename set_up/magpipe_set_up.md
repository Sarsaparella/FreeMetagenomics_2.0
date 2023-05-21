# How to start working with the `magpipe` pipeline

Here is the detailed guideline on how to clone the `magpipe` pipeline with all necessary arrangements.

## 1. Installing magpipe

This is the step that was described in the [README](https://github.com/SushiLab/magpipe/blob/master/README.md) 
file of the original repository.

```bash
conda create -n mapgipe
conda activate magpipe
conda install snakemake
```

```bash
git clone git@github.com:SushiLab/magpipe.git
cd magpipe/magpipe
python -m pip install -r requirements.txt -e .
```

## 2. Investigate the structure of the repository

**Structure of the `magpipe` repository as described in its README:**

- **[configs](https://github.com/SushiLab/magpipe/tree/master/configs)**: config files used by the launchers to run the pipeline;
- **[envs](https://github.com/SushiLab/magpipe/tree/master/envs)**: conda environments for the rules;
- **[figures](https://github.com/SushiLab/magpipe/tree/master/figures)**: R code behind the figures;
- **[launchers](https://github.com/SushiLab/magpipe/tree/master/launchers)**: bash scripts to run the pipeline (on sge etc);
- **[magpipe](https://github.com/SushiLab/magpipe/tree/master/magpipe)**: python module for the pipeline;
- **[resources](https://github.com/SushiLab/magpipe/tree/master/resources)**: some fixed input information on methods and datasets;
- **[rules](https://github.com/SushiLab/magpipe/tree/master/rules)**: the rules of the snakemake pipeline;
- **[sandbox](https://github.com/SushiLab/magpipe/tree/master/sandbox)**: scripts used for postprocessing and analyses;
- **[scripts](https://github.com/SushiLab/magpipe/tree/master/scripts)**: ad-hoc python scripts used by the pipeline;
- **[snakes](https://github.com/SushiLab/magpipe/tree/master/snakes)**: snakemake files.


It is better to explore these folders bu yourself, but the most important ones for us right now are **snakes**, **resources** and **configs**.

## 3. Preparations for running pipeline ( ￣ー￣)

### 3.1. Explore snake rules that `magpipe_snake.smk` in the `snakes` folder refers to.

❗Note! The `magpipe` pipeline starts from the depth.smk rule. This step is described in the original paper as:

> The resulting BAM files were processed using the jgi_summarize_bam_contig_depths script of MetaBAT2 
> (v.2.12.1)[55](https://www.nature.com/articles/s41586-022-04862-3#ref-CR55) to provide within- and 
> between-sample coverages for each scaffold.

### 3.2. Setting up configs:

The [configs](https://github.com/SushiLab/magpipe/tree/master/configs) folder contains the `test_config.yaml` file, 
that direct snakemakes to the desired step of the pipeline and also will give snakemake paths to environment files and useful directories.

- **targets** - choose targets for the magpipe pipeline by uncommenting chosen branch. Note, that you cannot run several phases at once, you need to choose one at a time. Also, be careful with indentation, if you mess up the pipeline won’t work :( ;
- **dataset_info** - a path to the file in the resources folder that contains information about the dataset;
- **external_genomes** - a path to the external genomes that will be added to the pipeline at the evaluation step (`evaluate.smk`);
- **output, envs, snake_workdir** - add these paths according to descriptions provided in the `test_config.yaml` file.


### 3.3. Setup resources

The `test_dataset.tsv` file is a table with columns:

1. Datasets: a name of the dataset;
2. Samples: a path to the list of samples == a path to the `samles.list` file, which is a text file with names of files with the filtered reads;
3. Debug: a single sample name for debugging;
4. Assemblies: a path to the directory with assemblies (formatted as `{sample}_scaffolds.min1000.fasta.gz`). Meaning for each sample there should be an assembly; 
    > Assemblies acquisition: see the [phase 0](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/phase_0)
    > on how to prepare suitable files.
5. Backmapping: a path to the bamfiles directory (called `path/*.bam`);
    > BAM files acquisition: see the [phase 0](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/phase_0)
    > on how to prepare suitable files.
6. diffcov: the number of samples mapped back to each assembly.


The alternative way how to prepare the same files from the same project can be found 
in the [Tara Ocean](https://merenlab.org/data/tara-oceans-mags/) project. Moreover, there are a lot of links with preprocessed files 
such as already-assembled scaffolds.

## 4. Preparations for running snakemake on cluster mode (SLURM) ☆ﾐ(o*･ω･)ﾉ

This step is based on the SchlossLab [tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial/blob/master/README.md).

Create the following files:
- `cluster.yaml`:
```yaml
__default__:
  jobname: "{rule}.{wildcards}"
  nodes: 1
  ntaskspernode: 1
  cpuspertask: 1
  mempercpu: "1000mb"
  time: "5:00"
  account: "YOUR_ACCOUNT"
  partition: "shared"
  output: "path/to/magpipe/snakes/cluster_config/snake_output.out"

count_words:
  cpuspertask: 2
  time: "10:00"
```
- `config.yaml`
```yaml
cluster-config: "path/to/magpipe/snakes/cluster_config/cluster.yaml"cluster: "sbatch --job-name={cluster.jobname} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntaskspernode} --cpus-per-task={cluster.cpuspertask} --mem-per-cpu={cluster.mempercpu} --time={cluster.time} --account={cluster.account} --partition={cluster.partition} --output={cluster.output} --export=ALL "jobs: 99verbose: truenotemp: true
```
- `run_snakemake.sh`
```bash
#!/bin/bash -l
# The -l above is required to get the full environment with modules

#Set the allocation to be charged for this job
#SBATCH -A <PROJ_NAME>

# The name of the script
#SBATCH -J snakemake

# The partition
#SBATCH -p shared

# 48 hours wall-clock time will be given to this job
#SBATCH -t 00:05:00

# Number of Nodes 
#SBATCH --nodes=1

# Number of MPI tasks per node
#SBATCH --ntasks-per-node=1
# Number of logical cores hosting OpenMP threads. Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=200


# Run the executable named my exe
snakemake -s magpipe_snake.smk --configfile path/to/magpipe/configs/test_config.yaml --use-conda --profile path/to/magpipe/snakes/cluster_config/ --latency-wait 20
```


## 5. Run snakemake and debug ＵＴｪＴＵ
```bash
sbatch run_snakemake.sh
```


Good luck! :shipit:




