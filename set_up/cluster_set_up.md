# How to create a cosy work-space on PDC

Here you can find a detailed guideline on how to install all modules and environments to start working 
with the [magpipe](https://github.com/SushiLab/magpipe/tree/master) pipeline by SushiLab on PDC.


## 1. Upload all necessary modules 
For example, if you work on PDC, you will 100% need `conda` and a text editor.

```bash
# if you working on PDC
module load PDC
module load anaconda3/2022.05
module load vim/8.2

# Useful commands for operating in a module system
module spider <module_name> # search for modules
module avail # show available modules
module list # show a list of all your modules
```

## 2. Initialise conda
If you are using `conda` for the first time on a cluster.

```bash
conda init
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set auto_activate_base false
```

## 3. Download and install mamba
Also, to optimise a work of `Snakemake` with several environments, 
[installation of `Mamba`](https://mamba.readthedocs.io/en/latest/installation.html) is recommended too.

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```
Additionally, it is recommended to check where are the folders for installed packages and environments.
If they are stored in a place where you can potentially run out of memory, please, change their location then.


Here are some useful commands you can use to check that:

```bash
conda config --show
conda/mamba info
conda config --add pkgs_dirs /new_path_pkgs/
```
## Create environment for phase 0

The original magpipe lacks the first step that prepairs input files for the pipeline. To prepares those files you will need to install additional software in a separate environment. 

You can find the file with requirements for this environment in the “phase_0” folder

```bash
conda create --name phase_0 --file phase0_requirements.txt
```

> Please, be aware that this is an example for the PDC cluster, however, we hope it is universal :)







