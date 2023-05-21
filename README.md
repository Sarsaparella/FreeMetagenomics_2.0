# Spatial mapping of potential enzymes encoded in the global ocean microbiome 
by Polina Guseva and Anna Zhurakovskaya


This repository contains a detailed explanation of how to use the [migpipe](https://github.com/SushiLab/magpipe) pipeline by SushiLab described in [Paoli _et al._ (2022)](https://www.nature.com/articles/s41586-022-04862-3).


The aim is to restore the whole pipeline on the PDC cluster to work from raw reads to BGCs (biosynthetic gene clusters) and annotated protein-coding genes of High Seas. The samples are taken from the published [database](https://www.microbiomics.io/ocean/supp_info/) from the original project.


1. To start working, set up your workspace (go to [./set_up/cluster_set_up.md](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/set_up/cluster_set_up.md)) and follow each step. If necessary tailor specifically to your cluster. 
3. _(Additional)_ If necessary sort your samples into maritime zones by the script from the [maritime_zonation](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/maritime_zonation) folder and extract only samples from High Seas. It is recommended first to have a look in the notebook and modify the code to your specific spreadsheet.
4. After that, transform raw reads to scaffolds and mapped reads as described in the [phase_0](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/phase_0) folder. Please, be aware, that the scaffold assembly step is an extremely memory-intensive process (:exclamation:). Even to run one sample through it is required ~36-48 hours with >=512 Gb RAM. Therefore, already-assembled scaffolds could be used.
5. Later, clone the original pipeline to your workspace and install all dependencies as described in [./set_up/magpipe_set_up.md](https://github.com/GusevaPolina/FreeMetagenomics/tree/main/set_up/magpipe_set_up.md). With completed phase 0, a path to scaffolds and mapped reads could be inserted into the `test_dataset.tsv` file (see step 3.3 of Setting up resources).
6. Run the main snakemake script as mentioned in the last step of magpipe_set_up.md. Debug each problem if necessary.


As a result, a list of BGCs and annotated protein-coding genes of prokaryotic organisms is provided. As an example of possible output from the Paoli _et al._ (2022) contains a lot of antibiotics, drugs, nutraceuticals, biofuels, etc. Among them, many new bioactive compounds and pathways were discovered.


Good luck! :shipit:
