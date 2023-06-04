#!/usr/bin/env python

"""
name: phase_0
description: This snakemake pipeline is part of the MAGs building pipeline called magpipe. Phase 0
dependency: phase0_requirements.txt
author: Polina Guseva and Anna Zhurakovskaya
rules:
    - prepare_for_one_sample 	# subsample and upload reads only for our sample
    - clean_reads		# clean from adapters, control and low-quality reads
    - norm_reads		# normalise reads for further assembly
    - assembly_reads
    - mapping_reads		# map reads to the assembly and sort only good maps
"""


rule prepare_for_one_sample:
    input:
        "samples.txt"
    output:
        "{sample_name}_samples.txt"
    params:
        sample_name="sample_name"
    shell:
        """
        # Upload all dependencies
        # conda create --name phase_0 --file phase0_requirements.txt
        
        echo "Starting of read downloading"

        python load_raw_reads.py {sample_name}

        echo "Reads are downloaded and all preparation for the new step are made"
        """


rule clean_reads:
    input:
        "{sample_name}_samples.txt",
        control = "phix.fa"
    output:
        "{sample_name}_R1_phix_removed.fastq",
        "{sample_name}_R2_phix_removed.fastq"
    shell:
        """
        echo "Starting cutting adapters"

        # generate a configuration file for processing sequencing data
              iu-gen-configs {wildcards.sample_name}_samples.txt

        # filter the quality of sequencing reads based on it
        iu-filter-quality-minoche {wildcards.sample_name}.ini --ignore-deflines

        echo "Adapters are cut"


        echo "Starting trimming low-quality reads"

        # remove low-quality reads
        bbduk.sh in1={wildcards.sample_name}-QUALITY_PASSED_R1.fastq in2={wildcards.sample_name}-QUALITY_PASSED_R2.fastq out1={wildcards.sample_name}_filtered_R1.fastq out2={wildcards.sample_name}_filtered_R2.fastq qtrim=rl trimq=14 maq=20 maxns=0 minlength=45

        # removing a control
        bbmap.sh in={wildcards.sample_name}_filtered_1.fastq out={wildcards.sample_name}_R1_phix_removed.fastq ref={control} nodisk
        bbmap.sh in={wildcards.sample_name}_filtered_2.fastq out={wildcards.sample_name}_R2_phix_removed.fastq ref={control} nodisk

        echo "Low-quality reads are trimmed"

        """


rule norm_reads:
    input:
        "{sample_name}_R1_phix_removed.fastq",
        "{sample_name}_R2_phix_removed.fastq"
    output:
        "{sample_name}_R1_normalized.fastq",
        "{sample_name}_R2_normalized.fastq"
    shell:
        """
        echo "Starting normalisation"

        # normalise reads
        bbnorm.sh in={wildcards.sample_name}_R1_phix_removed.fastq out={wildcards.sample_name}_R1_normalized.fastq target=40 mindepth=0
        bbnorm.sh in={wildcards.sample_name}_R2_phix_removed.fastq out={wildcards.sample_name}_R2_normalized.fastq target=40 mindepth=0

        echo "Normalisation is done"
        """


rule assembly_reads:
    input:
        "{sample_name}_R1_normalized.fastq",
        "{sample_name}_R2_normalized.fastq"
    output:
        "{sample_name}_contigs.fa"
    shell:
        """
        echo "Starting assembly. It is a TIME-CONSUMING process"
        mkdir {wildcards.sample_name}_assembly

        # assemble reads into scaffolds 
        megahit -1 {wildcards.sample_name}_R1_normalized.fastq -2 {wildcards.sample_name}_R2_normalized.fastq --min-contig-len 1000 -o {wildcards.sample_name}_assembly --num-cpu-threads 80

        # retrieve a file with scaffolds
        mv {wildcards.sample_name}_assembly/final.contigs.fa .
        mv final.contigs.fa {wildcards.sample_name}_contigs.fa

        echo "Assembly is done"
        """


rule mapping_reads:
    input:
        first_read = "{sample_name}_R1_phix_removed.fastq",
        second_read = "{sample_name}_R2_phix_removed.fastq",
        scaffold = "{sample_name}_contigs.fa"
    output:
        "{sample_name}_mapped_sorted.bam"
    shell:
        """
        echo "Starting mapping"

        # Index the reference
        bwa index {scaffold}

        # Map reads onto the reference
        bwa mem -a -o {wildcards.sample_name}_mapped.sam {scaffold} {first_read} {second_read}
        echo "Mapping is finished"

        echo "Starting sorting bam files"

        # Transform SAM into BAM
        samtools view -F 4 -bS -o {wildcards.sample_name}_mapped.bam {wildcards.sample_name}_mapped.sam

        # Remove reads with insufficient mapping
        samtools view -h {wildcards.sample_name}_mapped.bam | awk -F'\t' 'BEGIN {{OFS="\t"}} {{ if ($1 ~ /^@/ || ($6 ~ /^[0-9]+M$/ && $5 >= 45 && (($4 + $5) / length($10)) >= 0.8)) print }}' | samtools view -bS - > {wildcards.sample_name}_mapped_filtered.bam

        # Sort the final BAM file
        samtools sort -o {output} {wildcards.sample_name}_mapped_filtered.bam

        echo "Sorting of bam files is done. Phase 0 is finished. Please find all necessary files in the current directory"
        """
