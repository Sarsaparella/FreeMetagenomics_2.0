#!/bin/bash

###########################################################################################################
# Script Name: test_script_phase_0.sh
# Description: This script performs a test run for the ASE_68_05M sample from the Tara Ocean project
# Author: Polina Guseva
# Created: April 18, 2023
# Dependencies: illumina-utils, BBMap, MEGAHIT, bwa, samtools
###########################################################################################################

# download raw read
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/ERR598982/ERR598982_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599107/ERR599107_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/ERR598982/ERR598982_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599107/ERR599107_2.fastq.gz

# generate a configuration file for processing sequencing data
iu-gen-configs samples_test.txt
# filter the quality of sequencing reads based on it
iu-filter-quality-minoche ASE_68_05M.ini --ignore-deflines

# remove low-quality reads
./bbmap/bbduk.sh in1=ASE_68_05M-QUALITY_PASSED_R1.fastq in2=ASE_68_05M-QUALITY_PASSED_R2.fastq out1=ASE_68_05M_filtered_R1.fastq out2=ASE_68_05M_filtered_R2.fastq qtrim=rl trimq=14 maq=20 maxns=0 minlength=45
# removing a control
./bbmap/bbmap.sh in=ASE_68_05M_filtered_1.fastq out=ASE_68_05M_R1_phix_removed.fastq ref=phix.fa nodisk
./bbmap/bbmap.sh in=ASE_68_05M_filtered_2.fastq out=ASE_68_05M_R2_phix_removed.fastq ref=phix.fa nodisk

# merge forward and reverse reads
./bbmap/bbmerge.sh in1=ASE_68_05M_R1_phix_removed.fastq in2=ASE_68_05M_R2_phix_removed.fastq out=ASE_68_05M_merged.fastq outu=ASE_68_05M_unmerged.fastq ihist=histogram.txt minoverlap=16

# normalise reads
./bbmap/bbnorm.sh in=ASE_68_05M_R1_phix_removed.fastq out=ASE_68_05M_R1_normalized.fastq target=40 mindepth=0
./bbmap/bbnorm.sh in=ASE_68_05M_R2_phix_removed.fastq out=ASE_68_05M_R2_normalized.fastq target=40 mindepth=0
./bbmap/bbnorm.sh in=ASE_68_05M_merged.fastq out=ASE_68_05M_merged_normalized.fastq target=40 mindepth=0
./bbmap/bbnorm.sh in=ASE_68_05M_unmerged.fastq out=ASE_68_05M_unmerged_normalized.fastq target=40 mindepth=0

# assemble reads into scaffolds 
megahit -1 ASE_68_05M_R1_normalized.fastq  -2 ASE_68_05M_R2_normalized.fastq --min-contig-len 1000 -o assembly_folder --num-cpu-threads 80

# retrieve a file with scaffolds
mv assembly_folder/final.contigs.fa .

# index the reference
bwa index final.contigs.fa

# map reads onto the reference (after using megahit)
bwa mem -a -o ASE_68_05M_mapped.sam final.contigs.fa ASE_68_05M-QUALITY_PASSED_R1.fastq ASE_68_05M_R2-QUALITY_PASSED_R2.fastq

# transform SAM into BAM
samtools view -F 4 -bS -o ASE_68_05M_mapped.bam ASE_68_05M_mapped.sam
# remove reads with insufficient mapping
samtools view -h ASE_68_05M_mapped.bam | awk -F '\t' 'length($10) >= 45 && $3 >= 97 && ($4/length($10)) >= 0.8' | samtools view -bS - > ASE_68_05M_mapped_filtered.bam
# sort the final BAM file
samtools sort -o ASE_68_05M_mapped_sorted.bam ASE_68_05M_mapped_filtered.bam
