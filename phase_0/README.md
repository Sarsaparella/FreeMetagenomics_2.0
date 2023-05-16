# Phase 0: from raw reads to scaffolds and mapped reads

This part of the pipeline is missing in [the original repository](https://github.com/SushiLab/magpipe). 
> However, there is a description of this in the article ([Paoli _et al._ 2022](https://www.nature.com/articles/s41586-022-04862-3#Sec7)). It would be provided for each part.

So far, this phase could be split up into several steps:
1. removing adapters, a control as of the PhiX genome, and low-quality reads;
2. merging forward and reverse reads;
3. normalising obtained files (forward, reverse, merged, and unmerged reads)
4. performing an assembly of scaffolds;
5. mapping filtered reads to scaffolds.


## Step 1: removing the trash
> Sequencing reads from all metagenomes were quality filtered using BBMap (v.38.71) by removing sequencing adapters from the reads, 
> removing reads that mapped to quality control sequences (PhiX genome) and discarding low quality reads using the parameters trimq = 14, maq = 20, maxns = 0 and minlength = 45.

Following the original description, the first step could be performed as:
```ruby
# Quality filtering using BBMap
bbduk.sh in=input.fastq out=filtered.fastq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=14 maq=20 maxns=0 minlength=45
# Removing PhiX genome reads
bbduk.sh in=filtered.fastq out=filtered_phix_removed.fastq ref=phix_genome.fa k=31 hdist=1 stats=stats.txt
```

However, this requeres `adapters.fa` which is not provided. Therefore, it was replaced by:
```ruby
# Removing adapters using illumina-utils
pip install illumina-utils
iu-gen-configs samples.txt
iu-filter-quality-minoche sample_name.ini --ignore-deflines
```
Around 1 hour.

The `samples.txt` TAB-separated file from the TARA OCEAN project could be found by the following [link](http://merenlab.org/data/tara-oceans-mags/files/samples.txt). It looks like:
```
sample	r1	r2
sample_name	reads1_1.fastq.gz,reads2_1.fastq.gz	reads1_2.fastq.gz,reads2_2.fastq.gz
```

Afterwards, the clean reads still should be filtered:
```ruby
bbduk.sh in1=sample_name-QUALITY_PASSED_R1.fastq in2=sample_name-QUALITY_PASSED_R2.fastq out1=sample_name_filtered_R1.fastq out2=sample_name_filtered_R2.fastq qtrim=rl trimq=14 maq=20 maxns=0 minlength=45
```
Around 2.5 min.


And the control should be removed too:
```ruby
bbmap.sh in=sample_name_filtered_*.fastq out=*_phix_removed.fastq ref=phix.fa nodisk
```
Around 2 min.

## Step 2: merging
> Downstream analyses were performed on quality-controlled reads or if specified, merged quality-controlled reads (bbmerge.sh minoverlap = 16).
```ruby
bbmerge.sh in1=R1_phix_removed.fastq in2=R2_phix_removed.fastq out=merged.fastq outu=unmerged.fastq ihist=histogram.txt minoverlap=16
```
Around 4 min.

## Step 3: normalising
> Quality-controlled reads were normalized (bbnorm.sh target = 40, mindepth = 0)

All files should be normalised (forward, reverse, meerged, and unmerged). The code example is for `merged.fastq`:
```ruby
bbnorm.sh in=merged.fq out=merged_normalized.fq target=40 mindepth=0
```
Amound 20 min for 30Gb. Time increases linearly.

## Step 4: shit of shit
> ... they were assembled with metaSPAdes (v.3.11.1 or v.3.12 if required). The resulting scaffolded contigs (hereafter scaffolds) were finally filtered by length (≥1 kb).

## Step 5: mapping
> ... metagenomic samples were grouped into several sets and, for each sample set, the quality-controlled metagenomic reads from all samples were individually mapped against the scaffolds of each sample, resulting in the following numbers of pairwise readset to scaffold mappings
> Mapping was performed using the Burrows–Wheeler-Aligner (BWA) (v.0.7.17-r1188)54, allowing the reads to map at secondary sites (with the -a flag). 
> Alignments were filtered to be at least 45 bases in length, with an identity of ≥97% and covering ≥80% of the read sequence. 

```ruby
bwa index scaffold.fa
bwa mem -a -o mapped.sam scaffold.fa forward_normalized.fastq reverse_normalized.fastq
samtools view -bSh -F 4 -q 20 -f 0x2 -o mapped.bam mapped.sam
samtools sort -o mapped_sorted.bam mapped.bam
```
