#!/usr/bin/env python3

"""
Script Name: load_raw_reads.py
Description: This script uploads all raw reads for a provided sample name and prepares a *.csv file with names of raw reads
Author: Polina Guseva
Date: 04.06.2023

Dependencies:
    - pandas==1.5.3
"""

import os
import pandas as pd
import sys

sample_name = sys.argv[1]

# Read the sample list file
sample_df = pd.read_table("samples.txt")

# Filter the row based on the sample name
subsample_df = sample_df[sample_df['sample'] == sample_name]

# Write the DataFrame to the output TXT file
subsample_df.to_csv("{output}", sep="\t", index=False)

# Upload reads only for our sample
base_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
for name in ''.join(subsample_df.r1).split(',') + ''.join(subsample_df.r2).split(','):
	link = os.path.join(base_url, name[:6], name.split('_')[0], name)
  os.system(f"wget {link}")
  print(f"\n done: wget {link} for the {sample} sample \n")
  
