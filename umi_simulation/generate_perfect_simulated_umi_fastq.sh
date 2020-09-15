#!/bin/sh

# This script is designed for generating perfect simulated UMI in FASTQ format.

umi_count=10000
umi_length=10
fastq_file_prefix=/path/to/simulated_umi

python3 generate_simulated_umi_reads.py --umi_lengths ${umi_length} $umi_count $fastq_file_prefix
