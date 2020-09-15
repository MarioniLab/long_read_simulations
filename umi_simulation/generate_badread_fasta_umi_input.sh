#!/bin/sh

# This script is designed for generating simulated UMI in FASTA format as input for Badread simulator.

umi_count=20000
umi_length=10
fasta_file_prefix=/path/to/simulated_umi

python3 generate_simulated_umi_reads.py --output_type fasta --umi_lengths ${umi_length} $umi_count $fasta_file_prefix
