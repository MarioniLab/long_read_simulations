#!/bin/sh

# This script is designed for generating simulated reads in FASTA format as input for Badread simulator.

fasta_file=/path/to/repeat.fasta
fasta_output_file_prefix=/path/to/simulated/repeat
umi_length=22

python3 generate_simulated_reads.py --output_type fasta --add_adapter --umi_length ${umi_length} \
$fasta_file $fasta_output_file_prefix
