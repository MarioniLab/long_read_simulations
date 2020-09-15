#!/bin/sh

# This script is designed for generating perfect simulated reads in FASTQ format.

fasta_file=/path/to/repeat.fasta
fastq_file_prefix=/path/to/simulated/repeat

python3 generate_simulated_reads.py $fasta_file $fastq_file_prefix