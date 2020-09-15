#!/bin/sh

# This script is designed for evaluating alignment result of Sarlacc-processed simulated reads against their reference location.

full_reads=/path/to/organism/full_reads_list.txt
sam_files=$(ls /path/to/unprocessed/aligned/*.genome.sam)

for sam in ${sam_files}; do
  python3 evaluate_processed_aligned_reads.py -f ${full_reads} $sam ${sam%.genome.sam}.readnames.1.txt ${sam%.genome.sam}
done