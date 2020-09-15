#!/bin/sh

# This script is designed for evaluating alignment result of unprocessed simulated reads against their reference location.

sam_files=$(ls /path/to/unprocessed/aligned/*.genome.sam)

for sam in ${sam_files}; do
  python3 evaluate_unprocessed_aligned_reads.py $sam ${sam%.genome.sam}
done
