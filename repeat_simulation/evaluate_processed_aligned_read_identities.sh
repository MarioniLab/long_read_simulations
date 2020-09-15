#!/bin/sh

# This script is designed for evaluating the read identity of Sarlacc-processed simulated reads against
# their reference sequence based on the sequence alignment.

sam_files=$(ls /path/to/unprocessed/aligned/*.genome.sam)

for sam in ${sam_files}; do
  python3 evaluate_processed_aligned_read_identities.py $sam ${sam%.genome.sam}.readnames.1.txt ${sam%.genome.sam}_sequence_identity
done
