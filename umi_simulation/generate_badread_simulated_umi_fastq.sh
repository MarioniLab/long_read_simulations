#!/bin/bash

# This script is designed for simulation of ONT read using Badread. The identity parameters used for simulation are:
# - ONT - 92,96,2.5
# - ONT(96%) - 96,99.9,2.5

fasta_files=$(ls /path/to/simulated/*_umi.fasta)
quantity=(1x 5x 10x 15x)
trained_error_model=/path/to/trained_error_model
trained_qscore_model=/path/to/trained_qscore_model

for fasta in ${fasta_files}; do
  for q in ${quantity[@]}; do
    fastq_file=${fasta%.fasta}_${q}cov.fastq

    badread simulate \
    --reference $fasta \
    --quantity $q \
    --identity 92,96,2.5 \
    --length 100000,100 \
    --error_model ${trained_error_model} \
    --qscore_model ${trained_qscore_model} \
    --start_adapter_seq '' --end_adapter_seq '' \
    --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 > $fastq_file
  done
done

