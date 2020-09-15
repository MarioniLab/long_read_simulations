#!/bin/bash

# This script is designed for alignment of simulated reads with minimap2

species=human
fastq_list=$(ls /path/to/fastq/folder/*.fastq)

nproc=5

reference_folder=/path/to/reference/folder
index_type=(repeat)
if [ "$species" == "human" ]
then
  file_pattern=L1_GRCh38_
  index_files=(GRCh38_repeats_L1_2kbplus_unique.fa.mmi)
elif [ "$species" == "mouse" ]
then
  file_pattern=L1Md_mm10_
  index_files=(L1Md_mm10_ucsc_2kbplus_unique.fa.mmi)
fi

echo "Species set to ${species}"
echo "Species folder set to ${species_folder} and reference files set to ${index_files}"
echo "FASTQ files to be processed are ${fastq_list}"

for fastq in ${fastq_list}; do
  for i in "${!index_files[@]}"; do
    minimap2 -K 200M --eqx --MD -t ${nproc} -ax splice -uf --secondary=no -C5 \
    ${reference_folder}/${species_folder}/${index_files[$i]} ${fastq} > ${fastq%fastq}.${index_type[$i]}.sam
  done
done

