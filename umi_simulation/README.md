## Simulated UMI grouping evaluation 

This folder contains scripts for UMI grouping evaluation using simulated UMI and the extended result for simulated UMI grouping evaluations from CELLO-seq manuscript. 
For more detailed explanation of the simulated UMI grouping evaluation, please check out CELLO-seq manuscript.

The workflow for the simulated UMI grouping evaluation is as follow:

1. Simulate UMI using the helper scripts provided depending on the type of read:
   * Perfect reads (read without any sequence error and with Illumina high qualirty scores (Q40))
       1. Run `generate_perfect_simulated_umi_fastq.sh` to generate a FASTQ file containing perfect simulated UMI sequences.
   * ONT reads (read with current ONT read identity and sequence error rate)
       1. Run `generate_badread_fasta_umi_input.sh` to generate a FASTA file containing simulated UMI sequences.
       2. Run `generate_badread_simulated_umi_fastq.sh` with the FASTA file generated from previous step to generate a FASTQ file containing ONT simulated UMI sequence.
2. Run `create_umi_groups.R` to produce UMI groupings using the `expectedDist` function from R Sarlacc package [Github](https://github.com/MarioniLab/Sarlacc). 
3. Run `evaluate_umi_groupings_purity.R` to evaluate the purity of the UMI groups formed by the previous step.
