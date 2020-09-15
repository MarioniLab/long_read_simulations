## Long-read sequencing simulation evaluation of repeat elements alignment 

This folder contains scripts for evaluation of repeat elements alignment with reads generated using CELLO-seq protocol.
This folder also contains the extended result for L1 elements simulation evaluation from the CELLO-seq manuscript. 
For more detailed explanation of the simulated repeat elements evaluation, please check out CELLO-seq manuscript.

The workflow for the simulated repeat elements evaluation depends on the type of read:

1. Simulate reads from repeat elements using the helper scripts provided depending on the type of read:
   * Perfect reads (read without any sequence error and with Illumina high qualirty scores (Q40))
       1. Run `generate_perfect_simulated_fastq.sh` to generate a FASTQ file containing perfect simulated reads from repeat elements. 
   * ONT reads (read with current ONT read identity and sequence error rate)
       1. Run `generate_badread_fasta_input.sh` to generate a FASTA file containing simulated reads from repeat elements.
       2. Run `generate_badread_simulated_fastq.sh` with the FASTA file generated from previous step to generate a FASTQ file containing ONT simulated reads from repeat elements.
   * CELLO-seq processed reads (Deduplicated or error corrected ONT reads processed with CELLO-seq Sarlacc pipeline)
       1. Run `generate_badread_fasta_input.sh` to generate a FASTA file containing simulated reads from repeat elements.
       2. Run `generate_badread_simulated_fastq.sh` with the FASTA file generated from previous step to generate a FASTQ file containing ONT simulated reads from repeat elements.
       3. Run CELLO-seq Sarlacc data processing pipeline ([Github](https://github.com/MarioniLab/CELLOseq), along with modified scripts from the pipeline folder of this repository) to produce either deduplicated or error corrected FASTQ file.
2. Run `align_simulated_reads.sh` to align the simulated reads from repeat elements back to the FASTA containing all the repeat elements used for simulation. 
3. Evaluate the result of the simulated repeat element read alignment using the helper scripts provided depending on the type of read:
   * Perfect reads and ONT reads
       1. Run `evaluate_unprocessed_aligned_reads.sh` to evaluate the alignment result of unprocessed simulated reads against their reference location.
   * CELLO-seq processed reads 
       1. Run `extract_readnames_groupsizes.R` to extract the read name (containing the true reference location) and size for each UMI groups formed by Sarlacc pipeline.
       2. Run `evaluate_grouping_statistics.R` to generate statistics on the UMI groups formed by Sarlacc pipeline.
       1. Run `evaluate_processed_aligned_reads.sh` to evaluate the alignment result of Sarlacc-processed simulated reads against their reference location.
       2. Run `evaluate_processed_aligned_read_identities.sh` to evaluate the read identity of Sarlacc-processed simulated reads against their reference sequence based on the sequence alignment.
