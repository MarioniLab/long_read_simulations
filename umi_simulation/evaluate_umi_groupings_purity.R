library(sarlacc)
library(Biostrings)

# This script is designed to evaluate the purity of the UMI groups formed by create_umi_groups.R script

rds_file = commandArgs()
rds_file = rds_file[length(rds_file)]
rds_file_name = gsub(".rds", "", gsub(".*/", "", rds_file))

fastq_file = gsub(".rds", ".fastq", gsub("group_[^_]+_edist_", "", rds_file))
simulated_reads = readQualityScaledDNAStringSet(fastq_file)

# We need to fix the strand of the reads for badread simulated read
if (!grepl("perfect", rds_file_name)){
  reverse_reads_index =  which(grepl("-strand", names(simulated_reads)) == TRUE)
  reverse_reads = reverseComplement(simulated_reads[reverse_reads_index])
  
  corrected_reads = c(simulated_reads[-reverse_reads_index], reverse_reads)
  names(corrected_reads) = gsub(",.*", "", names(corrected_reads))
} else {
  corrected_reads = simulated_reads
  names(corrected_reads) = gsub(",.*", "", names(corrected_reads))
}

group_pregroups = readRDS(rds_file)

thresholds = c(2,4,6,8,10,12,14,16)

group_pregroups_mixed_counts = lapply(1:length(group_pregroups), function(group_thresh){
  group_composition_count = unlist(lapply(group_pregroups[[group_thresh]], function(x) length(table(names(corrected_reads[x])))))
  group_composition_count_summary = table(group_composition_count)
  
  pure_count = ifelse(names(group_composition_count_summary)[1] == 1, group_composition_count_summary[1], 0)
  mixed_count = length(group_composition_count) - pure_count
  
  summary_stat = data.frame(group_thresh=thresholds[group_thresh], stat=c("pure", "mixed"), counts=c(pure_count, mixed_count))
})
group_pregroups_mixed_counts = Reduce(rbind, group_pregroups_mixed_counts)

saveRDS(group_pregroups_mixed_counts, paste0(rds_file_name, "_mixed_count_stats.rds"))

