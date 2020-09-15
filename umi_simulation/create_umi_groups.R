library(sarlacc)
library(Biostrings)

# This script is designed to produce UMI groupings from simulated UMI with and without pregroupings

fastq_file = commandArgs()
fastq_file = fastq_file[length(fastq_file)]
fastq_file_name = gsub(".fastq", "", gsub(".*/", "", fastq_file))
simulated_reads = readQualityScaledDNAStringSet(fastq_file)

# We need to fix the strand of the reads for badread simulated read
if (!grepl("perfect", fastq_file_name)){
  reverse_reads_index =  which(grepl("-strand", names(simulated_reads)) == TRUE)
  reverse_reads = reverseComplement(simulated_reads[reverse_reads_index])

  corrected_reads = c(simulated_reads[-reverse_reads_index], reverse_reads)
  names(corrected_reads) = gsub(",.*", "", names(corrected_reads))
} else {
  corrected_reads = simulated_reads
  names(corrected_reads) = gsub(",.*", "", names(corrected_reads))
}

group_size = 100
thresholds = c(2,4,6,8,10,12,14,16)

# We will first run umiGroup with pre-grouping
read_names = gsub(",.*", "", names(simulated_reads))
unique_read_names = unique(read_names)
groups = as.character(rep.int(c(1:(length(unique_read_names)/group_size)), group_size))
pre_groups = groups[match(names(corrected_reads), unique_read_names)]
groupings_pregroup = lapply(thresholds, function(t){
  umiGroup(corrected_reads, threshold1 = t, groups = pre_groups)
})
saveRDS(groupings_pregroup, paste0("group_pregroup_edist_", fastq_file_name, ".rds"))

# Now we will run umiGroup without any pre-grouping
groupings_all = lapply(thresholds, function(t){
  umiGroup(corrected_reads, threshold1 = t) 
})
saveRDS(groupings_all, paste0("group_all_edist_", fastq_file_name, ".rds"))


