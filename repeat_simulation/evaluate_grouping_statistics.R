# This script is designed to report statistics on the UMI groups formed by CELLO-seq Sarlacc pipeline.

readnames_files = list.files(pattern="L1.*readnames\\.1\\.rds")

readnames_grouping_statistics = Reduce(rbind, lapply(readnames_files, function(readname_file) {
  organism = gsub("_simulated.*", "", readname_file)
  simulated_read_length = gsub("(.*pattern_|\\.fasta.*)", "", readname_file)
  
  readnames_grouping = readRDS(readname_file)
  readnames_grouping = lapply(readnames_grouping, function (x) gsub("_[0-9]+?bp_.*$", "", x))
  
  readnames_groupings_counts = lapply(readnames_grouping, table)
  
  data.frame(organism=organism,
             simulated_read_length=simulated_read_length,
             group_size=lengths(readnames_grouping), 
             unique_readnames_count=lengths(readnames_groupings_counts),
             largest_readnames_count=unlist(lapply(readnames_groupings_counts, max)), 
             unique=ifelse(lengths(readnames_groupings_counts)==1, TRUE, FALSE)
   )
}))

saveRDS(readnames_grouping_statistics, "simulated_readnames_grouping_statistics.rds")
