template.file <- "adaptor_filtering_simulation.Rmd"
template.content <- readLines(template.file)

sample.list <- list.files(pattern="L1_.*\\.fastq$")
for (sample in sample.list) {
  script <- gsub("FASTQ_FILE", deparse(sample), template.content)
  writeLines(text = script, con = paste0("demultiplex_simulation_", sample, ".Rmd"))
}