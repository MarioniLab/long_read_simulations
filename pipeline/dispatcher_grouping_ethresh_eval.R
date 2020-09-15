template.file <- "grouping_simulation.Rmd"
template.content <- readLines(template.file)

sample.list <- list.files(pattern="\\.rds$")
for (ethresh in c(3,5,7,9,11)) {
  for (sample in sample.list) {
    sample.prefix <- deparse(gsub(".rds", "", sample))
    script <- gsub("INDEX", sample.prefix, template.content)
    script <- gsub("THRESHOLD", ethresh, script)
    writeLines(text = script, con = paste0("grouping_", sample.prefix, ".Rmd"))
  }
}