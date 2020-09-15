# This script is designed to extract the read name and size of each UMI groups formed by CELLO-seq Sarlacc pipeline.

rds_files <- list.files(pattern="L1[M_].*readnames\\.1\\.rds$")

for (rds in rds_files) {
  r = readRDS(rds)
  # If the read groups are pure, we will extract the read name, otherwise we tag it as mixed group
  names = gsub(",.*", "",
               sapply(r, function(g)
                 ifelse(length(table(gsub(",.*", "", g))) == 1, g[1], "MIXED")
               )
  )
  write.table(names, file=gsub(".rds", ".txt", rds), col.names=FALSE, row.names=FALSE, quote=FALSE)

  len = lengths(r)
  write.table(data.frame(read_name=names, sizes=len), file=gsub(".rds", "_sizes.txt", i),
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}

