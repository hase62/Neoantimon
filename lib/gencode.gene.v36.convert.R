data <- scan("gencode.v36.annotation.gtf", "character", sep = "\n", skip = 6)

res <- matrix(nrow = length(data), ncol = 3, NA)
for(i in 1:length(data)){
  sep_i <- strsplit(gsub('"', '', data[i]), "\t|;| ")[[1]]
  chr_position <- paste(gsub("chr", "", sep_i[1]), ":", sep_i[4], "-", sep_i[5], sep = "")
  gene_name <- sep_i[match("gene_name", sep_i) + 1]
  gene_id <- sep_i[match("gene_id", sep_i) + 1]
  res[i, ] <- c(gene_name, gene_id, chr_position)
}

write.table(res, "gencode.gene.info.v36.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

