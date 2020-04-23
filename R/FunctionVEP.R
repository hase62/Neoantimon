convert_to_annovar_format_from_vep <- function(vep_file) {
  #Read vep data
  if(is.list(vep_file) | is.matrix(vep_file)){
    data <- as.matrix(vep_file)
  } else {
    data <- read_data(vep_file)
  }
  data <- data[grep("missense_variant|insertion|deletion|frameshift", data[, match("Consequence", colnames(data))]), ]

  #Execute ensembl
  print("Executing Transformation")
  ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
  values <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", values = data[, match("Gene", colnames(data))], mart= ensembl)

  #Make annovar format data
  index <- c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
  data_an <- matrix(ncol = length(index), nrow = nrow(data), ".")
  colnames(data_an) <- index
  data_an[, match(c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"), colnames(data_an))] <-
          t(apply2(data, 1, function(x) {
            x <- as.character(x)
            tmp <- strsplit(x[match("Location", colnames(data))], ":")[[1]]
            Chr <- tmp[1]
            Start <- tmp[2]
            End <- tmp[2]
            Ref_acid <- gsub("a|t|g|c", "", strsplit(x[match("Codons", colnames(data))], "/")[[1]][1])
            Ref <- ifelse(length(grep("STRAND=1", x[match("Extra", colnames(data))])) == 1, Ref_acid, trans_to[match(tolower(Ref_acid), trans_from)])
            Alt_acid <- gsub("a|t|g|c", "", strsplit(x[match("Codons", colnames(data))], "/")[[1]][2])
            Alt <- ifelse(length(grep("STRAND=1", x[match("Extra", colnames(data))])) == 1, Alt_acid, trans_to[match(tolower(Alt_acid), trans_from)])
            Func.refGene <- "exonic"
            Gene.refGene <- paste(unique(values[which(!is.na(match(values[, 2], x[match("Gene", colnames(data))]))), 3]), collapse = ";")
            ExonicFunc.refGene <- x[match("Consequence", colnames(data))]
            am <- paste("c.", Ref_acid, x[match("cDNA_position", colnames(data))], Alt_acid, sep = "")
            pr <- paste("p.", gsub("/", x[match("Protein_position", colnames(data))],
                        x[match("Amino_acids", colnames(data))]), sep = "")
            tmp <- values[which(!is.na(match(values[, 2], x[match("Gene", colnames(data))]))), c(3, 1)]
            if(is.null(tmp)) return("")
            if(nrow(tmp) == 0) return("")
            tmp <- tmp[tmp[, 2] != "",]
            AAChange.refGene <- paste(apply2(tmp, 1, function(y) paste(c(y, "exon_", am, pr), collapse = ":")), collapse = ",")

            c(Chr, Start, End, toupper(Ref), toupper(Alt), Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene)
          }))
  write.table(x = data_an, file = paste(vep_file, ".annovar_format.txt", sep = ""),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  return(paste(vep_file, ".annovar_format.txt", sep = ""))
}

annotation_by_vep <- function(input_file, vep_dir, cache_dir){
  print(paste(cash_dir, "will be used to", input_file))
  if(!file.exists(input_file)) {
    print(paste(input_file, "does not exist"))
    return(NULL)
  }
  if(!file.exists(cash_dir)) {
    print("No cache detected. Using --database, it takes too long time.")
    paste(vep_dir, " --i", input_file, " --o", input_file, ".conv.txt", " --database", sep = "")
  } else {
    paste(vep_dir, " --i", input_file, " --o", input_file, ".conv.txt", " --cache --dir_cache", cash_dir, sep = "")
  }
  return(paste(input_file, ".conv.txt", sep = ""))
}

