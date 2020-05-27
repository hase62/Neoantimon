convert_to_annovar_format_from_vep <- function(vep_file) {
  #Read vep data
  if(is.list(vep_file) | is.matrix(vep_file)){
    data <- as.matrix(vep_file)
    vep_file <- paste("data", round(runif(1) * 10000), sep = ".")
  } else {
    data <- read_data(vep_file)
  }
  data <- data[grep("missense_variant|insertion|deletion|frameshift", data[, match("Consequence", colnames(data))]), ]
  data[, match("Feature", colnames(data))] <- sapply(data[, match("Feature", colnames(data))],
                                                     function(x) ifelse(length(grep("NM", x) == 1), strsplit(x, "\\.")[[1]][1], x))

  #Execute ensembl
  print("Executing Transformation")
  ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
  values_rna_1 <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol", "ensembl_transcript_id_version"),
                        filters = "ensembl_transcript_id_version", values = data[, match("Feature", colnames(data))], mart= ensembl)
  values_rna_2 <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol", "refseq_mrna"),
                        filters = "refseq_mrna", values = data[, match("Feature", colnames(data))], mart= ensembl)
  colnames(values_rna_2)[4] <- "ensembl_transcript_id_version"
  values_rna <- rbind(values_rna_1, values_rna_2)
  values_rna <- values_rna[values_rna[, 1] != "",]
  if(nrow(values_rna) == 0){
    use_trans_id_version = FALSE
    values_dna <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"),
                        filters = "ensembl_gene_id", values = data[, match("Gene", colnames(data))], mart= ensembl)
    values_dna <- values_dna[values_dna[, 1] != "",]
  } else {
    use_trans_id_version = TRUE
  }

  #Make annovar format data
  index <- c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
  data_an <- matrix(ncol = length(index), nrow = nrow(data), ".")
  colnames(data_an) <- index
  data_an[, match(c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"), colnames(data_an))] <-
          t(apply2(data, 1, function(x) {
            x <- as.character(x)
            tmp <- strsplit(x[match("Location", colnames(data))], ":|-")[[1]]
            Chr <- tmp[1]
            Start <- tmp[2]
            End <- tmp[2]
            Ref_acid <- gsub("a|t|g|c", "", strsplit(x[match("Codons", colnames(data))], "/")[[1]][1])
            if(Ref_acid == "") {
              Ref <- "-"
            } else if(length(grep("STRAND=1", x[match("Extra", colnames(data))])) == 1){
              Ref <- Ref_acid
            } else {
              Ref <- paste(trans_to[match(tolower(strsplit(Ref_acid, "")[[1]]), trans_from)], collapse = "")
            }
            Alt_acid <- gsub("a|t|g|c", "", strsplit(x[match("Codons", colnames(data))], "/")[[1]][2])
            if(Alt_acid == "") {
              Alt <- "-"
            } else if(length(grep("STRAND=1", x[match("Extra", colnames(data))])) == 1){
              Alt <- Alt_acid
            } else {
              Alt <- paste(trans_to[match(tolower(strsplit(Alt_acid, "")[[1]]), trans_from)], collapse = "")
            }
            Func.refGene <- "exonic"
            if(use_trans_id_version){
              Gene.refGene <- paste(unique(values_rna[which(!is.na(match(values_rna[, 4], x[match("Feature", colnames(data))]))), 3]), collapse = ";")
              if(Gene.refGene == "") return(rep("", length(index) - 1))
            } else {
              Gene.refGene <- paste(unique(values_dna[which(!is.na(match(values_dna[, 2], x[match("Gene", colnames(data))]))), 3]), collapse = ";")
            }
            ExonicFunc.refGene <- x[match("Consequence", colnames(data))]
            if(Ref_acid == "-"){
              am <- paste("c.", gsub("-", "_", x[match("cDNA_position", colnames(data))]), "ins", sep = "")
              pr <- paste("p.", gsub("/",
                                     x[match("Protein_position", colnames(data))],
                                     paste(strsplit(x[match("Amino_acids", colnames(data))], "/")[[1]][1], "/ins", sep = "")), sep = "")
            } else if (Alt_acid == "-"){
              am <- paste("c.", gsub("-", "_", x[match("cDNA_position", colnames(data))]), "del", sep = "")
              pr <- paste("p.", gsub("/",
                                     x[match("Protein_position", colnames(data))],
                                     paste(strsplit(x[match("Amino_acids", colnames(data))], "/")[[1]][1], "/del", sep = "")), sep = "")

            } else {
              am <- paste("c.", Ref_acid, x[match("cDNA_position", colnames(data))], Alt_acid, sep = "")
              pr <- paste("p.", gsub("/", x[match("Protein_position", colnames(data))], x[match("Amino_acids", colnames(data))]), sep = "")
            }
            if(use_trans_id_version){
              tmp <- values_rna[which(!is.na(match(values_rna[, 4], x[match("Feature", colnames(data))]))), c(3, 1)]
            } else {
              tmp <- values_dna[which(!is.na(match(values_dna[, 2], x[match("Gene", colnames(data))]))), c(3, 1)]
            }
            if(is.null(tmp)) return(rep("", length(index) - 1))
            if(nrow(tmp) == 0) return(rep("", length(index) - 1))
            tmp <- tmp[tmp[, 2] != "",]
            AAChange.refGene <- paste(apply2(tmp, 1, function(y) paste(c(y, "exonX", am, pr), collapse = ":")), collapse = ",")
            c(Chr, Start, End, toupper(Ref), toupper(Alt), Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene)
          }))
  if(!is.matrix(data_an)) data_an <- t(data_an)
  write.table(x = data_an[nchar(data_an[, 1]) > 0, ], file = paste(vep_file, ".annovar_format.txt", sep = ""),
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

