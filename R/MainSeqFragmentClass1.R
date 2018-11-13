#'Calculate Neoantigen Candidates from A Given Sequence for MHC Class1
#'
#'@param input_sequence (Required) An input amino acid sequence
#'
#'@param input_nm_id (Required) An input amino acid sequence indicated as NM_ID
#'
#'@param group_ids flag to cluster the same group
#'
#'@param hla_file (Required) A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'
#'See by data(sample_hla_table_c1); sample_hla_table_c1;
#'
#'@param hla_types Set a list of HLA types
#'
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name (Default=input_file).
#'
#'@param hmdir Home directory for the analysis (Default = getwd()).
#'
#'@param job_id Job-Id to be attached in output files (Default = "NO_job_id").
#'
#'@param export_dir The directory will be stored results (Default = "paste("result", file_name_in_hla_table, job_id, sep=".")")
#'
#'@param peptide_length Peptide Length to be generated (Default = {8,9,10,11,12,13}).
#'
#'@param refflat_file refFlat file to be used in constructing peptide. (Default=paste(hmdir, "lib/refFlat.txt", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param refmrna_file refMrna file to be used in constructing peptide (Default=paste(hmdir, "lib/refMrna.fa", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param netMHCpan_dir The file directory to netMHCpan (Default="lib/netMHCpan-4.0/netMHCpan").
#'
#'@param reference_nm_id Corresponding original sequences that the input sequence is generated.
#'If franctions of peptides generated from the input are included in the indicated protein, such peptides are removed.
#'It can be indicated when gene_symbol is not NA.
#'
#'@param reference_gene_symbol Corresponding original sequences that the input sequence is generated.
#'If franctions of peptides generated from the input are included in the indicated protein, such peptides are removed.
#'It can be indicated when nm_id is not NA.
#'
#'@param reading_frame The starting frame of the input sequence (Default = 1)
#'
#'@return void (Calculated Neoantigen Files will be generated as .tsv files.):
#'
#'@return HLA:  HLA type used to calculate neoantigen.
#'
#'@return Pos:  The position of a fraction of peptide used to be evaluated from the full-length peptide.
#'
#'@return Gene:  Gene symbol used to be evaluated in NetMHCpan.
#'
#'@return Evaluated_Mutant_Peptide_Core:  The core peptide of the mutant peptide to be evaluated in NetMHCpan.
#'
#'@return Evaluated_Mutant_Peptide:  The mutant peptide to be evaluated.
#'
#'@return Mut_IC50: IC50 value for evaluated mutant peptide.
#'
#'@return Mut_Rank: Rank value for evaluated mutanat peptide.
#'
#'@return Chr: Chromosome Number of the mutation.
#'
#'@return NM_ID: NM_ID used to construct peptides from the mutation.
#'
#'@return Change: The annotation to be described in .vcf file.
#'
#'@return Ref: reference type nucleic acid base.
#'
#'@return Alt: alternative type nucleic acid base.
#'
#'@return Prob: A probability of reference nucleic acid base described in .vcf file.
#'
#'@return Mutation_Prob: A probability of alternative nucleic acid base described in .vcf file.
#'
#'@return Exon_Start: The exon start position of the corrsponding NM_ID.
#'
#'@return Exon_End: The exon end position of the corrsponding NM_ID.
#'
#'@return Mutation_Position: The mutation position of the corrsponding NM_ID.
#'
#'@return Total_Depth: The sum depth of the reference and alternative nucleic acid base.
#'
#'@return Tumor_Depth: The depth of the alternative nucleic acid base.
#'
#'@return Wt_Peptide: The full-length of the wild-type peptide.
#'
#'@return Mutant_Peptide: The full-length of the mutant peptide.
#'
#'@return Total_RNA: The expression amount of the corresponding RNA.
#'
#'@return Tumor_RNA_Ratio: The variant allele frequency of the corresponding RNA.
#'
#'@return Tumor_RNA: The modified amount of the corresponding RNA level based on RNA Reads.
#'
#'@return Tumor_RNA_based_on_DNA: The modified amount of the corresponding RNA level based on DNA Reads.
#'
#'@return MutRatio: The mean value of the cancer cell fraction probability.
#'
#'@return MutRatio_Min: The 1\% percentile of the cancer cell fraction probability.
#'
#'@return MutRatio_Max: The 99\% percentile of the cancer cell fraction probability.
#'
#'@export
MainSeqFragmentClass1<-function(input_sequence = NA,
                                input_nm_id = NA,
                                group_ids = NA,
                                hla_file = "here_is_a_table",
                                hla_types = NA,
                                file_name_in_hla_table = NA,
                                refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                                refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                                hmdir = getwd(),
                                job_id = "ID",
                                export_dir = paste("result", file_name_in_hla_table, job_id, "SeqFragment", sep="."),
                                netMHCpan_dir = paste(hmdir, "lib/netMHCpan-4.0/netMHCpan", sep="/"),
                                peptide_length = c(8, 9, 10, 11, 12, 13),
                                reference_nm_id = NA,
                                reference_gene_symbol = NA,
                                reading_frame = 1){

  #Check Required Files
  if(CheckRequiredFiles2(input_sequence = input_sequence,
                         input_nm_id = input_nm_id,
                         hla_file = hla_file,
                         hla_types = hla_types,
                         refflat_file = refflat_file,
                         refmrna_file = refmrna_file,
                         reference_nm_id,
                         reference_gene_symbol,
                         reading_frame)) return(NULL)

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  #Attach Group IDs
  if(is.null(group_ids[1])) group_ids <- seq(from = 1,
                                             to = ifelse(is.na(input_sequence[1]), 0, length(input_sequence)) +
                                               ifelse(is.na(input_nm_id[1]), 0, length(input_nm_id)))
  tmp <- c(input_sequence, input_nm_id)
  names(group_ids) <- tmp[!is.na(tmp)]

  #Generate FASTA and Mutation Profile
  job_id = paste(job_id, "SeqFragment", sep = "_")
  GenerateMutatedFragments(input_sequence = input_sequence,
                           input_nm_id = input_nm_id,
                           group_ids = group_ids,
                           hmdir = hmdir,
                           job_id = job_id,
                           refflat_file = refflat_file,
                           refmrna_file = refmrna_file,
                           max_peptide_length = max(peptide_length),
                           min_peptide_length = min(peptide_length),
                           reading_frame = reading_frame,
                           export_dir = export_dir,
                           reference_nm_id = reference_nm_id,
                           reference_gene_symbol = reference_gene_symbol)

  #Check Output
  output_peptide_txt_file <- paste(export_dir, "/", job_id, ".peptide.txt", sep="")
  if(!file.exists(output_peptide_txt_file)){
    print("Could not Generate Mutation File for Calculating Neoantigens. Finish.")
    return(NULL)
  }

  #NetMHCpan
  if(is.na(netMHCpan_dir) | !file.exists(netMHCpan_dir)) {
    print(paste("Did not find", netMHCpan_dir))
    return(NULL)
  }
  print(paste("Executing netMHCpan to", export_dir))
  ##SettingNetMHCpan(netMHCpan_dir)
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  #Get HLA-Type
  if(file.exists(hla_file)){
    hla <- t(sapply(scan(hla_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
    hit <- match(file_name_in_hla_table, hla[,1])
    if(is.na(hit)) {
      print(file_name_in_hla_table, "is not included in", hla_file)
      return (NULL)
    }
    return (NULL)
    hla_types <- hla[hit, -1]
  }

  #Execute NetMHCpan
  for(pep in c("peptide")){
    COUNT<-1
    output_f <- paste(export_dir, "/", job_id, ".", pep, ".", "fasta", sep="")
    USETEMP <- FALSE
    if(nchar(output_f) > 230) {
      output_f_new <- paste("temp.Neoantimon.", runif(1) * 1000000, "txt", sep = "")
      file.copy(from = output_f, to = output_f_new)
      output_f <- output_f_new
      USETEMP <- TRUE
    }
    for(hla_type in hla_types){
      paste("Calculating", pep, hla_type)
      system(paste(netMHCpan_dir,
                   " -BA ",
                   " -l ", paste(peptide_length, collapse = ","),
                   " -f ", output_f,
                   " -a HLA-", gsub("\\*", "", hla_type),
                   " > ", export_dir, "/", job_id, ".HLACLASS1.", COUNT, ".", pep, ".txt", sep=""))
      COUNT <- COUNT + 1
    }
    if(USETEMP) file.remove(output_f)
  }
  print("Merging Results...")
  result <- MergeINDELSVClass1(input_dir = export_dir,
                               file_prefix = job_id,
                               annotation_file = output_peptide_txt_file)

  print("Successfully Finished.")
  return(result)
}
