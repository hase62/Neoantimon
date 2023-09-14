#'Calculate A Set All Neoantigen Candidates from A Given Gene Symbol and nm_id for MHC Class2 (Not yet stably available)
#'
#'@param input_nm_id (Required) An input amino acid sequence indicated as NM_ID
#'
#'@param group_ids flag to cluster the same group
#'
#'
#'
#'@param hla_file A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'
#'See by data(sample_hla_table_c1); sample_hla_table_c1;
#'
#'
#'
#'
#'
#'@param hla_types Set a list of HLA types
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name (Default=input_file).
#'
#'@param hmdir Home directory for the analysis (Default = getwd()).
#'
#'@param job_id Job-id to be attached in output files (Default = "NO_job_id").
#'
#'@param export_dir The directory will be stored results (Default = "paste("result", file_name_in_hla_table, job_id, sep=".")")
#'
#'@param peptide_length Peptide Length to be generated (Default = {8,9,10,11,12,13}).
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@param refflat_file refFlat file to be used in constructing peptide. (Default=paste(hmdir, "lib/refFlat.txt", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param refmrna_file refMrna file to be used in constructing peptide (Default=paste(hmdir, "lib/refMrna.fa", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@param netMHCIIpan_dir The file directory to netMHCpan (Default="lib/netMHCIIpan-3.2/netMHCIIpan").
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'@param ignore_short Ignore to output results of Short Peptide Less Than min(peptide_length)
#'
#'@param reading_frame The starting frame of the input sequence (Default = 1)
#'
#'@param CalculateEL Whether Calculate EL by NetMHCpan or not.
#'
#'
#'
#'@return void (Calculated Neoantigen Files will be generated as .tsv files.):
#'
#'@return HLA:  HLA type used to calculate neoantigen.
#'
#'@return Pos:  The position of a fraction of peptide used to be evaluated from the full-length peptide.
#'
#'@return Gene:  Gene symbol used to be evaluated in NetMHCpan.
#'
#'@return Evaluated_Mutant_Peptide:  The mutant peptide to be evaluated.
#'
#'@return Evaluated_Mutant_Peptide_Core:  The core peptide of the mutant peptide to be evaluated in NetMHCpan.
#'
#'@return Mut_EL: EL value for evaluated mutant peptide.
#'
#'@return Mut_Rank: Rank value for evaluated mutanat peptide.
#'
#'
#'
#'
#'
#'
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
MainEntireRegionClass2<-function(input_nm_id,
                                group_ids = seq(1:length(input_nm_id)),
                                hla_file = "here_is_a_table",
                                hla_types = NA,
                                file_name_in_hla_table = NA,
                                refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                                refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                                hmdir = getwd(),
                                job_id = "ID",
                                export_dir = paste("result", job_id, "EntireRegion2", sep="."),
                                netMHCIIpan_dir = paste(hmdir, "lib/netMHCIIpan-3.1/netMHCIIpan", sep="/"),
                                peptide_length = c(15),
                                reading_frame = 1,
                                CalculateEL = FALSE,
                                ignore_short = TRUE){

  #Get HLA-Type
  if(file.exists(hla_file) & !is.na(hla_types[1])){
    print(paste("Using:", hla_file))
  }
  if(file.exists(hla_file)){
    hla_types <- getHLAtypes(hla_file, file_name_in_hla_table)
  } else {
    hla_types <- as.character(unlist(hla_types))
  }
  if(is.na(hla_types[1])) {
    print("Please indicate hla_file and file_name_in_hla_table, or hla_types appropriately.")
    return(NULL)
  }

  #Check Required Files
  if(CheckRequiredFiles2(input_sequence = NA,
                         input_nm_id = input_nm_id,
                         hla_types = hla_types,
                         refflat_file = refflat_file,
                         refmrna_file = refmrna_file,
                         reference_nm_id = NA,
                         reference_gene_symbol = NA,
                         reading_frame)) return(NULL)

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
  names(group_ids) <- input_nm_id

  #Generate FASTA and Mutation Profile
  job_id = paste(job_id, "EntireRegion", sep = "_")
  GenerateMutatedFragments(input_sequence = NA,
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
                           reference_nm_id = NA,
                           reference_gene_symbol = NA,
                           ignore_short = ignore_short)

  #Check Output
  output_peptide_prefix <- paste(export_dir, "/", job_id, sep="")
  output_peptide_txt_file <- paste(export_dir, "/", job_id, ".peptide.txt", sep="")
  if(!file.exists(output_peptide_txt_file)){
    print("Could not Generate Mutation File for Calculating Neoantigens. Finish.")
    return(NULL)
  }

  #NetMHCIIpan
  if(is.na(netMHCIIpan_dir)){
    print("netMHCIIpan is NA.")
    return(NULL)
  }
  if(!file.exists(netMHCIIpan_dir)) {
    print(paste("Did not find", netMHCIIpan_dir))
    return(NULL)
  }
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  if(CalculateEL){
    #Execute NetMHCpan
    ExeNetMHCpanClass2(output_peptide_prefix = output_peptide_prefix,
                       "peptide",
                       hla_types,
                       netMHCIIpan_dir,
                       peptide_length,
                       export_dir,
                       input_file = "",
                       job_id)

    #Merge Results
    result <- MergeFragmentsClass2(input_dir = export_dir,
                                 file_prefix = job_id,
                                 annotation_file = output_peptide_txt_file)
  } else {
    result <- scan(output_peptide_txt_file, "character", sep = "\t")
    result <- matrix(result, byrow = TRUE, ncol = 28)
    result <- result[match(unique(result[,4]), result[, 4]), c(12, 13, 14)]
    if(!is.matrix(result)) result <- t(result)
    result <- t(sapply(sort(unique(result[,1])),
                       function(x) apply3(result[!is.na(match(result[,1], x)), c(2, 3)], 2,
                                          function(y) median(as.numeric(y)))))
    result <- apply2(result, 1, function(x) as.numeric(x[2]) / as.numeric(x[1]))
  }

  print("Successfully Finished.")
  return(result)
}
