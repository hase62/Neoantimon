#'Calculate Neoantigen Candidates from A Given Sequence for MHC Class2
#'
#'@param input_sequence (Required) An input amino acid sequence
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
#'@param job_id Job-Id to be attached in output files (Default = "NO_job_id").
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
#'@param reference_nm_id Corresponding original sequences that the input sequence is generated.
#'If franctions of peptides generated from the input are included in the indicated protein, such peptides are removed.
#'It can be indicated when gene_symbol is not NA.
#'
#'@param reference_gene_symbol Corresponding original sequences that the input sequence is generated.
#'If franctions of peptides generated from the input are included in the indicated protein, such peptides are removed.
#'It can be indicated when nm_id is not NA.
#'

#'@param IgnoreShortPeptides Ignore to output results of Short Peptide Less Than min(peptide_length)
#'
#'
#'
#'
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
#'@return Mut_IC50: IC50 value for evaluated mutant peptide.
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
MainSeqFragmentClass2<-function(input_sequence = NA,
                                group_ids = seq(1:length(reference_nm_id)),
                                hla_file = "here_is_a_table",
                                hla_types = NA,
                                file_name_in_hla_table = NA,
                                refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                                refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                                hmdir = getwd(),
                                job_id = "ID",
                                export_dir = paste("result", job_id, "SeqFragment2", sep="."),
                                netMHCIIpan_dir = paste(hmdir, "lib/netMHCIIpan-3.1/netMHCIIpan", sep="/"),
                                peptide_length = c(15),
                                reference_nm_id = NA,
                                reference_gene_symbol = NA,
                                IgnoreShortPeptides = TRUE){

  #Get HLA-Type
  if(file.exists(hla_file) & !is.na(hla_types[1])){
    print(paste("Using:", hla_file))
  }
  if(file.exists(hla_file)){
    hla_types <- getHLAtypes(hla_file, file_name_in_hla_table)
  }
  if(is.na(hla_types[1])) {
    print("Please indicate hla_file and file_name_in_hla_table, or hla_types appropriately.")
    return(NULL)
  }

  #Check Required Files
  if(CheckRequiredFiles2(input_sequence = input_sequence,
                         input_nm_id = NA,
                         hla_types = hla_types,
                         refflat_file = refflat_file,
                         refmrna_file = refmrna_file,
                         reference_nm_id,
                         reference_gene_symbol,
                         1)) return(NULL)

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
  names(group_ids) <- input_sequence

  #Generate FASTA and Mutation Profile
  job_id = paste(job_id, "SeqFragment", sep = "_")
  GenerateMutatedFragments(input_sequence = input_sequence,
                           input_nm_id = NA,
                           group_ids = group_ids,
                           hmdir = hmdir,
                           job_id = job_id,
                           refflat_file = refflat_file,
                           refmrna_file = refmrna_file,
                           max_peptide_length = max(peptide_length),
                           min_peptide_length = min(peptide_length),
                           reading_frame = 1,
                           export_dir = export_dir,
                           reference_nm_id = reference_nm_id,
                           reference_gene_symbol = reference_gene_symbol,
                           IgnoreShortPeptides = IgnoreShortPeptides)

  #Check Output
  output_peptide_prefix <- paste(export_dir, "/", job_id, sep="")
  output_peptide_txt_file <- paste(export_dir, "/", job_id, ".peptide.txt", sep="")
  if(!file.exists(output_peptide_txt_file)){
    print("Could not Generate Mutation File for Calculating Neoantigens. Finish.")
    return(NULL)
  }

  #NetMHCIIpan
  if(is.na(netMHCIIpan_dir) | !file.exists(netMHCIIpan_dir)) {
    print(paste("Did not find", netMHCIIpan_dir))
    return(NULL)
  }

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

  print("Successfully Finished.")
  return(result)
}
