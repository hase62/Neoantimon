#'Calculate Neoantigen Candidates on INDELs for MHC Class2
#'
#'@param input_annovar_format_file An input vcf file annotated by ANNOVAR (http://annovar.openbioinformatics.org/en/latest/).
#'You can directly indicate a matrix, which is the same as annovar format vcf file, as input.
#'
#'See by data(sample_vcf.annovar); sample_vcf.annovar.txt;
#'
#'@param input_vep_format_file An input file annotated by Ensembl Variant Effect Predictor (VEP).
#'You can directly indicate a matrix, which is the same as annovar format VEP file, as input.
#'
#'See by data(sample_vcf.vep); sample_vcf.vep.txt;
#
#'@param input_vcf_format_file_and_vep An input vcf file and path to Ensembl Variant Effect Predictor (VEP).
#'Before using this option, please install vep according to the official cite ("https://asia.ensembl.org/info/docs/tools/vep/index.html").
#'
#'@param hla_file A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'
#'See by data(sample_hla_table_c1); sample_hla_table_c1;
#'
#'@param hla_types Set a list of HLA types
#'
#'@param nm_id_column The column number describing NM IDs in input_file such as
#'
#'"SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs" (Default=NA).
#'
#'
#'
#'
#'
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name.
#'
#'@param hmdir Home directory for the analysis (Default = getwd()).
#'
#'@param job_id Job-Id to be attached in output files (Default = "NO_job_id").
#'
#'@param export_dir The directory will be stored results (Default = "paste("result", file_name_in_hla_table, job_id, sep=".")")
#'
#'@param peptide_length Peptide Length to be generated (Default = {15} in HLA Class2).
#'
#'@param chr_column The column number describing chromosome number in input_file (Default=NA, but will automatically search "Chr" in header).
#'
#'@param mutation_start_column The column number describing mutation start Position in input_file (Default=NA, but will automatically search "Start" in header) .
#'
#'@param mutation_end_column The column number describing mutation end Position in input_file (Default=NA, but will automatically search "End" in header).
#'
#'@param mutation_ref_column The column number describing mutation Ref in input_file (Default=NA, but will automatically search "Ref" in header).
#'
#'@param mutation_alt_column The column number describing mutation Alt in input_file (Default=NA, but will automatically search "Alt" in header).
#'
#'@param depth_normal_column The column number describing the read count from normal cells (Default = NA).
#'
#'@param depth_tumor_column The column number describing the read count from tumor cells (Default = NA).
#'
#'@param ambiguous_between_exon The maximum number to permit the differences between Exon-Lengths from refFlat and refMrna (Default=0).
#'
#'@param ambiguous_codon The maximum number to permit the differences between inputfile- and refMrna-oriented translation start/end position (Default=0).
#'
#'@param refflat_file refFlat file to be used in constructing peptide. (Default=paste(hmdir, "lib/refFlat.txt", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param refmrna_file refMrna file to be used in constructing peptide (Default=paste(hmdir, "lib/refMrna.fa", sep="").
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param rnaexp_file A file including RNA expressions (Default=NA).
#'The 1st, 2nd and 3rd columns are "GeneSymbol Chr:Exonstart-Exonend (locus) ExpressionAmount", respectively.
#'The 1st row should be any header.
#'
#'See by data(sample_rna_exp); sample_rna_exp;
#'
#'@param rnabam_file RNA bam file to calculate variant allele frequency of RNA at each mutation (Default=NA).
#'
#'@param refdna_file refdna_file information to be used to calculate RNA VAF (Default=NA).
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param cnv_file A file including copy number variation to calculate cancer cell fraction probability (CCFP) (Default=NA).
#'The format is according to ASCAT output files.
#'The columns are "SNPName Chromosome Position LogR segmentedLogR BAF segmentedBAF CopyNumber MinorAllele RawCopyNumber"
#'The 1st row should be the above header.
#'
#'See data(sample_copynum); sample_copynum;
#'
#'@param purity Tumor purity or tumor contents ratio required to calculate CCFP (Default=1).
#'
#'@param netMHCIIpan_dir The file directory to netMHCpan (Default="lib/netMHCIIpan-3.2/netMHCpan").
#'
#'
#'
#'@param samtools_dir The file directory to samtools_0_x_x (Default="samtools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF.
#'
#'@param bcftools_dir The file directory to netMHCpan (Default="bcftools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF.
#'samtools 0_x_x includes bcftools in the directory.
#'
#'@param ignore_short Ignore to output results of short peptide less than min (peptide_length)
#'
#'@param SNPs Apply indivisual SNPs on peptides by indicate a vcf file.
#'
#'@param multiple_variants Reflect multiple variants on a peptide, e.g., SNVs on frameshift region.
#'
#'@param base_0 Apply 0-base format for mutation positions. 
#'
#'@return void (Calculated Neoantigen Files will be generated as .tsv files.):
#'
#'@return HLA: HLA type used to calculate neoantigen.
#'
#'@return Pos: The position of a fraction of peptide used to be evaluated from the full-length peptide.
#'
#'@return Gene Gene symbol used to be evaluated in NetMHCpan.
#'
#'@return Evaluated_Mutant_Peptide: The mutant peptide to be evaluated.
#'
#'@return Evaluated_Mutant_Peptide_Core: The core peptide of the mutant peptide to be evaluated in NetMHCpan.
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
#'@return P_I: Priority score using the EL.
#'Please use CalculatePriorityScores <- function(result, useRNAvaf = FALSE)
#'
#'@return P_R: Priority score using the percentage of rank affinity.
#'Please use CalculatePriorityScores <- function(result, useRNAvaf = FALSE)

#'@return P: Priority score implemented in MuPeXI (Bjerregaard et al. 2017).
#'Please use CalculatePriorityScores <- function(result, useRNAvaf = FALSE)

#'@export
MainINDELClass2<-function(input_annovar_format_file = NA,
                          input_vep_format_file = NA,
                          input_vcf_format_file_and_vep = NA,
                          hla_file = "here_is_a_table",
                          hla_types = NA,
                          file_name_in_hla_table = "sample",
                          refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                          refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                          hmdir = getwd(),
                          job_id = "ID",
                          export_dir = paste("result", job_id, "INDEL2", sep="."),
                          rnaexp_file = NA,
                          rnabam_file = NA,
                          cnv_file=NA,
                          purity = 1,
                          netMHCIIpan_dir = paste(hmdir, "lib/netMHCIIpan-3.1/netMHCIIpan", sep="/"),

                          refdna_file = NA,
                          samtools_dir = "samtools",
                          bcftools_dir = NA,
                          chr_column = NA,
                          mutation_start_column = NA,
                          mutation_end_column = NA,
                          mutation_ref_column = NA,
                          mutation_alt_column = NA,
                          nm_id_column = NA,
                          depth_normal_column = NA,
                          depth_tumor_column = NA,
                          ambiguous_between_exon = 0,
                          ambiguous_codon = 0,
                          peptide_length = c(15),
                          ignore_short = TRUE,
                          SNPs = NA,
                          multiple_variants = FALSE, 
                          base_0 = FALSE){

  #Obtain Data
  if(Read_files(input_annovar_format_file, input_vep_format_file, input_vcf_format_file_and_vep)) return(NULL)
  if(!is.na(input_vcf_format_file_and_vep)) input_vep_format_file <- annotation_by_vep(input_vcf_format_file_and_vep[1],
                                                                                       input_vcf_format_file_and_vep[2],
                                                                                       input_vcf_format_file_and_vep[3])
  if(is.null(input_vep_format_file)) {
    print("Failed to annotate by VEP.")
    return(NULL)
  }
  if(is.list(input_vep_format_file) | is.matrix(input_vep_format_file)) {
    input_annovar_format_file <- convert_to_annovar_format_from_vep(input_vep_format_file)
  } else if(!is.na(input_vep_format_file)) {
    input_annovar_format_file <- convert_to_annovar_format_from_vep(input_vep_format_file)
  }

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
  if(CheckRequiredFiles(input_file = input_annovar_format_file,
                        hla_types = hla_types,
                        refflat_file = refflat_file,
                        refmrna_file = refmrna_file)) return(NULL)

  flg <- CheckRequiredColumns(input_file = input_annovar_format_file,
                              chr_column = chr_column,
                              mutation_start_column = mutation_start_column,
                              mutation_end_column = mutation_end_column,
                              mutation_ref_column = mutation_ref_column,
                              mutation_alt_column = mutation_alt_column,
                              nm_id_column = nm_id_column,
                              depth_normal_column = depth_normal_column,
                              depth_tumor_column = depth_tumor_column)

  #Check and Set Required Columns
  if(length(flg)<=1) {
    print("There is no available column names.")
    return(NULL)
  }

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)













  #Generate FASTA and mutation Profile
  job_id <- paste(job_id, "INDEL", sep = "_")

  GenerateIndelSeq(input_file = input_annovar_format_file,
                   hmdir = hmdir,
                   job_id = job_id,
                   refflat_file = refflat_file,
                   refmrna_file = refmrna_file,
                   max_peptide_length = max(peptide_length),
                   chr_column = flg[1],
                   mutation_start_column = flg[2],
                   mutation_end_column = flg[3],
                   mutation_ref_column = flg[4],
                   mutation_alt_column = flg[5],
                   nm_id_column = flg[6],
                   depth_normal_column = flg[7],
                   depth_tumor_column = flg[8],
                   ambiguous_between_exon = ambiguous_between_exon,
                   ambiguous_codon = ambiguous_codon,
                   export_dir = export_dir,
                   ignore_short = ignore_short,
                   SNPs = SNPs,
                   multiple_variants = multiple_variants, 
                   base_0 = base_0)

  #Check Generated File exists
  if(is.list(input_annovar_format_file) | is.matrix(input_annovar_format_file)) input_annovar_format_file <- "data"
  output_peptide_prefix <- paste(export_dir, "/", rev(strsplit(input_annovar_format_file, "/")[[1]])[1], ".", job_id, sep="")
  output_peptide_txt_file <- paste(output_peptide_prefix, ".peptide.txt", sep="")
  if(!file.exists(output_peptide_txt_file)){
    print("Could not Generate Mutation File for Calculating Neoantigens. Finish.")
    return(NULL)
  }

  RNAExpression(rnaexp_file,
                output_peptide_txt_file,
                width = 2,
                samtools_dir,
                refdna_file,
                rnabam_file,
                bcftools_dir,
                indel = TRUE)

  CCFP.Calc(cnv_file,
            output_peptide_txt_file,
            purity)

  #NetMHCIIpan
  if(is.na(netMHCIIpan_dir)){
    print("netMHCIIpan is NA.")
    return(NULL)
  }
  if(!file.exists(netMHCIIpan_dir)) {
    print(paste("Did not find", netMHCIIpan_dir))
    return(NULL)
  }

  #Execute NetMHCpan
  ExeNetMHCpanClass2(output_peptide_prefix,
                     "peptide",
                     hla_types,
                     netMHCIIpan_dir,
                     peptide_length,
                     export_dir,
                     input_annovar_format_file,
                     job_id)

  #Merge Results
  result <- MergeINDELSVClass2(input_dir = export_dir,
                               file_prefix = paste(rev(strsplit(input_annovar_format_file, "/")[[1]])[1], job_id, sep = "."),
                               annotation_file = output_peptide_txt_file)
















  print("Successfully Finished.")
  return(result)
}











