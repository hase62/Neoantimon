#'Calculate Neoantigen Candidates on SNVs for MHC Class2
#'
#'@param input_file (Required) An input vcf file annotated by,
#'
#'e.g., ANNOVAR (http://annovar.openbioinformatics.org/en/latest/) or other softwares.
#'
#'See by data(sample_vcf); sample_vcf;
#'
#'@param hla_file (Required) A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'
#'See by data(sample_hla_table_c2); sample_hla_table_c2;
#'
#'
#'
#'
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
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name (Default=input_file).
#'
#'@param hmdir Home directory for the analysis (Default = getwd()).
#'
#'@param job_id Job-Id to be attached in output files (Default = "NO_job_id").
#'
#'@param export_dir The directory will be stored results (Default = "paste("result", file_name_in_hla_table, job_id, sep=".")")
#'
#'@param peptide_length Peptide Length to be generated (Default={15} in HLA Class2).
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
#'@param netMHCIIpan_dir The file directory to netMHCpan (Default="lib/netMHCIIpan-3.1/netMHCpan").
#'
#'@param samtools_dir The file directory to samtools_0_x_x (Default="samtools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF .
#'
#'@param bcftools_dir The file directory to netMHCpan (Default="bcftools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF .
#'samtools 0_x_x includes bcftools in the directory.
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
#'@return Mut_IC50: IC50 value for evaluated mutant peptide.
#'
#'@return Mut_Rank: Rank value for evaluated mutanat peptide.
#'
#'@return Evaluated_Wt_Peptide: The wild-type peptide to be evaluated.
#'
#'@return Wt_IC50: IC50 value for evaluated wild-type peptide.
#'
#'@return Wt_Rank: Rank value for evaluated wild-type peptide.
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
MainSNVClass2<-function(input_file,
                        hla_file = "here_is_a_table",
                        hla_types = NA,
                        file_name_in_hla_table = input_file,
                        refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                        refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                        hmdir = getwd(),
                        job_id = "ID",
                        export_dir = paste("result", file_name_in_hla_table, job_id, "SNV", sep="."),
                        rnaexp_file = NA,
                        rnabam_file = NA,
                        cnv_file=NA,
                        purity = 1,
                        netMHCIIpan_dir = paste(hmdir, "lib/netMHCIIpan-3.1/netMHCIIpan", sep="/"),
                        refdna_file = NA,
                        samtools_dir = NA,
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
                        peptide_length = c(15)){

  #Check Required Files
  if(CheckRequiredFiles(input_file = input_file,
                        hla_file = hla_file,
                        refflat_file = refflat_file,
                        refmrna_file = refmrna_file)) return(NULL)
  flg<-CheckRequiredColumns(input_file = input_file,
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
    return(NULL)
  } else {
    chr_column = flg[1]
    mutation_start_column = flg[2]
    mutation_end_column = flg[3]
    mutation_ref_column = flg[4]
    mutation_alt_column = flg[5]
    nm_id_column = flg[6]
    depth_normal_column = flg[7]
    depth_tumor_column = flg[8]
  }

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  #Generate FASTA and mutation Profile
  job_id = paste(job_id, "SNV", sep = "_")
  GenerateMutatedSeq(input_file = input_file,
                     hmdir = hmdir,
                     job_id = job_id,
                     refflat_file = refflat_file,
                     refmrna_file = refmrna_file,
                     max_peptide_length = max(peptide_length),
                     chr_column = chr_column,
                     mutation_start_column = mutation_start_column,
                     mutation_end_column = mutation_end_column,
                     mutation_ref_column = mutation_ref_column,
                     mutation_alt_column = mutation_alt_column,
                     nm_id_column = nm_id_column,
                     depth_normal_column = depth_normal_column,
                     depth_tumor_column = depth_tumor_column,
                     ambiguous_between_exon = ambiguous_between_exon,
                     ambiguous_codon = ambiguous_codon,
                     export_dir = export_dir)

  output_peptide_prefix <- paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, sep="")
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
                indel = FALSE)

  CCFP.Calc(cnv_file,
            output_peptide_txt_file,
            purity)

  #NetMHCIIpan
  if(is.na(netMHCIIpan_dir) | !file.exists(netMHCIIpan_dir)) {
    print(paste("Did not find", netMHCIIpan_dir))
    return(NULL)
  }
  print(paste("Executing netMHCIIpan to", export_dir))
  #SettingNetMHCIIpan(netMHCIIpan_dir)
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  #Get HLA-Type
  if(file.exists(hla_file)){
    hla<-t(sapply(scan(hla_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
    hit<-match(file_name_in_hla_table, hla[,1])
    if(is.na(hit)) {
      print(file_name_in_hla_table, "is not included in", hla_file)
      return (NULL)
    }
    hla_types<-hla[hit, -1]
  }
  
  #Execute NetMHCpan
  for(pep in c("peptide", "wtpeptide")){
    COUNT<-1
    output_f <- paste(output_peptide_prefix, pep, "fasta",sep=".")
    USETEMP <- FALSE
    if(nchar(output_f) > 230) {
      output_f_new <- paste("temp.Neoantimon.", runif(1) * 1000000, "txt", sep = "")
      file.copy(from = output_f, to = output_f_new)
      output_f <- output_f_new
      USETEMP <- TRUE
    }
    for(hla_type in hla_types){
      if(length(grep("DRB1", hla_type))==1) {
        system(paste(netMHCIIpan_dir,
                     " -length ", paste(peptide_length, collapse = ","),
                     " -f ", output_f,
                     " -a ", gsub("\\*","_", gsub("\\:","",hla_type)),
                     " > ", export_dir, "/", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt",
                     sep=""))
        COUNT <- COUNT + 1
      }

      if(length(grep("DPA1", hla_type))==1) {
        for(hla2 in hla_types[grep("DPB1", hla_types)]){
          system(paste(netMHCIIpan_dir,
                       " -length ", paste(peptide_length, collapse = ","),
                       " -f ", output_f,
                       " -choose -cha ", gsub("\\*|\\:","", hla_type),
                       " -choose -chb ", gsub("\\*|\\:","", hla2),
                       " > ", export_dir, "/", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt",
                       sep=""))
          COUNT <- COUNT + 1
        }
      }

      if(length(grep("DQA1", hla_type))==1) {
        for(hla2 in hla_types[grep("DQB1", hla_types)]){
          system(paste(netMHCIIpan_dir,
                       " -length ", paste(peptide_length, collapse = ","),
                       " -f ", output_f,
                       " -choose -cha ", gsub("\\*|\\:","", hla_type),
                       " -choose -chb ", gsub("\\*|\\:","", hla2),
                       " > ", export_dir, "/", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt",
                       sep=""))
          COUNT <- COUNT + 1
        }
      }
    }
    if(USETEMP) file.remove(output_f)
  }
  print("Merging Results...")
  result <- MergeSNVClass2(input_dir = export_dir,
                           file_prefix = job_id,
                           annotation_file = output_peptide_txt_file)

  print("Successfully Finished.")
  return(result)
}
