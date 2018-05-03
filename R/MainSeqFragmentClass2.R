#'Calculate Neoantigen Candidates on SV fusions for MHC Class2
#'
#'@param input_file (Required) An input vcf file (BND format) annotated by,
#'
#'e.g., ANNOVAR (http://annovar.openbioinformatics.org/en/latest/) or other softwares.
#'
#'See by data(sample_sv_bnd); sample_sv_bnd;
#'
#'@param hla_file (Required) A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'
#'See by data(sample_hla_table_c1); sample_hla_table_c1;
#'
#'@param refdna_file (Required) refdna_file information to be used to create SVs Region (Default=NA).
#'
#'See "https://github.com/hase62/Neoantimon"
#'
#'@param nm_id_column (Required if gene_symbol_column = NA) The column number describing NM IDs in input_file such as
#'
#'"SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs" (Default=NA).
#'
#'@param gene_symbol_column (Required if nm_id_column = NA) The column number describing gene symbol in input_file (Default=NA).
#'
#'@param mate_id_column (Required) The column indicating mateIDs or svIDs such as "SVMERGE1_1" (Default=NA).
#'
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name (Default=input_file).
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
#'@param mutation_alt_bnd_column The column number describing mutation Alt (BND format) in input_file (Default=NA, but will automatically search "Alt" in header).
#'
#'@param depth_normal_column The column number describing the read count from normal cells (Default = NA).
#'
#'@param depth_tumor_column The column number describing the read count from tumor cells (Default = NA).
#'
#'@param ambiguous_between_exon The maximum number to permit the differences between Exon-Lengths from refFlat and refMrna (Default=0).
#'
#'@param ambiguous_codon The maximum number to permit the differences between inputfile- and refMrna-oriented translation start/end position (Default=0).
#'
#'@param refflat_file refFlat file to be used in constructing peptide. (Default=paste(hmdir, "lib/refFlat.txt",sep="").
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
#'
#'
#'
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
#'@return void (Calculated Neoantigen Files will be generated as .tsv files.)
#'
#'@export
MainSeqFragmentClass2<-function(input_sequence,
                                hla_file,
                                file_name_in_hla_table,
                                refflat_file = paste(hmdir, "lib/refFlat.txt", sep="/"),
                                refmrna_file = paste(hmdir, "lib/refMrna.fa", sep="/"),
                                hmdir = getwd(),
                                job_id = "NO_job_id",
                                export_dir = paste("result", file_name_in_hla_table, job_id, "SeqFragment", sep="."),
                                netMHCIIpan_dir = paste(hmdir, "lib/netMHCIIpan-3.1/netMHCIIpan", sep="/"),
                                peptide_length = c(15),
                                nm_id,
                                gene_symbol,
                                reading_frame = 1){

  #Check Required Files
  if(CheckRequiredFiles2(input_sequence = input_sequence,
                         hla_file = hla_file,
                         refflat_file = refflat_file,
                         refmrna_file = refmrna_file,
                         nm_id,
                         gene_symbol,
                         reading_frame)) return(NULL)

  #Make Directory
  if(!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

  #Generate FASTA and Mutation Profile
  job_id = paste(job_id, "SeqFragment", sep = "_")
  GenerateMutatedFragments(input_sequence = input_sequence,
                           hmdir = hmdir,
                           job_id = job_id,
                           refflat_file = refflat_file,
                           refmrna_file = refmrna_file,
                           max_peptide_length = max(peptide_length),
                           min_peptide_length = min(peptide_length),
                           reading_frame = reading_frame,
                           export_dir = export_dir)

  #Check Output
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
  print(paste("Executing netMHCIIpan to", export_dir))
  SettingNetMHCIIpan(netMHCIIpan_dir)

  #Get HLA-Type
  hla <- t(sapply(scan(hla_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  hit <- match(file_name_in_hla_table, hla[,1])
  if(is.na(hit)) {
    print(file_name_in_hla_table, "is not included in", hla_file)
    return (NULL)
  }

  #Execute NetMHCpan
  hla_types <- hla[hit, -1]
  for(pep in c("peptide")){
    COUNT<-1
    for(hla_type in hla_types){
      if(length(grep("DRB1", hla_type))==1) {
        system(paste(netMHCIIpan_dir,
                     " -length ", paste(peptide_length, collapse = ","),
                     " -f ", paste(export_dir, "/", job_id, ".", pep, ".", "fasta", sep=""),
                     " -a ", gsub("\\*","_", gsub("\\:","",hla_type)),
                     " > ", export_dir, "/", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt",
                     sep=""))
        COUNT <- COUNT + 1
      }

      if(length(grep("DPA1", hla_type))==1) {
        for(hla2 in hla_types[grep("DPB1", hla_types)]){
          system(paste(netMHCIIpan_dir,
                       " -length ", paste(peptide_length, collapse = ","),
                       " -f ", paste(export_dir, "/", job_id, ".", pep, ".", "fasta", sep=""),
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
                       " -f ", paste(export_dir, "/", job_id, ".", pep, ".", "fasta", sep=""),
                       " -choose -cha ", gsub("\\*|\\:","", hla_type),
                       " -choose -chb ", gsub("\\*|\\:","", hla2),
                       " > ", export_dir, "/", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt",
                       sep=""))
          COUNT <- COUNT + 1
        }
      }
    }
  }
  result <- MainMergeINDELSVClass2(input_dir = export_dir,
                                   file_prefix = job_id,
                                   annotation_file = output_peptide_txt_file)

  print("Successfully Finished.")
  return(result)
}