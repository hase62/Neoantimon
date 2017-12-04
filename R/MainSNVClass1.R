#'Calculate Neoantigen Candidates on SNVs for MHC Class1
#'
#'@param input_file (Required) An input vcf file annotated by, e.g., ANNOVAR (http://annovar.openbioinformatics.org/en/latest/) or other softwares.
#'See by data(sample_vcf); sample_vcf;
#'
#'@param hla_file (Required) A tab separated file indicating HLA types.
#'The 1st column is input_file name, and the following columns indicate HLA types.
#'See by data(sample_hla_table_c1); sample_hla_table_c1;
#'
#'@param file_name_in_hla_table If the name (1st column) in HLA table is not the same as input_file, indicate the corresponding name (Default=input_file).
#'
#'@param hmdir Home directory for the analysis (Default=getwd()).
#'
#'@param job_id Job-Id to be attached in output files (Default="NO_job_id").
#'
#'@param peptide_length Peptide Length to be generated (Default={8,9,10,11,12,13}).
#'
#'@param chr_column The column number describing chromosome number in input_file (Default=1).
#'
#'@param mutation_start_column The column number describing mutation start Position in input_file (Default=2).
#'
#'@param mutation_end_column The column number describing mutation end Position in input_file (Default=3).
#'
#'@param mutation_ref_column The column number describing mutation Ref in input_file (Default=4).
#'
#'@param mutation_alt_column The column number describing mutation Alt in input_file (Default=5).
#'
#'@param nm_id_column The column number describing NM IDs in input_file (Default=10).
#'
#'@param depth_normal_column The column number describing the read count from normal cells (Default = NA). 
#'
#'@param depth_tumor_column The column number describing the read count from tumor cells (Default = NA). 
#'
#'@param ambiguous_between_exon The maximum number to permit the differences between Exon-Lengths from refFlat and refMrna (Default=0).
#'
#'@param ambiguous_codon The maximum number to permit the differences between inputfile- and refMrna-oriented translation start/end position (Default=0).
#'
#'@param refflat_file refFlat file to be used in constructing peptide. (Default=paste(hmdir,"lib/refFlat.txt",sep="").
#'
#'@param refmrna_file refMrna file to be used in constructing peptide (Default=paste(hmdir,"lib/refMrna.merge.fa",sep="").
#'The first column is NM_ID, the second column is transcript variant number, and the third column is the mRNA sequence.
#'This file is automaticalluy generated through the command in README.
#'
#'@param rnaexp_file A file including RNA expressions (Default=NA).
#'The 1st, 2nd and 3rd columns are "GeneSymbol Chr:Exonstart-Exonend(locus) Expression Amount", respectively.
#'The 1st row should be any header.
#'See by data(sample_rna_exp); sample_rna_exp;
#'
#'@param rnabam_file RNA bam file to calculate variant allele frequency of RNA at each mutation (Default=NA).
#'
#'@param refdna_file refdna_file information to be used to calculate RNA VAF (Default="lib/GRCh37.fa").
#'
#'@param cnv_file A file including copy number variation to calculate cancer cell fraction probability (CCFP) (Default=NA).
#'The format is according to ASCAT (https://www.crick.ac.uk/peter-van-loo/software/ASCAT) output files.
#'The columns are "Chromosome Position Log R segmented LogR BAF segmented BAF Copy number Minor allele Raw copy number"
#'The 1st row should be the above header.
#'See data(sample_copynum); sample_copynum;
#'
#'@param ccfp_dir The file directory to CCFP.pl (Default="lib/ccfp.jar").
#'
#'@param purity Tumor purity or tumor contents ratio required to calculate CCFP (Default=NA).
#'
#'@param netMHCpan_dir The file directory to netMHCpan (Default="lib/netMHCpan-3.0/netMHCpan").
#'
#'@param samtools_dir The file directory to samtools (Default="samtools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF .
#'
#'@param bcftools_dir The file directory to netMHCpan (Default="bcftools").
#'It shouled be indicated when you indicate RNA-bam and try to calculate RNA VAF .
#'samtools 0_x_x includes bcftools in the directory.
#'
#'@return void (Calculated Neoantigen Files will be generated as .tsv files.)
#'
#'@export
MainSNVClass1<-function(input_file, 
                        hla_file, 
                        file_name_in_hla_table = input_file,
                        hmdir = getwd(), 
                        job_id = "NO_job_id", 
                        rnaexp_file = NA, 
                        rnabam_file = NA,
                        cnv_file=NA, 
                        ccfp_dir = paste(hmdir, "lib/ccfp.jar", sep=""), 
                        purity = NA,
                        netMHCpan_dir = paste(hmdir, "lib/netMHCpan-3.0/netMHCIIpan", sep=""),
                        refdna_file = paste(hmdir, "lib/GRCh37.fa", sep=""),
                        refflat_file = paste(hmdir, "lib/refFlat.txt", sep=""),
                        refmrna_file = paste(hmdir, "lib/refMrna.merge.fa", sep=""),
                        samtools_dir = "samtools",
                        bcftools_dir = "bcftools",
                        chr_column = 1, 
                        mutation_start_column = 2,
                        mutation_end_column = 3, 
                        mutation_ref_column = 4, 
                        mutation_alt_column = 5,
                        nm_id_column = 10, 
                        depth_normal_column = NA, 
                        depth_tumor_column = NA,
                        ambiguous_between_exon = 0, 
                        ambiguous_codon = 0,
                        peptide_length = c(8,9,10,11,12,13)){

  if(!file.exists(input_file)) {
    print(paste("Did not find"), input_file)
    return(NULL)
  }
  if(!file.exists(hla_file)) {
    print(paste("Did not find"), hla_file)
    return(NULL)
  }
  if(!file.exists(refflat_file)) {
    print(paste("Did not find"), refflat_file)
    return(NULL)
  }
  if(!file.exists(refmrna_file)) {
    print(paste("Did not find"), refmrna_file)
    return(NULL)
  }

  #Generate FASTA and mutation Profile
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
                     ambiguous_codon = ambiguous_codon)

  output_peptide_txt_file<-paste(input_file, ".", job_id, ".peptide.txt", sep="")
  if(!file.exists(output_peptide_txt_file)){
    return(NULL)
  }

  #Attach RNAseq Data if Exist, Otherwise set NULL column
  skip=FALSE
  if(!is.na(rnaexp_file)){
    if(file.exists(rnaexp_file)){
      GenerateListForGetRNASeq(output_peptide_txt_file, width = 2)
      output_file_rna_list<-paste(output_peptide_txt_file, ".list.txt", sep="")
      print(paste(samtools_dir, "mpileup -l", output_file_rna_list, "-uf", refdna_file, rnabam_file,
                   ">", paste(output_peptide_txt_file, "list.mp", sep=".")))
      error<-tryCatch2(system(paste(samtools_dir, "mpileup -l", output_file_rna_list, "-uf", refdna_file, rnabam_file,
                   ">", paste(output_peptide_txt_file, "list.mp", sep="."))))
      if(error != 0) skip = TRUE

      if(!skip) {
        print(paste(bcftools_dir, "view -c", paste(output_peptide_txt_file, "list.mp", sep="."),
                   ">", paste(output_peptide_txt_file, "list.vcf", sep=".")))
        error<-tryCatch2(system(paste(bcftools_dir, "view -c", paste(output_peptide_txt_file, "list.mp", sep="."),
                   ">", paste(output_peptide_txt_file, "list.vcf", sep="."))))
      }
      if(error != 0) skip = TRUE
      else {
        print("Indicated RNA File does not exist.")
      }
  }
  GetRNAseq(output_peptide_txt_file = output_peptide_txt_file,
            rnaexp_file = rnaexp_file,
            output_file_rna_vcf = paste(output_peptide_txt_file, "list.vcf", sep="."))

  #mutation Rate
  if(ifelse(is.na(ccfp_dir), TRUE, !file.exists(ccfp_dir))) cnv_file<-NA
  if(ifelse(is.na(cnv_file), FALSE, file.exists(cnv_file))){
    GenerateListForCCFP(output_peptide_txt_file, cnv_file = cnv_file, purity = purity)
    if(file.exists(paste(output_peptide_txt_file,".cnv_file.txt",sep=""))){
      print(paste("java -jar ", hmdir, "/", gsub("\\./", "/", ccfp_dir), " ",
                   paste(output_peptide_txt_file, ".cnv_file.txt", sep=""), " > ",
                   paste(output_peptide_txt_file, ".cnv_file.estimate.txt", sep=""), sep=""))
      system(paste("java -jar ", hmdir, "/", gsub("\\./", "/", ccfp_dir), " ",
                   paste(output_peptide_txt_file, ".cnv_file.txt", sep=""), " > ",
                   paste(output_peptide_txt_file, ".cnv_file.estimate.txt", sep=""), sep=""))
      GetRatio(output_peptide_txt_file = output_peptide_txt_file,
               output_peptide_txt_cnc_estimate_file = paste(output_peptide_txt_file,".cnv_file.estimate.txt", sep=""))
    }
  }else{
    GetRatio(output_peptide_txt_file = output_peptide_txt_file,
             output_peptide_txt_cnc_estimate_file = NA)
  }

  #NetMHCpan
  if(is.na(netMHCpan_dir) | !file.exists(netMHCpan_dir)) {
    print(paste("Did not find"), netMHCpan_dir)
    return(NULL)
  }
  print("netMHCpan")
  netMHCpan_script<-scan(netMHCpan_dir, "character", sep="\n", blank.lines.skip = FALSE)
  if(length(grep("# Automatically Overwritten by Neoantimon", netMHCpan_script))==0){
    netMHCpan_par<-gsub("\\./","/", paste(rev(rev(strsplit(netMHCpan_dir, "/")[[1]])[-1]), collapse = "/"))
    if(!file.exists(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))){
      dir.create(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))
    }
    netMHCpan_script<-gsub("#setnv", "setenv", netMHCpan_script)
    
    netMHCpan_script[grep("setenv\tNMHOME", netMHCpan_script)]<-
      paste("setenv\tNMHOME ", getwd(), "/", netMHCpan_par, sep="")
    write.table(c(netMHCpan_script, "# Automatically Overwritten by Neoantimon"),
                netMHCpan_dir, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
  hla<-t(sapply(scan(hla_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  hla_types<-hla[match(file_name_in_hla_table, hla[,1]),-1]
  for(pep in c("peptide", "normpeptide")){
    COUNT<-1
    for(hla_type in hla_types){
      system(paste(netMHCpan_dir,
                   " -l ", paste(peptide_length, collapse = ","),
                   " -f ", paste(input_file, job_id, pep,"fasta",sep="."),
                   " -a HLA-", gsub("\\*","",hla_type),
                   " > ", input_file, ".", job_id, ".HLACLASS1.", COUNT, ".", pep, ".txt",
                   sep=""))
      COUNT <- COUNT + 1
    }
  }
}