#Updated on 2, Feb. 2017. 
------------------------------
#Preparation
------------------------------
(1)Download netMHCpan3.0 and netMHCIIpan 3.1 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan and http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan, respectively. 
Then, set {"netMHCpan-3.0a.Darwin.tar.gz" or "netMHCpan-3.0a.Linux.tar.gz"} and {netMHCIIpan-3.1a.Darwin.tar.gz or netMHCIIpan-3.1a.Linux.tar.gz} in "main". 

(2)Execute the script “Main_Preparatin.sh”. 

(3)See following texts and “SampleRun.R” to learn how to use this library. 

------------------------------
#Functions:
MainSNVClass1(), MainINDELClass1(), MainSNVClass2(), MainINDELClass2()
------------------------------

Arguments:
 (**Required)
 -input_file: An input file annotated by ANNOVAR (http://annovar.openbioinformatics.org/en/latest/) or other softwares to be analyzed. 

 -HLA_file: A tab separated file indicating HLA types. 
 	    The 1st column is input_file name, and the following columns indicate HLA types. 
	    Free to use any Header. 
 	    (e.g., 
	     input_file	A*02:01	A*32:01	B*15:17	B*51:01	C*07:01	C*15:02
	      other_file	A*01:01	A*02:01	B*15:17	B*51:01	C*07:01	C*15:02
	       ...)
			
 (**Not Required)
 (***Standard Profiles)
 -file_name_in_HLA_table (Default=input_file): If the name (1st column) in HLA table is not input_file, indicate the corresponding name. 

 -hmdir (Default=getwd()): Home directory for the analysis. 

 -job_ID (Default="NO_JOB_ID"): Job-Id to be attached in output files. 

 -peptide_length (Default=c(8,9,10,11,12,13)): Peptide Length to be generated. 
 
 -Chr_Column (Default=1): The column number describing Chromosome number in input_file. 
 
 -Mutation_Start_Column (Default=2): The column number describing Mutation Start Position in input_file. 
 
 -Mutation_End_Column (Default=3): The column number describing Mutation End Position in input_file. 
 
 -Mutation_Ref_Column (Default=4): The column number describing Mutation Ref in input_file. 
 
 -Mutation_Alt_Column (Default=5): The column number describing Mutation Alt in input_file. 
 
 -NM_ID_Column (Default=10): The column number describing NM IDs in input_file. 
 
 -ambiguous_between_exon (Default=0): The maximum number to permit the differences between Exon-Lengths from refFlat and refMrna. 
 
 -ambiguous_codon (Default=0):The maximum number to permit the differences between inputfile- and refMrna-oriented translation START/END position. 

 -refFlat_file (Default=paste(hmdir,"/../lib_int/refFlat.txt",sep=""): refFlat file to be used in constructing peptide. This file is automaticalluy generated by Main_Preparation.sh
 
 -refMrna_1 (Default=paste(hmdir,"/../lib_int/refMrna.merge.cut1.fa",sep=""): refMrna file to be used in constructing peptide. This file is automaticalluy generated by Main_Preparation.sh
 
 -refMrna_3 (Default=paste(hmdir,"/../lib_int/refMrna.merge.cut3.fa",sep=""): refMrna file to be used in constructing peptide. This file is automaticalluy generated by Main_Preparation.sh

 (***With RNA-seq Data) 
 -RNAseq_file (Default=NA): A file including RNA expressions. 
 	      The 1st, 2nd and 3rd columns are "GeneSymbol    Chr:ExonStart-ExonEnd(locus)	Expression Amount", respectively.
 	      The 1st row should be any hedear. 
	      (e.g., 
	       header
	        ABR	17:906757-1132315	15.0383 
		 AAGAB	15:67493370-67547533	15.8485
		  ...)

 -RNA_bam (Default=NA): RNA bam file to calculate variant allele frequency of RNA at each mutation. 
 
 -refDNA (Default="./../lib_int/GRCh37.fa"): refDNA information to be used to calculate RNA VAF. 
 
 -Genomon (Default=FALSE): 

 (***With Copy Number Data to Calculate CCFP) 
 -CNV (Default=NA): A file including copy number variation to calculate cancer cell fraction probability (CCFP). 
      The format is according to ASCAT (https://www.crick.ac.uk/peter-van-loo/software/ASCAT) output files. 
      The columns are "Chromosome      Position						      Log R  segmented LogR	BAF segmented BAF Copy number Minor allele Raw copy number"
      The 1st row should be the above hedear. 
      (e.g., 
      Chromosome	Position	Log R	segmented LogR	BAF	segmented BAF	Copy number	Minor allele	Raw copy number
      SNP_A-8575125	1		564621	-1.93486733232252	-0.0986229841341464  1		NA    1		0   -0.497321162431626
      SNP_A-8709646	1		721290	0.0991588550891442	-0.0986229841341464  0		NA    1		0   1.96316993015397
      ...)
 
 -Purity (Default=NA): Tumor purity or tumor contents ratio requierd to calculate CCFP. 
  
------------------------------
#Functions: Summarize The Results
MainMergeClass1(), MainMergeClass2()
------------------------------

Arguments:
 (**Required)
 -input_dir: Directory storing netMHCpan Results.
 
 -input_file_prefix: File prefix of netMHCpan Results. 
  If you have "sample_annovar.txt.NO_JOB_ID.HLACLASS1.1.peptide.txt", please set "sample_annovar", "sample_annovar.txt" and "sample_annovar.txt.NO_JOB_ID". 

 (**Not Required)
 (***Standard Profiles)
 -hmdir (Default=getwd()): Home directory for the analysis. 

 -Tumor_RNA_BASED_ON_DNA (Default=TRUE): In calculating tumor specific RNA expression, TRUE uses variant allele frequency on DNA. Otherwise, use VAF on RNA. 


