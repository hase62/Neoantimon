#!/bin/bash
#$ -S /bin/bash
#$ -cwd

COUNT=1
MAX_COUNT=1200
HOME=/share/icgc_pan/WORK/hasegawa/NeoAntimon
OUTPUT_DIR=/share/icgc_pan/WORK/hasegawa/NeoAntimon/Output.Indel.1

#
#This script require only one file including these information in each (1-7)column
#

#Format of An Input Table
#1: SampleID:   (ex.) EsoCRT_01A
#2: InputFile:  (ex.) /home/mimorin/Hirata/Neoantigen/filter.mutation.20160325/EsoCRT_01A_tumor_genomon_mutations.result.txt
#3: HLATable:   (ex.) /home/mimorin/Hirata/Neoantigen/EsoCRT.HLAGenotype.20160325.txt
#4: RNASeq:     (ex.) - / /home/mimorin/Hirata/Neoantigen/RNAseq.20160325/EsoCRT_01A.tab
#5: RNAbam:     (ex.) - / /home/mimorin/RNAseq/Tophat/EsoCRT_01A/sort_accepted_hits.bam
#6: Ascat:      (ex.) - / <fullpath>
#7: Genomon:    (ex.) - / Genomon
#8: CorrespondingNameInHLAtable:   (ex.) EsoCRT_01A_1, if you require. 
#9: Annotation: (ex.) ThisSampleIsAAA, if you require. 
#10: Purity
cat ${1} | while read line
do
  IFS=$'\t' read -r sample main hla rna rna_bam ascat genomon correspo annotation purity <<< "$line"
  if [ $MAX_COUNT -le $COUNT ] ; then
    break
  fi
  if [ -e ${OUTPUT_DIR}/${sample} ]; then
     echo "DIR. EXISTS ${sample}"
  else
     mkdir ${OUTPUT_DIR}/${sample}
     cd ${OUTPUT_DIR}/${sample}

     #CorrespondingNameInHLAtable
     if [ ${#correspo} -lt 1 ]
     then
     correspo=${sample}
     fi     
     #Annotation
     if [ ${#annotation} -lt 1 ]
     then
     annotation=${sample}
     fi
     #1:  Home Directory
     #2:  Annovar Output
     #3:  Job id
     #4:  Names in HLA-Table
     #5:  Annotation
     #6:  HLA-Table
     #7:  TopHat Output
     #8:  RNA BAM-File
     #9:  Genomon
     #10: Ascat
     #11: Purity (if required)
     echo "./../../Exe_indel_class1.sh ${HOME} ${main} ${sample} ${correspo} ${annotation} ${hla} ${rna} ${rna_bam} ${genomon} ${ascat} ${purity}"
     qsub -l ljob ./../../Exe_indel_class1.sh ${HOME} ${main} ${sample} ${correspo} ${annotation} ${hla} ${rna} ${rna_bam} ${genomon} ${ascat} ${purity}
     cd ${OUTPUT_DIR}
     COUNT=`expr $COUNT + 1`
  fi
  echo "---"
done
cd ${HOME}
