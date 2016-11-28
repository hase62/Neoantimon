#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#INPUT
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

#Echo Basic Information
echo ${2} > ${3}.summary.txt
for key in intronic intergenic exonic ncRNA_intronic ncRNA_exonic downstream upstream UTR3 UTR5 splicing
do
echo ${key} >> ${3}.summary.txt
grep -e $'\t'${key} ${2} | wc >> ${3}.summary.txt
done
echo "exonic - synonymous" >> ${3}.summary.txt
grep -e $'\t'exonic ${2} | grep -e $'\t'synonymous | wc >> ${3}.summary.txt
echo "exonic - nonsynonymous" >> ${3}.summary.txt
grep -e $'\t'exonic ${2} | grep -e $'\t'nonsynonymous | wc >> ${3}.summary.txt
echo "exonic - stopgain" >> ${3}.summary.txt
grep -e $'\t'exonic ${2} | grep -e $'\t'stopgain | wc >> ${3}.summary.txt
echo "exonic - unknown" >> ${3}.summary.txt
grep -e $'\t'exonic ${2} | grep -e $'\t'unknown | wc >> ${3}.summary.txt

#Get HLA-type
echo "grep -w ^${4} ${6}"
TEXT=`grep -w ^${4} ${6}`
IFS=$'\t' read -a HLAs <<< "$TEXT"
sec=`date '+%S%M'`

#Save HLA-type
echo "HLA Type" > ${3}.HLAtype.txt
echo ${TEXT} >> ${3}.HLAtype.txt
echo ${5} >> ${3}.HLAtype.txt

#Generate FASTA and Mutation Profile
if [ -s ${2} ]
then
echo "ok"
else
exit 0
fi

if [ ${#9} -ge 3 ]
then
echo "Not Implemented Yet. "
exit 0
else
cp ${2} ${3}.output.tsv
cp ${2%/*}/*consensus.indel.vcf ${3}.in.vcf
echo "R --vanilla --slave --args < ${1}/GenerateIndelSeq.r ${3}.output.tsv ${3}.in.vcf ${1} 13"
R --vanilla --slave --args < ${1}/GenerateIndelSeq.r ${3}.output.tsv ${3}.in.vcf ${1} 13
fi

#No Indels Detected
if [ -e ${3}.output.tsv.peptide.indel.txt ]
then
   echo "OK"
else
 #NetMHCpan
 COUNT=${#HLAs[@]}
 COUNT=$((COUNT - 1))
 while [ $COUNT -ge 1 ];
 do
  HLA=${HLAs[${COUNT}]}
  echo ${HLA} >> ${3}.HLAtype.txt
  COUNT=$((COUNT - 1))
 done
 exit 0
fi

#Attach RNAseq Data if Exist, otherwise set NULL Column
a="${7##*/}"
if [ ${#a} -ge 8 -a -e ${7} -a -e ${8} ]
then
b="${a%.*}"
echo "${a}"
echo "${b}"
echo "R --vanilla --slave --args < ${1}/GenerateListForGetRNASeq_indel.r ${3}.output.tsv.peptide.indel.txt ${3}.list.txt samtools view -L ${3}.list.txt -h ${8} > ${3}.sam"
R --vanilla --slave --args < ${1}/GenerateListForGetRNASeq_indel.r ${3}.output.tsv.peptide.indel.txt ${3}.list.txt
samtools view -L ${3}.list.txt -h ${8} > ${3}.sam
cut -f3,4,6 ${3}.sam > ${3}.cut.sam.temp
grep -v ^\* ${3}.cut.sam.temp > ${3}.cut.sam 
rm ${3}.sam
rm ${3}.cut.sam.temp
echo "R --vanilla --slave --args < ${1}/GetRNAseq_indel.r ${3}.output.tsv.peptide.indel.txt ${7} ${3}.cut.sam"
R --vanilla --slave --args < ${1}/GetRNAseq_indel.r ${3}.output.tsv.peptide.indel.txt ${7} ${3}.cut.sam
else
echo "R --vanilla --slave --args < ${1}/GetRNAseq_indel.r ${3}.output.tsv.peptide.indel.txt"
R --vanilla --slave --args < ${1}/GetRNAseq_indel.r ${3}.output.tsv.peptide.indel.txt
fi

#Mutation Rate
if [ ${#10} -ge 5 ]
then
mkdir seqware-results
mkdir seqware-results/0
mkdir seqware-results/0/ascat
if [ -e ${10} ]
then
 cp ${10} seqware-results/0/ascat/aaa.copynumber.txt
fi
if [ -e ${10%/*}/*.copynumber.txt ]
then
 cp ${10%/*}/*.copynumber.txt seqware-results/0/ascat/
fi
purity=""
if [ ${#11} -ge 5 ]
then
purity=${11}
fi
fi
if [ -e seqware-results/0/ascat/*.copynumber.txt ]
then
echo "R --vanilla --slave --args < ${1}/GenerateListForCCFP.r ${3}.output.tsv.peptide.indel.txt ${purity}"
R --vanilla --slave --args < ${1}/GenerateListForCCFP.r ${3}.output.tsv.peptide.indel.txt ${purity}
echo "perl ${1}/MUT/perl/CCFP.pl ${3}.output.tsv.peptide.indel.txt.cnv.txt > ${3}.output.tsv.peptide.indel.txt.cnv.estimate.txt"
perl ${1}/MUT/perl/CCFP.pl ${3}.output.tsv.peptide.indel.txt.cnv.txt > ${3}.output.tsv.peptide.indel.txt.cnv.estimate.txt
echo "R --vanilla --slave --args < ${1}/GetRatio.r ${3}.output.tsv.peptide.indel.txt ${3}.output.tsv.peptide.indel.txt.cnv.estimate.txt"
R --vanilla --slave --args < ${1}/GetRatio.r ${3}.output.tsv.peptide.indel.txt ${3}.output.tsv.peptide.indel.txt.cnv.estimate.txt
else
echo "Dose not exist ${10}"
echo "R --vanilla --slave --args < ${1}/GetRatio.r ${3}.output.tsv.peptide.indel.txt"
R --vanilla --slave --args < ${1}/GetRatio.r ${3}.output.tsv.peptide.indel.txt
fi

#NetMHCpan
j=peptide.indel
COUNT=0
COUNT1=${#HLAs[@]}
COUNT1=$((COUNT1 - 1))
while [ $COUNT1 -ge 1 ];
do
 HLA=${HLAs[${COUNT1}]} 
 echo ${HLA} >> ${3}.HLAtype.txt
 if [ `echo ${HLA}|grep DRB1` ]; then
  COUNT=`expr ${COUNT} + 1`
  tmp=${1}/netMHCIIpan-3.1/tmp/$RANDOM$RANDOM$RANDOM$RANDOM
  mkdir $tmp
  echo "" > ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt
  echo "#!/bin/bash" > Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  echo "#$ -S /bin/bash" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  echo "#$ -cwd" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  echo "" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  echo "${1}/netMHCIIpan-3.1/netMHCIIpan -length 15 -f ${3}.output.tsv.${j}.fasta -a ${HLA} -tdir ${tmp} >> ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  sed -i -e "s/\*/_/g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  sed -i -e "s/\://g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
  qsub ./Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
 elif [ `echo ${HLA}|grep DPA1` ]; then
  COUNT2=${#HLAs[@]}
  COUNT2=$((COUNT2 - 1))
  while [ $COUNT2 -ge 1 ];
  do
   HLA2=${HLAs[${COUNT2}]}
   if [ `echo ${HLA2}|grep DPB1` ]; then
    COUNT=`expr ${COUNT} + 1`
    tmp=${1}/netMHCIIpan-3.1/tmp/$RANDOM$RANDOM$RANDOM$RANDOM
    mkdir $tmp
    echo "" > ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt
    echo "#!/bin/bash" > Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "#$ -S /bin/bash" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "#$ -cwd" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "${1}/netMHCIIpan-3.1/netMHCIIpan -length 15 -f ${3}.output.tsv.${j}.fasta -choose -cha ${HLA} -chb ${HLA2} -tdir ${tmp} >> ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    sed -i -e "s/\*//g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    sed -i -e "s/\://g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    qsub ./Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
   fi
   COUNT2=$((COUNT2 - 1))
  done
 elif [ `echo ${HLA}|grep DQA1` ]; then
  COUNT2=${#HLAs[@]}
  COUNT2=$((COUNT2 - 1))
  while [ $COUNT2 -ge 1 ];
  do
   HLA2=${HLAs[${COUNT2}]}
   if [ `echo ${HLA2}|grep DQB1` ]; then
    COUNT=`expr ${COUNT} + 1`
    tmp=${1}/netMHCIIpan-3.1/tmp/$RANDOM$RANDOM$RANDOM$RANDOM
    mkdir $tmp
    echo "" > ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt
    echo "#!/bin/bash" > Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "#$ -S /bin/bash" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "#$ -cwd" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    echo "${1}/netMHCIIpan-3.1/netMHCIIpan -length 15 -f ${3}.output.tsv.${j}.fasta -choose -cha ${HLA} -chb ${HLA2} -tdir ${tmp} >> ${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.txt" >> Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    sed -i -e "s/\*//g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    sed -i -e "s/\://g" Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh
    qsub ./Exe.${3}.peptide.fasta.HLACLASS2.${COUNT}.${j}.sh  
   fi
   COUNT2=$((COUNT2 - 1))
  done
 fi
 COUNT1=$((COUNT1 - 1))
done

