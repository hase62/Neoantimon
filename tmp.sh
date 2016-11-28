#!/bin/bash
#$ -S /bin/bash
#$ -cwd


#echo ${TEXT}
#IFS=$'\t' read -r ID <<< "$TEXT"

TEXT=`cat ${1}`
echo ${TEXT}
IFS=$'\t' read -a HLAs <<< "$TEXT"

#NetMHCpan
COUNT=${#HLAs[@]}
COUNT=$((COUNT - 1))
while [ $COUNT -ge 1 ];
do
 HLA=${HLAs[${COUNT}]}
 #echo ${HLA} >> ${3}.HLAtype.txt
 if [ -n "${HLA}" ]
 then
  if [ ${#HLA} -ge 10 ]
  then
   COUNT=$((COUNT - 1))
   continue
  fi
  echo ${HLA}
  for j in peptide normpeptide
  do
   echo peptide.fasta.HLACLASS1.${COUNT}.${j}.sh
  done
 fi
COUNT=$((COUNT - 1))
done
