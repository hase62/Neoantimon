#!/bin/bash
#$ -S /bin/bash
#$ -cwd

for i in `ls ${1}/*/*.peptide.fasta.HLA.[0-9].*.extracted.txt.*peptide.txt `
do
if [ -s ${i} ]; then
:
else
#if test 40 -le ${#i}; then
echo ${i}
#mv "${i%/*}" tmp2/NetMHCPan_Results
#rm -r ${i%/*}
#mv ${i%/*} tmp/${i%/*}
fi
done

for i in `ls ${1}/*/*.peptide.fasta.HLA*.[0-9].peptide.indel.txt `
do
if [ -s ${i} ]; then
:
else
echo ${i}
#mv "${i%/*}" tmp2/NetMHCPan_Results
#rm -r ${i%/*}
fi
done
