#!/bin/bash
#$ -S /bin/bash
#$ -cwd

for i in `seq ${1} ${2}`
do
#echo ${i}
qdel ${i}
done

