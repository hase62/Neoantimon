#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#e.g., qsub -t 1-20 Merger_snv.sh Output
R --vanilla --slave --args < Merge_snv.class1.r ${1} ${SGE_TASK_ID} ${SGE_TASK_ID}
#R --vanilla --slave --args < Merge_snv.class1.r ${1} ${2} ${3}

