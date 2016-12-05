#!/bin/bash
#$ -S /bin/bash
#$ -cwd

R --vanilla --slave --args < Merge_indel.r ${1} ${SGE_TASK_ID} ${SGE_TASK_ID}

