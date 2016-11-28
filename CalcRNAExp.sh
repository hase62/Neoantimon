#!/bin/bash
#$ -S /bin/bash
#$ -cwd

R --vanilla --slave --args < CalcRNAExp.r ${SGE_TASK_ID} ${1}


