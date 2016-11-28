#!/bin/bash
#$ -S /bin/bash
#$ -cwd

R --vanilla --slave --args < StatSnvPersonMatrix.class2.r ${SGE_TASK_ID} ${1}

