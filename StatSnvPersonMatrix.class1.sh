#!/bin/bash
#$ -S /bin/bash
#$ -cwd

R --vanilla --slave --args < StatSnvPersonMatrix.class1.r ${SGE_TASK_ID} ${1}
