#!/bin/bash
#$ -S /bin/bash
#$ -cwd

R --vanilla --slave --args < /share/icgc_pan/WORK/hasegawa/NeoAntimon/StatSnvMatrix.class1.r

