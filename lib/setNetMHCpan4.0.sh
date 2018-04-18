#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCpan-4.0a.Linux.tar.gz ]; then
 tar zxvf netMHCpan-4.0a.Linux.tar.gz
fi
if [ -e netMHCpan-4.0a.Darwin.tar.gz ]; then
 tar zxvf netMHCpan-4.0a.Darwin.tar.gz
fi
cd netMHCpan-4.0
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/packages/netMHCpan/4.0/netMHCpan-4.0:${cdir}:g" netMHCpan
sed -i -e "s:#setenv:setenv:g" netMHCpan
wget http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.tar.gz
tar -xvf data.tar.gz
