#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCpan-3.0a.Linux.tar.gz ]; then
 tar zxvf netMHCpan-3.0a.Linux.tar.gz
fi
if [ -e netMHCpan-3.0a.Darwin.tar.gz ]; then
 tar zxvf netMHCpan-3.0a.Darwin.tar.gz
fi
cd netMHCpan-3.0
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/packages/netMHCpan/3.0/netMHCpan-3.0:${cdir}:g" netMHCpan
sed -i -e "s:#setenv:setenv:g" netMHCpan
wget http://www.cbs.dtu.dk/services/NetMHCpan-3.0/data.tar.gz
tar -xvf data.tar.gz
