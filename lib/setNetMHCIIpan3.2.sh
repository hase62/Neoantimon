#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCIIpan-3.2.Linux.tar.gz ]; then
 tar zxvf netMHCIIpan-3.2.Linux.tar.gz
fi
if [ -e netMHCIIpan-3.2.Darwin.tar.gz ]; then
 tar zxvf netMHCIIpan-3.2.Darwin.tar.gz
fi
cd netMHCIIpan-3.2
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/bio/src/netMHCIIpan-3.2:${cdir}:g" netMHCIIpan
sed -i -e "15s:^:setenv  TMPDIR \$\{NMHOME\}/tmp:" netMHCIIpan
if [ -e ./../netMHCIIpan-3.2.Linux.tar.gz ]; then
 wget http://www.cbs.dtu.dk/services/NetMHCIIpan-3.2/data.Linux.tar.gz
 tar -xvf data.Linux.tar.gz
fi
if [ -e ./../netMHIICpan-3.2.Darwin.tar.gz ]; then
 wget http://www.cbs.dtu.dk/services/NetMHCIIpan-3.2/data.Darwin.tar.gz
 tar -xvf data.Darwin.tar.gz
fi
