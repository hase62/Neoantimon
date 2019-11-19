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
wget http://www.cbs.dtu.dk/services/NetMHCIIpan-3.2/data.tar.gz
tar -xvf data.tar.gz
