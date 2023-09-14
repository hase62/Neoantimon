#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCIIpan-4.2b.Linux.tar.gz ]; then
 tar zxvf netMHCIIpan-4.2b.Linux.tar.gz
fi
if [ -e netMHCIIpan-4.2b.Darwin.tar.gz ]; then
 tar zxvf netMHCIIpan-4.2b.Darwin.tar.gz
fi
cd netMHCIIpan-4.2
mkdir tmp
cdir=`pwd`
sed -i -e "s:/tools/src/netMHCIIpan-4.2:${cdir}:g" netMHCIIpan
sed -i -e "15s:^:setenv  TMPDIR \$\{NMHOME\}/tmp:" netMHCIIpan

