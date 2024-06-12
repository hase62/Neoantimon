#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCIIpan-4.3e.Linux.tar.gz -a ! -e netMHCIIpan-4.3 ]; then
 tar zxvf netMHCIIpan-4.3e.Linux.tar.gz
elif [ -e netMHCIIpan-4.3e.Darwin_x86_64.tar.gz -a ! -e netMHCIIpan-4.3 ]; then
 tar zxvf netMHCIIpan-4.3e.Darwin_x86_64.tar.gz
elif [ -e netMHCIIpan-4.3e.Darwin_arm64.tar.gz -a ! -e netMHCIIpan-4.3 ]; then
 tar zxvf netMHCIIpan-4.3e.Darwin_arm64.tar.gz
fi
cd netMHCIIpan-4.3
mkdir tmp
cdir=`pwd`
sed -i -e "s:/tools/src/netMHCIIpan-4.3:${cdir}:g" netMHCIIpan
sed -i -e "18s:^:setenv  TMPDIR \$\{NMHOME\}/tmp:" netMHCIIpan

