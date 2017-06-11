#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#At first, download netMHCIIpan-3.1a.Linux.tar (LINUX) or netMHCIIpan-3.1a.Darwin.tar.gz (MAC) from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan

#Install netMHCIIpan3.1
if [ -e netMHCIIpan-3.1a.Linux.tar.gz ]; then
 tar zxvf netMHCIIpan-3.1a.Linux.tar.gz
elif [ -e netMHCIIpan-3.1a.Darwin.tar.gz ]
 tar zxvf netMHCIIpan-3.1a.Darwin.tar.gz
fi
cd netMHCIIpan-3.1
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/packages/netMHCIIpan/3.1/netMHCIIpan-3.1:${cdir}:g" netMHCIIpan
sed -i -e "15s:^:setenv  TMPDIR \$\{NMHOME\}/tmp:" netMHCIIpan
wget http://www.cbs.dtu.dk/services/NetMHCIIpan-3.1/data.tar.gz
#gunzip -c data.tar.gz | tar xvf -
tar -xvf data.tar.gz
