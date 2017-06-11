#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#At first, Download netMHCpan-3.0a.Linux.tar.gz (LINUX) or netMHCpan-3.0a.Darwin.tar.gz (MAC) from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan

#Install netMHCpan3.0
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
