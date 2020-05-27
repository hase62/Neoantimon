#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCpan-4.1.Linux.tar.gz ]; then
 tar zxvf netMHCpan-4.1.Linux.tar.gz
fi
if [ -e netMHCpan-4.1.Darwin.tar.gz ]; then
 tar zxvf netMHCpan-4.1.Darwin.tar.gz
fi
cd netMHCpan-4.1
mkdir tmp
cdir=`pwd`
sed -i -e "s:/net/sund-nas.win.dtu.dk/storage/services/www/packages/netMHCpan/4.1/netMHCpan-4.1:${cdir}:g" netMHCpan
sed -i -e "s: /tmp:${cdir}/tmp:g" netMHCpan

wget --no-check-certificate https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz
tar -xvf data.tar.gz
