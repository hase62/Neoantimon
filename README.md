#Updated on 22th, Feb. 2017. 
==============================
##Preparation
------------------------------
*Download netMHCpan*
Download netMHCpan3.0 and netMHCIIpan 3.1 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan and http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan, respectively. Then, prepare to use netMHCpan/netMHCIIpan by extracting *tar.gz. 
```
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
```

(2) 
install.packages("devtools")
library(devtools)
install_github('hase62/Neoantimon')
library(Neoantimon)
vignette("SampleCodeForNeoantimon")

(3)See following texts and “SampleRun.R” to learn how to use this library. 
