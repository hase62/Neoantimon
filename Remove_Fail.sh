#!/bin/bash
#$ -S /bin/bash
#$ -cwd

for i in `grep token ${1}/*/Exe_*.sh.e*  | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep unexpected ${1}/*/Exe_*.sh.e*  | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep Exception ${1}/*/Exe_*.sh.e*  | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep while ${1}/*/Exe_*.sh.e*  | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep Unable ${1}/*/Exe_*.sh.e*  | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep \: ${1}/*/Exe_*.sh.e* | grep -v JAVA | grep -v Warning | grep -v 1\: | grep -v pkill | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

for i in `grep Can ${1}/*/Exe.*.sh.e* | cut -f1 -d ":"`
do
if test 40 -le ${#i};
then
echo ${i}
echo ${i%/*}
rm -r ${i%/*}
fi
done

