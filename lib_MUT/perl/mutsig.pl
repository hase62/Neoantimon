#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and push(@INC, $&);
require MUT;

my $usage = "usage: $0 mutations.maf coverage.txt covariates.txt output.txt \n";
@ARGV == 4 or die $usage;
my $tmpFilePrefix = "tmp${$}"; 
chomp(my $pwd = `pwd`);

my $sgeMem = 8;
my $sgeMemArg = "";
if($sgeMem > 2){
    $sgeMemArg = "-l s_vmem=${sgeMem}G -l mem_req=${sgeMem}";
}

my $mutsig = get_home_dir()."/mutsig/MutSigCV_1.4/run_MutSigCV.sh";
my $matlab = get_home_dir()."/mutsig/MATLAB_Compiler_Runtime/v81";

my $command = join(" ",($mutsig, $matlab, @ARGV));
print_SGE_script("$command", "${tmpFilePrefix}.pl");
while(system("qsub $sgeMemArg  -b y -cwd  -o ${tmpFilePrefix}.out  -e ${tmpFilePrefix}.err   $pwd/${tmpFilePrefix}.pl")!=0){
	sleep(10);
}
wait_for_SGE_finishing("${tmpFilePrefix}");
#`rm -r ${tmpFilePrefix}.*`;
