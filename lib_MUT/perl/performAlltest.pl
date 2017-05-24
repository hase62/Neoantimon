#!/usr/local/bin/perl                                                                                  
use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();


my $home = get_home_dir();

my $geneset = "$home/data/curatedCp.gmt";
my $PPI = "$home/data/HPRD_BIND_BioGrid_IntAct_MINT.human.sif";
my $cutoff = 5;
my $rmSignificantGenesInGsTest = 0;
my $HETcutoff = 2;
my $skipHST = 0;

my $usage = qq(usage: $0 [options] mutfile
    -g <arg> geneset file (default: $geneset)
    -p <arg> PPI file  (default: $PPI)
    -c <arg> count cutoff  (default: $cutoff)
    -r <arg> pvalue cutoff for removing significant genes in geneset test
    -h <arg> HET pvalue cutoff (default: $HETcutoff)
    -s       skip HST
);

my $argv = join(" ", @ARGV);
if($argv =~ s/-g\s+(\S+)//){
    $geneset = $1;
}
if($argv =~ s/-p\s+(\S+)//){
    $PPI = $1;
}
if($argv =~ s/-c\s+(\S+)//){
    $cutoff = $1;
}
if($argv =~ s/-r\s+(\S+)//){
    $rmSignificantGenesInGsTest = $1;
}
if($argv =~ s/-s//){
    $skipHST = 1;
}
@ARGV = split(" ", $argv);


my $infile = shift or die $usage;
my $suffix = get_filename_without_suffix($infile);

my $outDir = $suffix;
unless(-d $outDir){
    `mkdir  $outDir`;
}

{
    my @script = qw(MRT.pl DMET.pl getNSmutCount.pl HST.pl);
    if($skipHST){
        @script = qw(MRT.pl DMET.pl getNSmutCount.pl);
    }
    my @outfile = qw(MRT.txt  DMET.txt Count.txt HST.txt);

    for(my $i=0; $i<@script; $i++){
        my $command = "perl $home/perl/$script[$i] -c $cutoff   $infile > $outDir/$outfile[$i]";
	warn($command."\n");
	system($command);
    }
    my $command = "perl $home/perl/combineResultsAsTsv.pl ".join(" ", map { "$outDir/$_" } @outfile)." > $outDir/geneTest.tsv";
    warn($command."\n");
    system($command);
}


{
    my @script = qw(MRTdm.pl DMETdm.pl);
    my @outfile = qw(MRTdm.txt  DMETdm.txt);

    for(my $i=0; $i<@script; $i++){
        my $command = "perl $home/perl/$script[$i] -c $cutoff   $infile > $outDir/$outfile[$i]";
        warn($command."\n");
	system($command);
    }
    my $command = "perl $home/perl/combineResultsAsTsv.pl  ".join(" ", map { "$outDir/$_" } @outfile)." > $outDir/dmTest.tsv";
    warn($command."\n");
    system($command);
}



{
    my @script = qw(MRTgs.pl DMETgs.pl);
    my @outfile = qw(MRTgs.txt DMETgs.txt);
    
    if($rmSignificantGenesInGsTest){
	my @significantGenes;
	foreach my $f ("$outDir/MRT.txt", "$outDir/DMET.txt"){
	    open(IN, "$f");
	    my @gene;
	    while(<IN>){
		chomp;
		my @tmp = split("\t");
		if($tmp[1] > $rmSignificantGenesInGsTest){
		    push(@gene, $tmp[0]);
		}
	    }
	    push(@significantGenes, join(":", @gene));
	}
	for(my $i = 0; $i < 2; $i++){
	    if($significantGenes[$i]){
		$script[$i] .= " -r ".$significantGenes[$i];
	    }
	}
    }
    
    for(my $i=0; $i<@script; $i++){
        my $command = "perl $home/perl/$script[$i] -c $cutoff  $infile $geneset> $outDir/$outfile[$i]";
        warn($command."\n");
	system($command);
    }
    my $command = "perl $home/perl/combineResultsAsTsv.pl  ".join(" ", map { "$outDir/$_" } @outfile)." > $outDir/gsTest.tsv";
    warn($command."\n");
    system($command);
}



my $command = "perl $home/perl/HET.pl -c $HETcutoff  $outDir/MRT.txt $PPI > $outDir/HET.txt";
warn($command."\n");
system($command);

$command = "perl $home/perl/getCytoscapeInputFileFromHET.pl $outDir/HET.txt $outDir/HET";
warn($command."\n");
system($command);

$command = "perl perl/getNSmutMat.pl -c $cutoff $infile > $outDir/NSmutMat.tab";
warn($command."\n");
system($command);

