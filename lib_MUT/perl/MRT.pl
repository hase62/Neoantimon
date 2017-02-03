#!/usr/local/bin/perl 

use strict;
use warnings;


$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();

###############
#set variables
###############

my @gene;
my %gene;
my @recurrentMutatedGene;
my %recurrentMutatedGene;

my $recurrentCutoff = 3;
my $rmMultiMut = 0;

my %geneLength;
my $geneLengthSum;

my $mutSum;	
my %mutCount; # $mutCounnt{$gene}
my %poissonParameter; # $poissonParameter{$gene}
my %pvalue; #$pvalue{$gene}

my $home = get_home_dir();

my $usage = "usage: $0 [-c recurence_cutoff -m (remove_multiple_mutations_within_a_sample)] mutfile\n";
my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-m//){
    $rmMultiMut  = 1;
}
@ARGV = split(" ", $argv);

my $mutFile = shift or die $usage;
my $geneLengthFile = "${home}/data/cdsLength.tsv";

###############
#read mutFile
###############

open(IN,$mutFile);

my  %seen; 
while(<IN>){
    chomp;    
    my @tmp =  split("\t");
    my $sample = $tmp[0];
    my $type = $tmp[6];
    $type or next;
    isNonsilent($type) or  next;
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    my $gene = $tmp[1];
    if($rmMultiMut and $seen{$sample}{$gene}){
	next;
    }
    $seen{$sample}{$gene}=1;
    $mutCount{$gene}++;  
    $mutSum++;
    $gene{$gene}=1;
}


@gene = sort keys %gene;
@recurrentMutatedGene = grep {$mutCount{$_} and $mutCount{$_} >= $recurrentCutoff} @gene; 
@recurrentMutatedGene = sort {$mutCount{$b} <=> $mutCount{$a}} @recurrentMutatedGene;  
map {$recurrentMutatedGene{$_}=1} @recurrentMutatedGene;


#################
#set gene length
#################

open(IN, $geneLengthFile);
while(<IN>){
    chomp;
    my @tmp =  split("\t");
    $geneLength{$tmp[0]} = $tmp[1];
    $geneLengthSum +=  $tmp[1];
}


##############################
#get  poisoon parameter
##############################

@recurrentMutatedGene = grep {$geneLength{$_}} @recurrentMutatedGene;
map {$poissonParameter{$_} = $mutSum*$geneLength{$_}/$geneLengthSum}  @recurrentMutatedGene;

####################
#print R input file
####################

my $parameterFile = "tmp${$}.para.tsv";
my $statisticsFile = "tmp${$}.stat.tsv";

open(OUT, ">$parameterFile");
print OUT join("\n",map {$poissonParameter{$_}} @recurrentMutatedGene)."\n";
close(OUT);

open(OUT, ">$statisticsFile");
print OUT join("\n",map {$mutCount{$_}} @recurrentMutatedGene)."\n";
close(OUT);


####################
#print R script file
####################

my $RscriptFile = "tmp${$}.script.R";
my $RoutFile = "tmp${$}.out.txt";

open(OUT, ">$RscriptFile");
print OUT  "P <- scan(\"$parameterFile\")\n";
print OUT "S <- scan(\"$statisticsFile\")\n";
print OUT<<"EOF";
pvalue<-1-ppois(S-1, P)
pvalue<--log10(pvalue)
pvalue[pvalue>100]<-100
EOF
print OUT "write(pvalue,\"$RoutFile\",sep =\"\\n\")\n"; 

####################
#get P value
####################

`R --vanilla < $RscriptFile >& /dev/null`;
open(IN, $RoutFile);
chomp(my @pvalue = <IN>);
for(my $i=0; $i<@recurrentMutatedGene; $i++){
     $pvalue{$recurrentMutatedGene[$i]}=$pvalue[$i];
}

system("rm tmp${$}.*");

####################
#print result
####################

foreach(sort {$pvalue{$b} <=> $pvalue{$a}} @recurrentMutatedGene){
    print $_."\t".$pvalue{$_}."\t".$mutCount{$_}."/".$geneLength{$_}."\n";
}


