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
my $nonsilent = 0;
my %bgGene;

my %trial;        #  $trial{$gene}
my $binomialParameter;
my %success;      # $success{$gene}
my %pvalue; #$pvalue{$gene}

my $home = get_home_dir();

my $usage = qq(usage: $0 [options] mutfile
        -c recurence_cutoff
        -m (remove_multiple_mutations_within_a_sample)
        -n (count_nonsilent_mutaions)
        -b background_gene (gene1:gene2....)
);

my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-m//){
    $rmMultiMut = 1;
}
if($argv =~ s/-n//){
    $nonsilent  = 1;
}
if($argv =~ s/-b\s+(\S+)//){
    map {$bgGene{$_}=1} split(":",$1);
}
@ARGV = split(" ", $argv);

my $mutFile = shift or die $usage;

###############
#read mutFile
###############

open(IN,$mutFile);

my %seen;
while(<IN>){
    chomp;    
    my @tmp =  split("\t");
    my $sample = $tmp[0];
    my $type = $tmp[6];
    $type or next;
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    my $gene = $tmp[1];
    $gene{$gene}=1;
    if($rmMultiMut and $seen{$sample}{$gene}){
	next;
    }
    $seen{$sample}{$gene}=1;
    $trial{$gene}++;
    if($nonsilent?isNonsilent($type):isDamaging($type)){
	$success{$gene}++
    }
}

@gene = sort keys %gene;
@recurrentMutatedGene = grep {$trial{$_} and $trial{$_} >= $recurrentCutoff} @gene;
@recurrentMutatedGene = sort {$trial{$b} <=> $trial{$a}} @recurrentMutatedGene;
map {$recurrentMutatedGene{$_}=1} @recurrentMutatedGene;
map {$success{$_}?():($success{$_}=0)} @gene;
map {$trial{$_}?():($trial{$_}=0)} @gene;    

################################
#get binomial parameter
################################

my @bgGene = %bgGene?(keys %bgGene):(grep {$trial{$_} == 1} @gene);
my $trialSum;
map {$trialSum += $trial{$_}} @bgGene;
my $successSum;
map {$successSum += $success{$_}}   @bgGene;
$binomialParameter = $successSum/$trialSum;

warn "binomialParameter: $binomialParameter\n";

####################
#print R input file
####################

my $trialFile = "tmp${$}.trial.tsv";
my $successFile = "tmp${$}.success.tsv";
my $geneFile = "tmp${$}.gene.txt";


open(OUT, ">$trialFile");
print OUT join("\n",map {$trial{$_}} @recurrentMutatedGene)."\n";
close(OUT);

open(OUT, ">$successFile");
print OUT join("\n",map {$success{$_}} @recurrentMutatedGene)."\n";
close(OUT);

open(OUT, ">$geneFile");
print OUT join("\n",@recurrentMutatedGene)."\n";
close(OUT);


####################
#print R script file
####################

my $RscriptFile = "tmp${$}.script.R";
my $RoutFile = "tmp${$}.out.txt";

open(OUT, ">$RscriptFile");
print OUT "S <- scan(\"$successFile\")\n";
print OUT"T <- scan(\"$trialFile\")\n";
print OUT"p <- $binomialParameter\n";
print OUT<<"EOF";
pvalue<-1-pbinom(S-1,T, p)
pvalue<--log10(pvalue)
pvalue[pvalue>100]<-100
EOF
print OUT "write(pvalue,\"$RoutFile\",sep =\"\\n\")\n"; 

####################
#run R script
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
    print $_."\t".$pvalue{$_}."\t".$success{$_}."/".$trial{$_}."\n";
}


