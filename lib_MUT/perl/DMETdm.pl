#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();

###############
#set variables
###############

my @domain;
my %domain;
my @recurrentMutatedDomain;
my %recurrentMutatedDomain;
my $recurrentCutoff = 3;
my $rmMultiMut = 0;
my $nonsilent = 0;
my %bgGene;

my @gene;
my %gene;
my @recurrentMutatedGene;
my %recurrentMutatedGene;


my %trial;        #  $trial2{$gene}
my %success;      # $success2{$gene}

my %trial2;        #  $trial2{$domain}
my $binomialParameter;
my %success2;      # $success2{$domain}

my %pvalue; #$pvalue{$domain}111
my %dm2gene;

my $home = get_home_dir();


my $usage = qq(usage: $0 [options] mutfile gsfile
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
    $rmMultiMut  = 1;
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
    my $gene = $tmp[1] or next;
    if($rmMultiMut and $seen{$sample}{$gene}){
        next;
    }
    $seen{$sample}{$gene}=1;
    $gene{$gene}=1;
    $trial{$gene}++;
    if($nonsilent?isNonsilent($type):isDamaging($type)){
        $success{$gene}++;
    }
    $tmp[3] or next;
    my @dm = split(",", $tmp[3]);
    map {$trial2{$_}++} @dm;
    map {$domain{$_}=1} @dm;
    if($nonsilent?isNonsilent($type):isDamaging($type)){
	map {$success2{$_}++} @dm;
	foreach(@dm){
	    $dm2gene{$_}{$tmp[1]}++;
	} 
    }
}

@domain = sort keys %domain;
@gene = sort keys %gene;
@recurrentMutatedDomain = grep {$trial2{$_} and $trial2{$_} >= $recurrentCutoff} @domain;
@recurrentMutatedDomain = sort {$trial2{$b} <=> $trial2{$a}} @recurrentMutatedDomain;
map {$recurrentMutatedDomain{$_}=1} @recurrentMutatedDomain;
map {$success2{$_}?():($success2{$_}=0)} @domain;
map {$trial2{$_}?():($trial2{$_}=0)} @domain;
map {$success{$_}?():($success{$_}=0)} @gene;
map {$trial{$_}?():($trial{$_}=0)} @gene;

################################
#get binomial parameter
################################

my @bgGene = %bgGene?(keys %bgGene):(grep {$trial{$_} == 1} @gene);
my $trialSum;
map {$trialSum += $trial{$_}}  @bgGene;
my $successSum;
map {$successSum += $success{$_}}  @bgGene;
$binomialParameter = $successSum/$trialSum;

warn "binomialParameter: $binomialParameter\n";

####################
#print R input file
####################

my $trialFile = "tmp${$}.trial.tsv";
my $successFile = "tmp${$}.success2.tsv";
my $domainFile = "tmp${$}.gene.txt";


open(OUT, ">$trialFile");
print OUT join("\n",map {$trial2{$_}} @recurrentMutatedDomain)."\n";
close(OUT);

open(OUT, ">$successFile");
print OUT join("\n",map {$success2{$_}} @recurrentMutatedDomain)."\n";
close(OUT);

open(OUT, ">$domainFile");
print OUT join("\n",@recurrentMutatedDomain)."\n";
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

for(my $i=0; $i<@recurrentMutatedDomain; $i++){
    $pvalue{$recurrentMutatedDomain[$i]}=$pvalue[$i];
}

system("rm tmp${$}.*");

####################
#print result
####################

foreach my  $d (sort {$pvalue{$b} <=> $pvalue{$a}} @recurrentMutatedDomain){
    print $d."\t".$pvalue{$d}."\t".$success2{$d}."/".$trial2{$d}."\t";
    my @tmp =  grep {$dm2gene{$d}{$_}} keys %{$dm2gene{$d}};
    @tmp = sort {$dm2gene{$d}{$b} <=> $dm2gene{$d}{$a}} @tmp;
    print join(",", map {$_."(".$dm2gene{$d}{$_}.")"} @tmp)."\n";
}

