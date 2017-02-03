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
my @recurrentMutatedGs;
my %recurrentMutatedGs;
my $recurrentCutoff = 10;
my $rmMultiMut = 0;

my @gs;
my %gene2gs;
my %gs2gene;

my %geneLength;
my %gsLength;
my $geneLengthSum;

my %removedGene;

my $mutSum;	
my %mutCount;   # $mutCounnt{$gene}
my %mutCount2;  # $mutCounnt{$gs}
my %poissonParameter; # $poissonParameter{$gene}
my %pvalue;     # $pvalue{$gene}

my $home = get_home_dir();

my $usage = "usage: $0 [-c recurence_cutoff -r removed_gene (gene1:gene2...) -m (remove_multiple_mutations_within_a_sample)] mutfile gsFile\n";
my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-r\s+(\S+)//){
    my @tmp  = split(/[^\w]/,$1);
    map {$removedGene{$_} = 1} @tmp;
}
if($argv =~ s/-m//){
    $rmMultiMut  = 1;
}
@ARGV = split(" ", $argv);

my $mutFile = shift or die $usage;
my $gsFile = shift or die $usage;
my $geneLengthFile = "${home}/data/cdsLength.tsv";

###############
#read gsFile
############### 

open(IN, $gsFile);
while(<IN>){
    chomp;
    my @tmp = split("\t");
    my $id = shift(@tmp);
    shift(@tmp);
    push(@gs, $id);
    $gs2gene{$id} = [@tmp]; 
    foreach(@tmp){
        unless($gene2gs{$_}){
            $gene2gs{$_} = [$id];
        }else{
            push(@{$gene2gs{$_}}, $id);
        }
    }
}

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
   isNonsilent($type) or  next;
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    my $gene = $tmp[1];
    $removedGene{$gene} and next;
    if($rmMultiMut and $seen{$sample}{$gene}){
        next;
    }
    $seen{$sample}{$gene}=1;
    $mutCount{$gene}++;  
    $mutSum++;
    $gene{$gene}=1;
}

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

@gene = sort keys %gene;

foreach my $g (@gene){
    foreach my $gs (@{$gene2gs{$g}}){
	if($geneLength{$g}){
	    $gsLength{$gs} +=  $geneLength{$g};
	    $mutCount2{$gs} +=  $mutCount{$g};
	}
    }   
}

@recurrentMutatedGs = grep {$mutCount2{$_} and $mutCount2{$_} >= $recurrentCutoff} @gs;
@recurrentMutatedGs = sort {$mutCount2{$b} <=> $mutCount2{$a}} @recurrentMutatedGs;
map {$recurrentMutatedGs{$_}=1} @recurrentMutatedGs;

##############################
#get  poisoon parameter
##############################

map {$poissonParameter{$_} = $mutSum*$gsLength{$_}/$geneLengthSum}  @recurrentMutatedGs;

####################
#print R input file
####################

my $parameterFile = "tmp${$}.para.tsv";
my $statisticsFile = "tmp${$}.stat.tsv";

open(OUT, ">$parameterFile");
print OUT join("\n",map {$poissonParameter{$_}} @recurrentMutatedGs)."\n";
close(OUT);

open(OUT, ">$statisticsFile");
print OUT join("\n",map {$mutCount2{$_}} @recurrentMutatedGs)."\n";
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
for(my $i=0; $i<@recurrentMutatedGs; $i++){
     $pvalue{$recurrentMutatedGs[$i]}=$pvalue[$i];
}

system("rm tmp${$}.*");

####################
#print result
####################

foreach(sort {$pvalue{$b} <=> $pvalue{$a}} @recurrentMutatedGs){
    print $_."\t".$pvalue{$_}."\t".$mutCount2{$_}."/".$gsLength{$_}."\t";
    my @tmp =  grep {$mutCount{$_}} @{$gs2gene{$_}};
    @tmp = sort {$mutCount{$b} <=> $mutCount{$a}} @tmp;
    print join(",", map {$_."(".$mutCount{$_}.")"} @tmp)."\n";
}
