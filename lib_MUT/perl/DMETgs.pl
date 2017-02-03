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
my $nonsilent = 0;
my %bgGene;

my @gs;
my %gene2gs;
my %gs2gene;

my %removedGene;

my %trial;        #  $trial{$gene}
my %trial2;       #  $trial2{$gs}  
my $binomialParameter;
my %success;      # $success{$gene}
my %success2;     # $success2{$gene}
my %pvalue; #$pvalue{$gs}

my $home = get_home_dir();

my $usage = qq(usage: $0 [options] mutfile gsfile
       -c recurence_cutoff
       -r renoved_gene (gene1:gene2...)
       -m (remove_multiple_mutations_within_a_sample)
       -n (count_nonsilent_mutaions)
       -b background_gene (gene1:gene2....) 
);
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
if($argv =~ s/-n//){
    $nonsilent  = 1;
}
if($argv =~ s/-b\s+(\S+)//){
    map {$bgGene{$_}=1} split(":",$1);
}
@ARGV = split(" ", $argv);

my $mutFile = shift or die $usage;
my $gsFile = shift or die $usage;


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
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    my $gene = $tmp[1];
    $removedGene{$gene} and next;
    if($rmMultiMut and $seen{$sample}{$gene}){
        next;                                                                                                                 
    }
    $seen{$sample}{$gene}=1;
    $gene{$gene}=1;
    $trial{$gene}++;
    if($nonsilent?isNonsilent($type):isDamaging($type)){
	$success{$gene}++;
    }
}


@gene = sort keys %gene;
foreach my $g (@gene){
    foreach my $gs (@{$gene2gs{$g}}){
	$trial2{$gs} +=  $trial{$g};
	if($success{$g}){
	    $success2{$gs} += $success{$g};
	}
    }
}

@recurrentMutatedGs = grep {$trial2{$_} and $trial2{$_} >= $recurrentCutoff} @gs;
@recurrentMutatedGs = sort {$trial2{$b} <=> $trial2{$a}} @recurrentMutatedGs;
map {$recurrentMutatedGs{$_}=1} @recurrentMutatedGs;
map {$success2{$_}?():($success2{$_}=0)} @recurrentMutatedGs;
map {$success{$_}?():($success{$_}=0)} @gene;
map {$trial{$_}?():($trial{$_}=0)} @gene;
    
################################
#get binomial parameter1
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
my $successFile = "tmp${$}.success.tsv";
my $geneFile = "tmp${$}.gene.txt";


open(OUT, ">$trialFile");
print OUT join("\n",map {$trial2{$_}} @recurrentMutatedGs)."\n";
close(OUT);

open(OUT, ">$successFile");
print OUT join("\n",map {$success2{$_}} @recurrentMutatedGs)."\n";
close(OUT);

open(OUT, ">$geneFile");
print OUT join("\n",@recurrentMutatedGs)."\n";
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

for(my $i=0; $i<@recurrentMutatedGs; $i++){
    $pvalue{$recurrentMutatedGs[$i]}=$pvalue[$i];
}

system("rm tmp${$}.*");

####################
#print result
####################

foreach(sort {$pvalue{$b} <=> $pvalue{$a}} @recurrentMutatedGs){
    print join("\t", ($_, $pvalue{$_}, $success2{$_}."/".$trial2{$_}))."\t";
    my @tmp =  grep {$success{$_}} @{$gs2gene{$_}};
    @tmp = sort {$success{$b} <=> $success{$a}} @tmp;
    print join(",", map {$_."(".$success{$_}.")"} @tmp)."\n";
}

