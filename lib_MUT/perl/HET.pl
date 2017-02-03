#!/usr/local/bin/perl 

use strict;
use warnings;

my %gene2neighbor;
my %pvalue;
my %egde2pvalue;
my @gene;

my $cutoff = 2;

my $usage = "usage: $0  [-c cutoff] pvaluefile ppiFile\n";
my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $cutoff = $1;
}
@ARGV = split(" ", $argv);

my $pvalueFile = shift or die $usage;
my $ppiFile = shift or die $usage;

my %seen;
open(IN, $ppiFile);
while(<IN>){
    chomp;
    my @tmp = split("\t");
    @tmp = sort @tmp[0,2];
    my $id1 = $tmp[0];
    my $id2 = $tmp[1];
    $seen{$id1}{$id2} and next;
    $seen{$id1}{$id2} = 1;
    $id1 eq $id2 and next;
    if($gene2neighbor{$id1}){
        push(@{$gene2neighbor{$id1}}, $id2);
    }else{
        $gene2neighbor{$id1} =  [$id2];
    }
    if($gene2neighbor{$id2}){
        push(@{$gene2neighbor{$id2}}, $id1);
    }else{
        $gene2neighbor{$id2} =  [$id1];
    }
}

open(IN,$pvalueFile);
while(<IN>){
    chomp;
    my @tmp = split("\t");
    $tmp[1] < $cutoff and  next;
    $pvalue{$tmp[0]} = $tmp[1];
}

@gene =  keys %pvalue;
 foreach my $g1 (grep {defined($pvalue{$_})}  @gene){
     foreach my $g2 (grep {defined($pvalue{$_})}  @{$gene2neighbor{$g1}}){
	 #my $p = $pvalue{$g2} - log(scalar(@{$gene2neighbor{$g1}}))/log(10);
	 #my $p1 = $pvalue{$g2} - log(scalar(@{$gene2neighbor{$g1}}))/log(10);
	 #my $p2 = $pvalue{$g1} - log(scalar(@{$gene2neighbor{$g2}}))/log(10);
	 my $p1 = $pvalue{$g2};
         my $p2 = $pvalue{$g1};
	 my $p = ($p1+$p2)/2;
	 $p < 0 and $p = 0;
	 ($g1, $g2) = sort(($g1, $g2));
	 my $id = $g1."[".$pvalue{$g1}."]"."\t".$g2."[".$pvalue{$g2}."]";
	 if(!$egde2pvalue{$id} or $p < $egde2pvalue{$id}){
	     $egde2pvalue{$id} = $p;
	 }
     }
} 

foreach(sort {$egde2pvalue{$b} <=> $egde2pvalue{$a}} keys %egde2pvalue){
    $egde2pvalue{$_} < $cutoff and next; 
    print  $_."\t".$egde2pvalue{$_}."\n";
}


