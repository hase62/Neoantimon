#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();

my %mutCount;
my %seen;
my $recurrentCutoff = 3;
my $damaging = 0;
my $countSample = 0;

my $usage = "usage: $0 [-s (count_sample) -c recurence_cutoff -d (get damaging mutation)] mutfile\n";
my $argv = join(" ", @ARGV);

if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-d//){
    $damaging = 1;
}
if($argv =~ s/-s//){
    $countSample = 1;
}

@ARGV = split(" ", $argv);
my $mutFile = shift or die $usage;

open(IN,$mutFile);

while(<IN>){
    chomp;    
    my @tmp =  split("\t");
    my $sample = $tmp[0];
    my $context = $tmp[5];
    my $type = $tmp[6];
    $type or next;
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    my $gene = $tmp[1];
    if($countSample  and $seen{$gene}{$sample}){
	next;
    }
    $seen{$gene}{$sample}=1;
    if($damaging?isDamaging($type):isNonsilent($type)){
	$mutCount{$gene}++;
    }
}

my @gene = sort {$mutCount{$b} <=> $mutCount{$a}} grep {$mutCount{$_} >= $recurrentCutoff} keys %mutCount;
foreach my $gene (@gene){
    print join("\t", ($gene, $mutCount{$gene}))."\n";
}
