#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();

my %mutMatrix;
my %mutCount;
my %sample;
my $recurrentCutoff = 3;
my $damaging = 0;

my $usage = "usage: $0 [-c recurence_cutoff -d (get damaging mutation)] mutfile\n";
my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-d//){
    $damaging = 1;
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
    if($damaging?isDamaging($type):isNonsilent($type)){
	$mutMatrix{$gene}{$sample}++;
	if($mutMatrix{$gene}{$sample}==1){
	    $mutCount{$gene}++;
	}
    }
    $sample{$sample}=1;
}

my @sample = sort keys %sample;
my @gene = sort {$mutCount{$b} <=> $mutCount{$a}} grep {$mutCount{$_} >= $recurrentCutoff} keys %mutCount;
print join("\t", ("", @sample))."\n";
foreach my $gene (@gene){
	print join("\t", ($gene, map {$mutMatrix{$gene}{$_}?1:0} @sample))."\n";
}
