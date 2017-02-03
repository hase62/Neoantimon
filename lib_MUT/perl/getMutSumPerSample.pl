#!/usr/local/bin/perl 

use strict;
use warnings;
                                                                              
my %count; 
my %sum;

my $usage = "usage: $0 [-i (include indel) ]  mutfile\n";
my $indel = 0;
my $count = 0;
my $argv = join(" ", @ARGV);
if($argv =~ s/-c//){
    $count = 1;
}
@ARGV = split(" ", $argv) or die $usage;

while(<>){
    chomp;    
    my @tmp =  split("\t");
    $tmp[2] =~ /^\d+$/ or next;
    my $sample = $tmp[0];
    my $context = $tmp[5];
    if(!$indel){
	$context eq "indel" and  next;
    }
    $count{$context}{$sample}++;
    $sum{$sample}++;
}

my @sample = sort {$sum{$b} <=> $sum{$a}} keys %sum;

foreach my $s (@sample){
	print $s."\t".$sum{$s}."\n";
}

