#!/usr/local/bin/perl 

use strict;
use warnings;
                                                                              
my %count; 
my %sum;

my $usage = "usage: $0 [-i (include indel) -c (show count)]  mutfile\n";
my $indel = 0;
my $count = 0;
my $argv = join(" ", @ARGV);
if($argv =~ s/-i//){
    $indel = 1;
}
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

my @sample = sort keys %sum;
my @context = sort keys %count; 

if(!$count){
    foreach my $s (@sample){
	my $sum = $sum{$s};
	map { $count{$_}{$s} = ($count{$_}{$s}?$count{$_}{$s}:0)/$sum } @context;
    }
}

print join("\t", @sample)."\n";
foreach my $c (@context){
    print join("\t", ($c, map {$count{$c}{$_}?$count{$c}{$_}:0} @sample))."\n";
}

