#!/usr/local/bin/perl

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;

my $infile = shift or die;
my $suffix = shift;
$suffix  or  $suffix = get_filename_without_suffix($infile);

my %nodePvalue;
my %edgePvalue;

open(IN, $infile);
open(OUT, ">${suffix}.sif");
while(<IN>){
    /^(\S+)\[(\S+?)\]\t(\S+)\[(\S+?)\]\t(\S+)/;
    $nodePvalue{$1} = $2;
    $nodePvalue{$3} = $4;
    $edgePvalue{"$1\tpp\t$3"} = $5;
    print OUT "$1\tpp\t$3\n";
}

open(OUT, ">${suffix}.noa.tsv");
print OUT "node\tpvalue\n";
foreach(keys %nodePvalue){
    print OUT $_."\t".$nodePvalue{$_}."\n";
}

open(OUT, ">${suffix}.eda.tsv");
print OUT "node1\tinteraction\tnode2\tpvalue\n";
foreach(keys %edgePvalue){
    print OUT $_."\t".$edgePvalue{$_}."\n";
}

close(OUT)
