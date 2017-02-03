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
    $edgePvalue{"$1 (pp) $3"} = $5;
    print OUT "$1\tpp\t$3\n";
}

open(OUT, ">${suffix}.noa");
print OUT "nodePvalue\n";
foreach(keys %nodePvalue){
    print OUT $_." = ".$nodePvalue{$_}."\n";
}

open(OUT, ">${suffix}.eda");
print OUT "edgePvalue\n";
foreach(keys %edgePvalue){
    print OUT $_." = ".$edgePvalue{$_}."\n";
}

close(OUT)
