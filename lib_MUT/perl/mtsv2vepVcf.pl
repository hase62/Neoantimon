#!/usr/bin/perl 

use strict;
use warnings;

print "##fileformat=VCFv4.0\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

while(<>){
    my @tmp = split("\t");
    $tmp[1] =~ /^\d+$/ or next;
    $tmp[2] =~ /^\d+$/ or next;
    my $id =  join(":", ($tmp[0],$tmp[1],$tmp[2]));
    print join("\t", ($tmp[1], $tmp[2], $id, $tmp[3], $tmp[4], ".", ".", "."))."\n";
}
