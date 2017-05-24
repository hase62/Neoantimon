#!/usr/bin/perl

use strict;
use warnings;

my $usage = "usage: $0  result_file1  result_file2 ...\n";
@ARGV or die $usage;
my %data;
foreach(@ARGV){
    open(IN, $_);
    while(<IN>){
	chomp;
	my @tmp = split("\t");
	if(!defined($data{$tmp[0]})){
	    $data{$tmp[0]} = [$tmp[1]];
	}else{
	    push(@{$data{$tmp[0]}}, $tmp[1]);
	} 
    }
}

my @field = ("ID");
for(my $i = 1; $i <= @ARGV; $i++){
    push(@field, "VALUE".$i);
}
print join("\t", @field)."\n";
foreach my $k (keys %data){
    if(@{$data{$k}}==@ARGV){ 
	print join("\t", ($k, @{$data{$k}}))."\n";
    }
}
