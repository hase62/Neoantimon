#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;                                                                                        
set_environment();

my $usage = "usage: $0 infile\n";

my $argv = join(" ", @ARGV);
my $header;
if($argv =~ s/-h//){
    $header = 1;
}
@ARGV = split(" ", $argv);


chomp(my $tmp = <>);
my @field = map {lc($_) }split("\t", $tmp);
my %index;
for(my $i = 0; $i < @field; $i++){
	    $index{$field[$i]} = $i;
}

my @tmp = grep {/barcode/} @field or die;
my $barcodeField =  $tmp[0];
@tmp = grep {/chrom/} @field or die;
my $chromosomeField = $tmp[0];
@tmp = grep {/end/} @field or die;
my $endField = $tmp[0];

#print join("\t", qw(sample chr pos normal_base tumor_base context type))."\n";
while(<>){
    chomp;
    my @tmp = split("\t");
    my $build = $tmp[$index{"ncbi_build"}];
    my ($id, $chr, $pos) = map { $tmp[$index{$_}] } ($barcodeField, $chromosomeField, $endField);
    $id =~ /(\w{4}-\w{2}-\w{4})-01/ or next;
    $id = $1;
    $chr =~ s/^chr//;
    #$chr = ($chr eq "X")?23:$chr;
    #$chr = ($chr eq "Y")?24:$chr;
    if($build == 36){
	($chr, $pos) = convert36to37($chr, $pos) or next; 
    }
    my $from =  $tmp[$index{"reference_allele"}];
    my $to = $tmp[$index{"tumor_seq_allele1"}];
    if($from eq $to){
	$to = $tmp[$index{"tumor_seq_allele2"}];
    }
    if($from eq $to){
	next;
    }
    unless($from =~  /^[\w-]$/ and $to =~  /^[\w-]$/){
	next;
    }
    print join("\t",   ($id, $chr, $pos, $from, $to))."\n";
}
