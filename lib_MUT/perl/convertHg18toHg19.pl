#!/usr/local/bin/perl                                                                                      
                                                                                                          
use strict;                                                                                               
use warnings;                                                                                             

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;                                                                                        
set_environment();

my $usage = "usage: $0 [-h  (if header exsists)] mutfile\n";
my $argv = join(" ", @ARGV);
my $header;
if($argv =~ s/-h//){
    $header = 1;
}
@ARGV = split(" ", $argv);
my $infile = shift(@ARGV) or die $usage;

require MUTensembl;
                                                                                                          
open(IN, $infile);
if($header){
    my $tmp = <IN>;                                                                                          
    print $tmp;
}

while(<IN>){                                                                                                
    chomp;
    my @tmp  = split("\t"); 
    my ($chr, $pos) = @tmp[1,2];
    $chr =~ s/^chr//;                                                                                     
    $chr = ($chr eq "X")?23:$chr;                                                                           
    $chr = ($chr eq "Y")?24:$chr;
    ($chr, $pos) = convert36to37($chr, $pos);
    $tmp[1] = $chr;
    $tmp[2] = $pos;
    print join("\t", @tmp)."\n";
}
