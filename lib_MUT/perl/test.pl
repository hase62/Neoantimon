#!/usr/local/bin/perl                                                                                      
                                                                                                          
use strict;                                                                                               
use warnings;                                                                                             

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;                                                                   
set_environment();
          
unshift(@INC, get_home_dir()."/perl/ensembl/modules");
require Bio::EnsEMBL::Registry;                                                                               

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
    chomp(my $tmp = <IN>);                                                                                          
    print $tmp."\tcontext\ttype\n";
}

while(<IN>){                                                                                                
    chomp;
    my ($sample,  $chr, $mut_pos,  $from, $to, undef, $type)  = split("\t");                                                
    #$from =~ /^[ATGC-]$/ or next;                                                                          
    #$to =~ /^[ATGC-]$/ or next;                                                                            
    $chr =~ s/^chr//;                                                                                     
    #$chr = ($chr eq "X")?23:$chr;                                                                           
    #$chr = ($chr eq "Y")?24:$chr;
    my $context = getMutContext($chr, $mut_pos, $from, $to) or next;
    #my $type = getMutType($chr, $mut_pos, $from, $to) or next;
    print join("\t", ($sample,  $chr, $mut_pos, $from, $to, $context, $type))."\n";
}


