#!/usr/local/bin/perl                                                                                      
                                                                                                          
use strict;                                                                                               
use warnings;                                                                                             
use IO::File;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;                                                                   
set_environment();


my $annovarDir = get_home_dir()."/annovar";
my $tmpFile = "annovar.tmp${$}";

my $usage =  "usage: $0 [-h  (if header exsists)] mutfile\n";
my $argv = join(" ", @ARGV);
my $header;
if($argv =~ s/-h//){
    $header = 1;
}
@ARGV = split(" ", $argv);
my $infile = shift(@ARGV) or die $usage;

open(IN, $infile);
if($header){
    chomp(my $tmp = <IN>);
    print $tmp."\ttype\n";
}
open(OUT, ">${tmpFile}.in");

while(<IN>){
    chomp;
    my @tmp = split("\t");
    my ($s, $c, $pos, $from, $to) = @tmp[0..4];
    my $out = join("\t", ($c,$pos, $pos +length($from)-1, $from, $to));
    print OUT $out."\n";
}

#my $annovar = "perl $annovarDir/annotate_variation.pl -out ${tmpFile}.out  -buildver hg19 ${tmpFile}.in $annovarDir/humandb";

my $annovar = "perl $annovarDir/table_annovar.pl ${tmpFile}.in   $annovarDir/humandb -buildver hg19 -out ${tmpFile}.out -remove ";
$annovar .=  "-protocol refGene -operation g ";

system($annovar);

open(IN, "${tmpFile}.out.hg19_multianno.txt");

my %typemap = ( 
    "frameshift deletion"=>"frameshift",
    "frameshift insertion"=>"frameshift",
    "nonframeshift deletion"=>"missense",
    "nonframeshift insertion"=>"missense",
    "nonsynonymous SNV"=>"missense",
    "stopgain"=>"nonsense",
    "stoploss"=>"missense",
    "synonymous SNV"=>"synonymous",
    "unknown"=>""
    );


chomp(my $tmp = <IN>);
my @field = split("\t",$tmp);
my %field2index;
map {$field2index{$field[$_]} = $_} (0..$#field);

my %annovar;
while(<IN>){
    chomp;
    my @tmp =split("\t");
    my $f = $tmp[$field2index{"Func.refGene"}];
    my $ef =  $tmp[$field2index{"ExonicFunc.refGene"}];
    my $gene = $tmp[$field2index{"Gene.refGene"}];
    my $aac = $tmp[$field2index{"AAChange.refGene"}];
    $ef or $ef = "";
    $gene  or $gene = "";
    $aac or $aac = "";

    my $muttype = "";
    if($muttype =~ /splicing/){
	$muttype = "splice";
    }elsif($ef){
	$muttype = $typemap{$ef};
    }
    if(!defined($muttype) and $ef =~  /;/){
	my @tmp = split(";", $ef);
	$muttype = $typemap{$tmp[0]};
    } 
    if($gene =~ /[^\(,]+/){
	$gene  = $&;
    }

    if($aac =~ /p\.([A-Za-z]+\d+[A-Za-z]+)/){
	$aac = $1;
    }else{
	$aac = "";
    }
    my $tmp = join("\t", @tmp[0..4]);
    $annovar{$tmp} = join(";", ($muttype,$gene,$aac));
}

open(IN, $infile) or die $infile;
if($header){
    chomp(my $tmp = <IN>);
}

while(<IN>){
    chomp;
    my @tmp = split("\t");
    my ($s, $c, $pos, $from, $to) = @tmp[0..4];
    my $out = join("\t", ($c,$pos, $pos +length($from)-1, $from, $to));
    print  $_."\t".$annovar{$out}."\n";
}



`rm $tmpFile.*`;
