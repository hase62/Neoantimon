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
my $status = 0;
my $db = get_home_dir()."/polyphen/polyphen-2.2.2-whess-2011_12.sqlite";

open(IN, $infile);

if($header){
    my $tmp = <IN>;     
   print $tmp;
}

while(<IN>){
    chomp;
    my ($sample, $chr, $pos, $from, $to, $context, $type) = split("\t");
    $type or next;
    my @type = split(";",$type);
    if(@type>4){
	$type = join(";",@type[0..3]);
    }
    my $in = join("\t", ($sample, $chr, $pos, $from, $to, $context, $type));
    $chr =~ s/^chr//;
    $chr =~ s/23/X/;
    $chr =~ s/24/Y/;
    $chr = "chr".$chr;
    if($type =~ /^missense/){
	my $select = "select prediction from features where chrom = \"$chr\"  and chrpos = \"$pos\" and refa = \"${from}${to}\"";
	my $command = "sqlite3 $db '$select'";
	if(system("$command > tmp${$}.txt")){
	    warn ("ERR: sqlite connection failed!");
	    $status = 1;
	    next;
	}
	chomp(my @out = `cat tmp${$}.txt`);
	`rm tmp${$}.txt`;
	my $out;
	if(grep {$_ eq "probably damaging"} @out){
	    $out = "probably damaging";
	}elsif(grep {$_ eq "possibly damaging"} @out){
	    $out = "possibly damaging";
	}elsif(grep {$_ eq "benign"} @out){
	    $out = "benign";
	}else{
	    $out = "";
	}
	$out =~ s/ //;
	print join(";", $in,$out)."\n";
    }else{
	print $in."\n";;
    }
}
exit($status);
