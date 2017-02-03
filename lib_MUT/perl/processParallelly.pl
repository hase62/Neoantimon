#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
require Cwd;
set_environment();

my $usage = "usage: $0 [-n number_of_child -h (if header exsists)] child_script infile\n";
my $n = 100;

my $argv = join(" ", @ARGV);
if($argv =~ s/-n\s+([\d.]+)//){
    $n = $1;
}
my $header;
my $header2;
if($argv =~ s/-h//){
    $header = 1;
}
@ARGV = split(" ", $argv);

@ARGV==2 or die $usage;

my $child = Cwd::abs_path(shift @ARGV);
my $inFile = shift @ARGV;

open(IN, $inFile);

if($header){
    $header2 = <IN>;     
}

my $tmp = `wc $inFile`;
$tmp =~ /\d+/;
$tmp = $&;
$n = int($tmp/$n);

my $i = 0;
my $j = 1;
open(OUT, ">tmp${$}.$j.in");
if($header){
    print OUT $header2;
}
foreach(<IN>){
    if($i==$n){
	$j++;
	open(OUT, ">tmp${$}.$j.in");
	if($header){
	    print OUT $header2;
	}
	$i=0;
    }
    print OUT;
    $i++;
}
close(OUT);

my %flag;
my $mem = 2;
L: while(1){ 
    for(my $i = 1; $i <= $j; $i++){
	$flag{$i} and  next;
	my $script = "$child ".($header?"-h":"")." tmp${$}.$i.in";
	print_SGE_script("perl $script", "tmp${$}.$i.pl");
	while(system("qsub -N tmp${$}.$i  -l s_vmem=${mem}G -l mem_req=${mem} -cwd  -o tmp${$}.$i.out -e tmp${$}.$i.err tmp${$}.$i.pl > /dev/null")){
	    sleep(10);
	}
    }
    wait_for_SGE_finishing("tmp${$}");
    
    for(my $i = 1; $i <=  $j; $i++){
	my $tmp = `tail -1  tmp${$}.$i.err`;
	if($tmp =~ /^ended/){
	    $flag{$i} =1 
	}
    }
    my @tmp = grep {!$flag{$_}} (1..$j);
    foreach (@tmp){
	warn "the ${_}-th job failed!\n";
    }
    if(@tmp){
	if($mem == 32){
	    die "not enough memory!";
	}
	$mem *= 2;
	next;
    }else{
	last;
    }
}

for(my $i = 1; $i <=  $j; $i++){
    open(IN, "tmp${$}.$i.out");
    if($header){
	$header2 = <IN>;
	if($i==1){
	    print $header2;
	}
    }
    while(<IN>){
	print;
    }
}

`rm tmp${$}.*`;
