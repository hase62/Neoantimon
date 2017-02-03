#!/usr/local/bin/perl                                                                                                                     

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
require Cwd;
set_environment();

my $usage = "usage: $0 [-h (if header exsists)]  infile\n";
my $header;
my $argv = join(" ", @ARGV);
if($argv =~ s/-h//){
    $header = 1;
}
@ARGV = split(" ", $argv);

@ARGV  or die $usage;

my $home = get_home_dir();

my $infile = shift;
    
system("perl ${home}/perl/processParallelly.pl -n 10 ".($header?"-h":"")." ${home}/perl/addContextAndType.pl $infile  > tmp${$}.1.tsv");
system("perl ${home}/perl/processParallelly.pl -n 10 ".($header?"-h":"")." ${home}/perl/addPolyPhen.pl tmp${$}.1.tsv  > tmp${$}.2.tsv");
print `cat tmp${$}.2.tsv`;
`rm tmp${$}.*`;


