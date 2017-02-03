#!/usr/local/bin/perl
use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();

my $javaHeap = 2048; 

my $command =  "java  -Xms${javaHeap}m -Xmx${javaHeap}m clone.CCFP ";
exit(system($command.join(" ", @ARGV)));
