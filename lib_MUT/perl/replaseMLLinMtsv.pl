#!/use/bin/perl 

use strict;
use warnings;


my %map = (
	"KMT2A"=>"MLL",	
	"KMT2D"=>"MLL2",
	"KMT2B"=>"MLL4",
	"KMT2E"=>"MLL5",
	"KMT2C"=>"MLL3"
	);	

while(<>){
	chomp;	
	my @tmp = split("\t");	
	my @tmp2 = split(";", $tmp[6]);

	if($map{$tmp2[1]}){
		$tmp2[1] = $map{$tmp2[1]};
	}
	$tmp[6] = join(";",@tmp2);
	print join("\t", (@tmp))."\n";
}