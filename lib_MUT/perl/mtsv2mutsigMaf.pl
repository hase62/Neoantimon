use strict;
use warnings;



print join("\t", qw(patient	gene	effect	categ))."\n";
while(<>){
	chomp;
	my @tmp = split("\t");
	my @tmp2 = split(";", $tmp[6]);
	my $smp = $tmp[0];
	my $gene = $tmp2[1] or next;

	

	my $type1 = "nonsilent";
	if($tmp2[0] eq "synonymous"){
		$type1 = "silent";
	}



	my $type2 = $tmp[5];
	
	if($tmp[5] =~ /\(C>T\)G/){ #1 CpG transitions
		$type2 = 1
	}elsif($tmp[5] =~ /\(C>[AG]\)G/){ #2 CpG transversions
		$type2 = 2
	}elsif($tmp[5] =~ /\(C>T\)/){ #3 C:G transitions
		$type2 = 3
	}elsif($tmp[5] =~ /\(C>[AG]\)/){ #4 C:G transversions
		$type2 = 4
	}elsif($tmp[5] =~ /\(A>G\)/){ #5 A:T transitions
		$type2 = 5
	}elsif($tmp[5] =~ /\(A>[CT]\)/){ #6 A:T transversions
		$type2 = 6
	}elsif($tmp[5] eq "indel"){#7 null+indel mutations
		$type2 = 7
	}else{ 
		warn $tmp[5];
	}

	print join("\t", ($smp,$gene,$type1,$type2))."\n";	

}
