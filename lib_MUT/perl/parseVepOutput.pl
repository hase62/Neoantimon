#!/usr/bin/perl

use strict;
use warnings;

while(<>){
    if(/^#/){
       next;
   }
    last;
}
<>;



my %aaChange = (
		ala => "A",
		asx => "B",
		cys => "C",
		asp => "D",
		glu => "E",
		phe => "F",
		gly => "G",
		his => "H",
		ile => "I",
		lys => "K",
		leu => "L",
		met => "M",
		asn => "N",
		pro => "P",
		gln => "Q",
		arg => "R",
		ser => "S",
		thr => "T",
		sec => "U",
		val => "V",
		trp => "W",
		xaa => "X",
		tyr => "Y",
		glx => "Z"
	    );

my @id;
my %data;
my %tier;
my %gene;
while(<>){
    chomp;
    my @tmp = split("\t");
    my $id = $tmp[0];
    my $type = $tmp[6];
    
    #synonymous, nonsense, frameshift, missense, UTR, splice, intron, intergenic, undef
    if($type =~ /frameshift/){
	$type = "frameshift";
    }elsif($type =~ /splice/){
	$type = "splice"; 
    }elsif($type =~ /stop_gained/){
	$type = "nonsense";
    }elsif($type =~ /missense/ or $type =~ /stop_lost/ or $type =~ /initiator_codon_variant/){
	$type ="missense";
    }elsif($type =~ /UTR/){
	$type ="UTR";
    }elsif($type =~ /intron/){
	$type ="intron";
    }elsif($type =~ /intergenic/ or $type =~ /stream/){
	$type ="intergenic";
    }elsif($type =~ /synonymous/){
        $type ="synonymous";
    }else{
	$type = "undef";
    }

    my $annot = $tmp[13];
    my $symbol = "";
    my $polyphen = "",
    my $CADD = "";
    if($annot =~ /SYMBOL=([^\;]+)/){
	if(grep { $type eq $_ } ("synonymous", "nonsense", "frameshift", "missense", "splice", "UTR")){
	    $symbol = $1;
	}
    }
    if($annot =~ /PolyPhen=(\w+)/){
	$polyphen = $1;
	$polyphen =~ s/_//;
    }
    if($annot =~ /CADD_WG_SNV=([^\;]+)/){
	$CADD = $1;
	my @cadd = split(",",  $CADD);
	$CADD = $cadd[$#cadd];
    }
    my $aapos = $tmp[9];
    if($aapos =~ /\d+/){
	$aapos = $&;
    }
    my $aa = $tmp[10];
    #print join("\t", ($id, $type, $symbol, $polyphen, $CADD, $aa, $aapos))."\n";
    #print $annot."\n";
    #next;
    if(!$symbol){
	$aa = "";
    }elsif($aa =~ /\//){
	my  @tmp2 = split("/", $aa);
	$aa = $tmp2[0].$aapos.$tmp2[1];
	$aa =~ s/\*/_/;
    }elsif($type eq "frameshift" and $annot =~/\.([A-Za-z]+)$aapos/ and $aaChange{lc($1)}){
	$aa = $aaChange{lc($1)}.$aapos."fs";
    }else{
	$aa = "";
    }
    
    my  $out =  join("\t", ($id, $type, $symbol, $aa, $polyphen, $CADD))."\n";
    
    #tier 
    # 6. cannonical & damaging
    # 5. cannonical & nonsyn
    # 4. damaging
    # 3. nonsyn
    # 2. cannonical
    # 1. other

    my  $nonsyn = 0;
    if($symbol and grep { $type eq $_ } ("nonsense", "frameshift", "missense", "splice")){
	$nonsyn = 1;
    }
    my $damaging = 0;
    if($symbol and grep { $type eq $_ } ("nonsense", "frameshift", "splice")){
        $damaging  = 1;
    }

    my $t = 0;
    if(($annot =~ /CANONICAL=YES/) and $damaging){
	$t = 6;
    }elsif(($annot =~ /CANONICAL=YES/) and $nonsyn){
	$t = 5;
    }elsif($nonsyn){
	$t = 4;
    }elsif($nonsyn){
	$t = 3;
    }elsif($annot =~ /CANONICAL=YES/){
	$t = 2;
    }else{
	$t = 1;
    }

    if(!$data{$id}){
	$data{$id} = $out;
	$tier{$id} = $t;
	$gene{$id} = $symbol;
	push(@id, $id);
    }elsif($t > $tier{$id}){
	$data{$id} = $out;
	$tier{$id} = $t;
	$gene{$id} = $symbol;
    #}    
    }elsif($t == $tier{$id}  and $data{$id} ne $out){
    	#warn $t."\n";
	#warn $data{$id};
    	#warn $out."\n";
	if(length($symbol)< length($gene{$id})){
	    $data{$id} = $out;
	    $tier{$id} = $t;
	    $gene{$id} = $symbol;
	}
    }
}

for(@id){
    #print join("\t", ($_, $gene{$_}, $tier{$_}))."\n";
    print $data{$_};
}
