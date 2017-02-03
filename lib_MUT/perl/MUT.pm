#!/usr/local/bin/perl 

use strict;
use warnings;

####home dir###
#my $HOME = $ENV{HOME}."/MUT";
my $HOME = "/Users/takanorihasegawa/Git/Neoantimon/lib_MUT";
sub get_home_dir{
    return $HOME;
}

sub set_environment{
    my $isN =  ($ENV{HOSTNAME}=~/^n/)?1:0;
    chomp(my @lib = `ls  $HOME/java/lib`);
    my $javapath = "/usr/local/package/java/current6/bin";
    unless(grep {$_ eq $javapath} split(":", $ENV{PATH})){
        $ENV{PATH} = $javapath.":".$ENV{PATH};
    }
    my $rpath = $isN?"/usr/local/package/r/current/bin":"/usr/local/package/r/2.15.1_gcc/bin";
    unless(grep {$_ eq $rpath} split(":", $ENV{PATH})){
        $ENV{PATH} = $rpath.":".$ENV{PATH};
    }
    $ENV{CLASSPATH} = $HOME."/java/bin:$HOME/java/lib/".join(":".$HOME."/java/lib/",@lib);   
    $ENV{PERL5LIB} = "$HOME/perl";
    $ENV{ENSEMBL_REGISTRY} =  $HOME."/data/ensembl_init";
    #$ENV{R_LIBS} = $HOME."/R";
}

sub get_filename_without_suffix{
    if($_[0] =~ /^([^.\/]+)([^\/]*)$/){
        return $1;
    }elsif($_[0] =~ /((.*?\/)*([^.]+))/){
        return $3;
    }else{
        return;
    }
}

sub get_dir_and_filename_without_suffix{
    if($_[0] =~ /^([^.\/]+)([^\/]*)$/){
        return $1;
    }elsif($_[0] =~ /((.*?\/)*([^.]+))/){
        return $&;
    }else{
        return;
    }
}

sub print_SGE_script{
    my $command = shift;
    my $fh;
    my @env = grep {defined($ENV{$_})} qw(PATH PERL5LIB R_LIBS CLASSPATH LD_LIBRARY_PATH BOWTIE_INDEXES);
  my @out = (
      '#! /usr/local/bin/perl',
      '#$ -S /usr/local/bin/perl',
      #'#$ -v '.join (",",map {$_."=".$ENV{$_}} @env),
      );
    foreach(@env){
        push(@out , '$ENV{'.$_.'}="'.$ENV{$_}.'";');
    }
  push(@out, (
           'warn "command : '. $command.'\n";',
           'warn "started @ ".scalar(localtime)."\n";',
           "if(system (\"$command\" )){",
           'die "failed @ ".scalar(localtime)."\n";',
           "}else{",
           'warn "ended @ ".scalar(localtime)."\n";',
            "}"
       ));
  
    if(@_){
        open($fh, ">$_[0]");
        print $fh join("\n",@out)."\n";
        `chmod a+x $_[0]`;
    }else{
        print  join("\n",@out)."\n";
    }
}

sub  wait_for_SGE_finishing{
    my $script = shift;
    my $cutoff;
    if(@_){
        $cutoff = shift;
    }else{
        $cutoff = 1;
    }
    $script = substr($script,0,10);
    while(1){
        while(system("qstat > /dev/null") != 0){
            sleep(10);
        }
        my $out = `qstat| grep $script | wc`;
        $out =~ /\d+/;
        if($& < $cutoff ){
            return;
        }else{
            sleep(10);
        }
    }
}



#####################################
#
#get mutation  context: 
#
#0:  [^T](C>G)N
#1:  T(C>G)N
#2:  [^T](C>A)N
#3:  T(C>A)N
#4:  [^T](C>T)[^G]
#5:  T(C>T)[^G]
#6:  [^T](C>T)G
#7:  T(C>T)G
#8:  N(A>T)N 
#9:  N(A>C)N
#10: N(A>G)N
#11: indel
#
#####################################

sub getContext0{    
    my @context = (
	'[^T](C>G)N',
	'T(C>G)N',
	'[^T](C>A)N',
	'T(C>A)N',
	'[^T](C>T)[^G]',
	'T(C>T)[^G]',
	'[^T](C>T)G',
	'T(C>T)G',
        'N(A>T)N',
	'N(A>C)N',
	'N(A>G)N',
	'indel'
	);

    my @regexp  = (
	'[^T]\(C>G\)[ATGC]',    
	'T\(C>G\)[ATGC]',       
	'[^T]\(C>A\)[ATGC]',    
	'T\(C>A\)[ATGC]',       
	'[^T]\(C>T\)[^G]',      
	'T\(C>T\)[^G]',         
	'[^T]\(C>T\)G',         
	'T\(C>T\)G',            
	'[ATGC]\(A>T\)[ATGC]',  
	'[ATGC]\(A>C\)[ATGC]',  
	'[ATGC]\(A>G\)[ATGC]',
	'[ATGC]\((->[ATGC]|[ATGC]>-)\)[ATGC]'
	);

    my $fromTriplet  = shift;
    my $to = shift;
    if($2 eq "T" or $2 eq "G"){
	$fromTriplet =~ tr/ATGC/TACG/;
	$fromTriplet = reverse($fromTriplet);
	$to =~ tr/ATGC/TACG/;
    }
    $fromTriplet  =~ /([ATGC-])([ATGC-])([ATGC-])/ ;
    my $context = "${1}(${2}>${to})${3}";
    my @index = grep {$context =~ /$regexp[$_]/} (0..$#regexp);
    
    if(@index==1){ 
	return $context[$index[0]];
    }else{
	return;
    }
}

sub getContext{
    my $fromTriplet  = shift;
    my $to = shift;
    if($2 eq "T" or $2 eq "G"){
        $fromTriplet =~ tr/ATGC/TACG/;
        $fromTriplet = reverse($fromTriplet);
        $to =~ tr/ATGC/TACG/;
    }
    $fromTriplet  =~ /([ATGC-])([ATGC-])([ATGC-])/ ;
    if($2 eq "-" or $to eq "-"){
	return "indel";
    }else{
	return "${1}(${2}>${to})${3}";
    }1
}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;

    my(%genetic_code) = (

    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
	);

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        die "Bad codon \"$codon\"!!\n";
    }
}

sub isNonsilent{
    my %nonSilentType = (
	"missense" => 1,
	"nonsense" => 1,
	"frameshift" => 1,
	"splice" => 1,
	);
    
    my $type  =  shift;
    my @tmp = split(";", $type);
    if($nonSilentType{$tmp[0]}){
	return 1;
    }else{
	return 0;
    }
}

sub isDamaging{
    my %nonSilentType = (
        "missense" => 1,
        "nonsense" => 1,
        "frameshift" => 1,
        "splice" => 1,
        );

    my $type  =  shift;
    my @tmp = split(";", $type);
    if($tmp[0] eq "missense"){
	if($type =~  "probablydamaging"){
	#if($tmp[4] and $tmp[4] eq "probablydamaging"){
	#if($tmp[5] and $tmp[5] >= 25){ #use CADD score  
	 return 1;
	}else{
	    return 0;
	}
    }elsif($nonSilentType{$tmp[0]}){
	return 1;
    }else{
	return 0;
    }
}

1;
