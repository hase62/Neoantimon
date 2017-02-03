#!/usr/local/bin/perl 

use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;
set_environment();


###############
#set variables
###############

my @gene;
my %gene;
my @recurrentMutatedGene;
my %recurrentMutatedGene;

my $recurrentCutoff = 5;
my $rmMultiMut = 0;
my $maxCount = 100;

my $useCosmic = 0;
my $useOnlyMissense = 0;
my %cosmic;

my %pos;
my %mutCount;
my %aaLength;
my %pvalue;
my %out;
my $status;

my $home = get_home_dir();

my $usage = "usage: $0 [-c recurence_cutoff  -m (remove_multiple_mutations_within_a_sample) -C (use_cosmic) -M max_count -n (use only missense)] mutfile\n";
my $argv = join(" ", @ARGV);
if($argv =~ s/-c\s+([\d.]+)//){
    $recurrentCutoff = $1;
}
if($argv =~ s/-m//){
    $rmMultiMut  = 1;
}
if($argv =~ s/-C//){
    $useCosmic  = 1;
}
if($argv =~ s/-n//){
    $useOnlyMissense  = 1;
}
if($argv =~ s/-M\s+(\d+)//){
    $maxCount  = $1;
}

@ARGV = split(" ", $argv);

my $mutFile = shift or die $usage;
my $geneLengthFile = "${home}/data/cdsLength.tsv";
my $RnmcFile = "${home}/R/myNMC.R";
my $cosmicFile = "${home}/data/CosmicCodingMuts_v64_02042013_noLimit.vcf";

###############                                                                                                              
#read mutFile                                                                                                                
###############                                                                                                              

open(IN,$mutFile);

my  %seen; 
while(<IN>){
    chomp;
    my @tmp =  split("\t");
    my $sample = $tmp[0];
    my $context = $tmp[5];
    my $type = $tmp[6];
    $type or next;
    @tmp  = split(";",$type);
    @tmp >= 2 or next;
    $type = $tmp[0];
    my $gene = $tmp[1];
    if($useOnlyMissense){
	$type eq "missense" or next;
    }else{
	isNonsilent($type) or  next;
    }
    $mutCount{$gene}++;
    $gene{$gene}=1;
    if($rmMultiMut and $seen{$sample}{$gene}){
        next;
    }
    $seen{$sample}{$gene}=1;
    $tmp[2] or next;
    $tmp[2] =~ /^\w(\d+)\w+$/ or next;
    my $pos = $1;
    if($pos{$gene}){
	push(@{$pos{$gene}}, $pos);
    }else{
	$pos{$gene} = [$pos];
    }
}


@gene = sort keys %gene;
@recurrentMutatedGene = grep {$mutCount{$_} and $mutCount{$_} >= $recurrentCutoff and $pos{$_}} @gene;
@recurrentMutatedGene = sort {$mutCount{$b} <=> $mutCount{$a}} @recurrentMutatedGene;
map {$recurrentMutatedGene{$_}=1} @recurrentMutatedGene;


#################                                                                                                           
#set aa length                                                                                                             
#################                                                                                                            

open(IN, $geneLengthFile);
while(<IN>){
    chomp;
    my @tmp =  split("\t");
    $aaLength{$tmp[0]} = $tmp[1]/3;
}
@recurrentMutatedGene =  grep {$aaLength{$_}} @recurrentMutatedGene;


#################
#get Cosmic 
#################

if($useCosmic){
    open(IN, $cosmicFile);
    while(<IN>){
	/^#/ and next;
	chomp;
	my @tmp = split("\t");
        @tmp = split(";", $tmp[7]);
        $tmp[0] =~ /GENE=(\S+)/  or next;
	my $gene = $1;
	$tmp[3] =~ /AA=p\.(\S)(\d+)(\S)/ or next;
	$1 eq  $3 and next;
	if($cosmic{$gene}){
	    push(@{$cosmic{$gene}}, $2);
	}else{
	    $cosmic{$gene} = [$2];
	}
    }
}


#################
#run NMC by SGE
#################

sub sample {
    my @b;
    my $a;
    if(ref($_[0])){
	@b = @{shift(@_)};
	$a = shift;
	$a or $a = @b;
    }else{
	@b   = @_[0..@_-2];
	$a = pop;
    }
    my @c;
    for(my $i = 0; $i < $a; $i++){
	push(@c, splice(@b, int(rand(@b)), 1));
    }
    return @c;
}

my %pos0 = %pos;

#if(1){
for my $gene (@recurrentMutatedGene){

    my $tmpFilePrefix = "tmp${$}.$gene";
    my $RinputFile = "$tmpFilePrefix.in";
    my $RoutputFile  = "$tmpFilePrefix.out";
    my $RscriptFile  = "$tmpFilePrefix.R";
    my $scriptFile  = "$tmpFilePrefix.pl";
    my $errFile  = "$tmpFilePrefix.err";
    my @pos = @{$pos{$gene}} or next;
    if(@pos > $maxCount){
	@pos = sample(@pos, $maxCount);
    }else{
	if($cosmic{$gene}){
	    warn "added cosmic mtations for ${gene}.\n";
 	    if(@pos +  @{$cosmic{$gene}}  <= $maxCount){
		push(@pos, @{$cosmic{$gene}});
	    }else{
		push(@pos, sample(@{$cosmic{$gene}}, $maxCount-@pos));
	    }
	}
    }
    $pos{$gene} = [@pos];

    open(OUT, ">$RinputFile");
    #die join(" ",scalar(@{$pos{$gene}}));
    for my $pos (@{$pos{$gene}}){
	my @tmp = ("0") x $aaLength{$gene};
	unless($pos <=  @tmp){
	    warn "warn: pos and length are inconsistent ($pos".'@'."$gene)\n";
	    next;
	} 
	$tmp[$pos-1] = 1;
	print OUT join("\t", @tmp)."\n";
    }
    close(OUT);
    
    open(OUT, ">$RscriptFile");
    print OUT "M <- read.table(\"$RinputFile\")\n";
    print OUT "source(\"$RnmcFile\")\n";
    print OUT "N<-nmc(M)\n";
    print OUT "if(is.null(N)){N<-1}\n";
    print OUT "if(is.matrix(N)){N<-N[1,]}\n";
    print OUT "write(N, \"$RoutputFile\", sep=\"\\n\")\n";
    close(OUT);

    print_SGE_script("R --vanilla < $RscriptFile", "${tmpFilePrefix}.pl");
    while(system("qsub  -cwd  -o /dev/null -e ${tmpFilePrefix}.err   ${tmpFilePrefix}.pl > /dev/null")!=0){
	sleep(10);
    }
}
   
wait_for_SGE_finishing("tmp${$}");
#}

#################
#get pvalue
#################

$status=0;
for my $gene (@recurrentMutatedGene){
    my $tmpFilePrefix = "tmp${$}.$gene";
    #my $tmpFilePrefix = "tmp6022.$gene";
    my $RoutputFile  = "$tmpFilePrefix.out";
    my $errFile  = "$tmpFilePrefix.err";
    if(-f $RoutputFile  and -f $errFile){
	chomp(my @err = `cat $tmpFilePrefix.err`);
	chomp(my @out = `cat $tmpFilePrefix.out`);
	if(grep {/^ended/} @err){
	    if(@out==5){
		my($length,  $start, $end, $count, $pvalue) = @out;
		if($pvalue==0){
		    $pvalue = 100;
		}elsif($pvalue==1){
		    $pvalue = 1;
		}else{
		    $pvalue = -log($pvalue)/log(10);
		}
		$pvalue{$gene} = $pvalue;
		if($useCosmic){
		    my $count0 = scalar(grep {$_ >= $start and $_ <= $end} @{$pos0{$gene}});
		    $out{$gene} = join("\t", ($gene, $pvalue, $count."/".scalar(@{$pos{$gene}})."(".$count0."/".scalar(@{$pos0{$gene}}).")",  $length."/".$aaLength{$gene}, $start."-".$end))."\n";
		}else{
		    $out{$gene} = join("\t", ($gene, $pvalue, $count."/".scalar(@{$pos{$gene}}),  $length."/".$aaLength{$gene}, $start."-".$end))."\n";
		}
	    }else{
		$pvalue{$gene} = 0;
		$out{$gene} = join("\t", ($gene, 0, "NA", "NA", "NA"))."\n";
	    }
	}else{
	    $status = 1;
	    $pvalue{$gene} = -1;
	    $out{$gene} = join("\t", ($gene, "NA", "NA", "NA", "NA"))."\n";
	}
    }else{
	$status = 1;
	$pvalue{$gene} = -1;
	$out{$gene} = join("\t", ($gene, "NA", "NA", "NA", "NA"))."\n";
    }
}

foreach my $gene (sort {$pvalue{$b} <=> $pvalue{$a}} keys %pvalue){
    print $out{$gene};
}



unless($status){
    `rm tmp${$}.*`; 
}

exit($status);
