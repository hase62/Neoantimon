#!/usr/local/bin/perl

use strict;
use warnings;

unshift(@INC, get_home_dir()."/perl/ensembl/modules");
unshift(@INC, get_home_dir()."/perl/bioperl-live");
require Bio::EnsEMBL::Registry;
require Bio::EnsEMBL::Feature;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all();
my $slice_adaptor = $registry->get_adaptor('Human' ,'Core','Slice') or die;
my $gene_adaptor = $registry->get_adaptor( "Human", "Core", "Gene") or die;
my $transcript_adaptor = $registry->get_adaptor( "Human", "Core", "Transcript") or die;


sub getMutContext{
    my ($chr, $mut_pos, $from, $to) = @_;
    if($from =~  /-/ or $to =~ /-/){
        return "indel";
    }
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $mut_pos-1, $mut_pos+1) or return;
    my $seq = $slice->seq;
    $seq =~ /(\w)(\w)(\w)/ or die $seq;
    unless($2 eq $from){
        #warn "ref bases are inconsistentr!\n";
        $seq = ${1}.${from}.${3};
    }
    return getContext($seq, $to);
}

sub convert36to37{
    my ($chr, $pos) = @_[0,1];
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $pos, $pos, '1', 'NCBI36') or return;
    my $feature = new Bio::EnsEMBL::Feature(
	-start  => 1,
	-end    => 1,
	-strand => 1,
	-slice  => $slice
	);
    $feature = $feature->transform('chromosome', 'GRCh37') or return;
    $chr = $feature->slice->seq_region_name();
    $pos =$feature->start();
    return ($chr, $pos);
}



sub getMutType{
    my ($chr, $mut_pos, $from, $to) = @_;
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $mut_pos, $mut_pos) or return;
    my @genes = @{$gene_adaptor->fetch_all_by_Slice($slice)};
    unless($from ne "-" and $from eq  $slice->seq()){                                                                                
	#warn "ref bases are inconsistentr!\n";                                                                                      
    }  
    @genes = map {$_->transform('toplevel')} @genes;
    if(@genes >=1){ 
	my $gene;
	my $transcript;
	my @exons;
	for(my $i = 0; $i < @genes; $i++){
	    my $tmp = $genes[$i]->canonical_transcript();
	    my @tmp =  @{$tmp->get_all_Exons()};
	    if(grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @tmp){
		$gene = $genes[$i];
		$transcript = $tmp;
		last;
	    }	
	    my @transcripts =  @{$genes[$i]->get_all_Transcripts()} or next;
	    for(my $j = 0; $j < @transcripts; $j++){
		my @tmp =  @{$transcripts[$j]->get_all_Exons()}; 
		if(grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @tmp){
		    if(!$transcript or $transcripts[$j]->length() > $transcript->length()){ 
			$gene = $genes[$i];
			$transcript = $transcripts[$j]; 
		    }
		}
	    }
	}
	
	if($transcript){
	    @exons =   grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @{$transcript->get_all_Exons()};
	}
   

	if(@exons){
	    my $gene_id = $gene->external_name();
	    unless($transcript->is_canonical()){
		$gene_id .= "_".$transcript->stable_id();
	    }
	    #check CDS
	    #warn $gene->stable_id()." ".$gene->strand()."\n";
	    my $exon;
	    my $exon_start_in_cdna; 
	    foreach(@exons){
		$exon = shift(@exons); 
		$exon_start_in_cdna  = $exon->cdna_coding_start($transcript) and last; #in cDNA coordinate
	    }
	    unless($exon_start_in_cdna){
		return "UTR;$gene_id";
	    }
	    my $mut_pos_in_cdna;
	    if($gene->strand() == 1){
		    my $start = $exon->start();   #in genome coordinate 
		    $mut_pos_in_cdna = $exon_start_in_cdna  + $mut_pos - $start;
	    }else{
		my $end = $exon->end();     #in genome coordinate 
		$mut_pos_in_cdna = $exon_start_in_cdna  + $end - $mut_pos; 
		$from =~ tr/ATGC/TACG/;
		$to =~ tr/ATGC/TACG/; 
	    }
	    my $transcript_seq = $transcript->spliced_seq();
	    my $cds_start_in_cdna = $transcript->cdna_coding_start();
	    my $cds_end_in_cdna = $transcript->cdna_coding_end();
	    
	    if( $mut_pos_in_cdna >=  $cds_start_in_cdna  and  $mut_pos_in_cdna  <=  $cds_end_in_cdna-3  ){
		#CDS mutation
	        
		my $mut_pos_in_cds = $mut_pos_in_cdna - $cds_start_in_cdna + 1;
		my $mut_codon_pos;
		my $mut_aa_ps; 
		if(($mut_pos_in_cds+2)/3 == int(($mut_pos_in_cds+2)/3)){
		    $mut_codon_pos = 1;
		    $mut_aa_ps = ($mut_pos_in_cds+2)/3;
		}elsif(($mut_pos_in_cds+1)/3 == int(($mut_pos_in_cds+1)/3)){
		    $mut_codon_pos = 2;
		    $mut_aa_ps = ($mut_pos_in_cds+1)/3; 
		}else{
		    $mut_codon_pos = 3;
		    $mut_aa_ps = ($mut_pos_in_cds)/3; 
		}
		
		my $cds_length = $cds_end_in_cdna - $cds_start_in_cdna + 1;
		my $coding_seq = substr($transcript_seq, $cds_start_in_cdna-1 , $cds_length);
		my $mut_codon_seq_from;
		my $mut_codon_seq_to;
		
		#warn join("\t", ($mut_pos_in_cdna,  $cds_start_in_cdna, $cds_end_in_cdna, $cds_length, length($transcript_seq)))."\n";
	
		my $aa_from;
		my $aa_to;
		if($from eq "-" or $to eq "-"){
		    if($mut_codon_pos == 1){
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-1 , 3);
		    }elsif($mut_codon_pos == 2){
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-2 , 3);
		    }else{
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-3 , 3);
		    }
		    $aa_from = codon2aa($mut_codon_seq_from);
		    $aa_to = "fs";
                }else{
		    if($mut_codon_pos == 1){
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-1 , 3);
			$mut_codon_seq_from =~ /(\w)(\w)(\w)/;
			$mut_codon_seq_from = $from.$2.$3;
			$mut_codon_seq_to = $to.$2.$3;
		    }elsif($mut_codon_pos == 2){ 
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-2 , 3);
			$mut_codon_seq_from =~ /(\w)(\w)(\w)/;
			$mut_codon_seq_from = $1.$from.$3;
			$mut_codon_seq_to = $1.$to.$3;
		    }else{
			$mut_codon_seq_from = substr($coding_seq, $mut_pos_in_cds-3 , 3);
			$mut_codon_seq_from =~ /(\w)(\w)(\w)/;
			$mut_codon_seq_from = $1.$2.$from;
			$mut_codon_seq_to = $1.$2.$to;
		    }
		    
		    $aa_from = codon2aa($mut_codon_seq_from);
		    $aa_to = codon2aa($mut_codon_seq_to);
		}
		
		my $translation = $transcript->translation();
		my $pfeatures = $translation->get_all_ProteinFeatures();
		my %dom;
		while ( my $pfeature = shift @{$pfeatures} ) {
		    my $logic_name = $pfeature->analysis()->logic_name();
		    my $start  = $pfeature->start(); 
		    my $end = $pfeature->end(); 
		    my $ac = $pfeature->interpro_ac() or next;
		    my $desc = $pfeature->idesc() or next;
		    if( $mut_aa_ps >= $start  and $mut_aa_ps <= $end){
			$dom{"$desc($ac)"}=1;
		    }
		}
		
		my $ret;
		if($aa_from eq $aa_to){
		    $ret  =  "synonymous";
		}elsif($aa_to eq "_"){
		    $ret = "nonsense";
		}elsif($aa_to eq "fs"){
		    $ret = "frameshift";
		}else{
		    $ret = "missense";
		}
		$ret .= ";$gene_id;${aa_from}${mut_aa_ps}${aa_to};".join(",", sort keys %dom);
		return $ret;
	    }else{
		#UTR mutation
		return  "UTR;$gene_id";
	    }
	    
	}else{
	    #check splice mutation based on GT-AG rule 	    
	    my $gene;
	    my $transcript;
	    my @introns;
	    for(my $i = 0; $i < @genes; $i++){ 
		my @transcripts =   @{$genes[$i]->get_all_Transcripts()} or next;
		for(my $j = 0; $j < @transcripts; $j++){
		    my @tmp =  @{$transcripts[$j]->get_all_Introns()};
		    if(grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @tmp){
			if(!$transcript or $transcripts[$j]->length() > $transcript->length()){ 
			    $gene = $genes[$i];
			    $transcript = $transcripts[$j]; 
			}
		    }
		}
	    }
	    unless($transcript){
		return "undef";
	    }
	    my $gene_id = $gene->external_name();
	    unless($transcript->is_canonical()){
                $gene_id .= "_".$transcript->stable_id();
            }
	    @introns =  grep {$_->start <= $mut_pos  and $_->end >= $mut_pos}  @{$transcript->get_all_Introns()}; 
	    my $intron = $introns[0];
	    my $mut = 0;
	    if($intron->strand == 1){
		if($mut_pos ==  $intron->start){
		    if($to ne "G"){
			$mut = 1;
		    }
		}elsif($mut_pos ==  ($intron->start)+1){
		    if($to ne "T"){
			$mut = 1;
		    }
		}elsif($mut_pos ==  ($intron->end)-1){
		    if($to ne "A"){
			$mut = 1;
		    }
		}elsif($mut_pos ==  $intron->end){ 
		    if($to ne "G"){
			$mut = 1;
		    }
		}
	    }else{
		if($mut_pos ==  $intron->start){
		    if($to ne "C"){
                        $mut = 1;
                    }
		}elsif($mut_pos ==  ($intron->start)+1){  
		    if($to ne "T"){
                        $mut = 1;
                    }
		}elsif($mut_pos ==  ($intron->end)-1){ 
		    if($to ne "A"){
                        $mut = 1;
                    }
		}elsif($mut_pos ==  $intron->end){
		    if($to ne "C"){
                        $mut = 1;
                    }
		}
	    }
	    if($mut){
		return "splice;$gene_id";
	    }else{
		return "intron;$gene_id";
	    }
	}
    }else{
	    return  "intergenic";
    }
}


sub getDomain{
    my ($chr, $mut_pos, $from, $to) = @_;
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $mut_pos, $mut_pos) or return;
    my @genes = @{$gene_adaptor->fetch_all_by_Slice($slice)};
    unless($from ne "-" and $from eq  $slice->seq()){
        #warn "ref bases are inconsistentr!\n";
    }
    @genes = map {$_->transform('toplevel')} @genes;
    unless(@genes >=1){
        return;
    }
    my $gene;
    my $transcript;
    my @exons;
    for(my $i = 0; $i < @genes; $i++){
        my $tmp = $genes[$i]->canonical_transcript();
        my @tmp =  @{$tmp->get_all_Exons()};
        if(grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @tmp){
            $gene = $genes[$i];
            $transcript = $tmp;
            last;
        }
        my @transcripts =  @{$genes[$i]->get_all_Transcripts()} or next;
        for(my $j = 0; $j < @transcripts; $j++){
            my @tmp =  @{$transcripts[$j]->get_all_Exons()};
            if(grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @tmp){
                if(!$transcript or $transcripts[$j]->length() > $transcript->length()){
                    $gene = $genes[$i];
                    $transcript = $transcripts[$j];
                }
            }
        }
    }
    
    if($transcript){
        @exons =   grep {$_->start <= $mut_pos and $_->end >= $mut_pos}  @{$transcript->get_all_Exons()};
    }
    
    unless(@exons){
        return;
    }
    my $gene_id = $gene->external_name();
    unless($transcript->is_canonical()){
        $gene_id .= "_".$transcript->stable_id();
    }
    #check CDS
    #warn $gene->stable_id()." ".$gene->strand()."\n";
    my $exon;
    my $exon_start_in_cdna;
    foreach(@exons){
        $exon = shift(@exons);
	if($exon->cdna_coding_start($transcript)){
	    $exon_start_in_cdna  = $exon->cdna_start($transcript) and last; #in cDNA coordinate
	}
    }
    unless($exon_start_in_cdna){
        return;
    }
    my $mut_pos_in_cdna;
    if($gene->strand() == 1){
        my $start = $exon->start();   #in genome coordinate
        $mut_pos_in_cdna = $exon_start_in_cdna  + $mut_pos - $start;
    }else{
        my $end = $exon->end();     #in genome coordinate
        $mut_pos_in_cdna = $exon_start_in_cdna  + $end - $mut_pos;
        $from =~ tr/ATGC/TACG/;
        $to =~ tr/ATGC/TACG/;
    }
    my $transcript_seq = $transcript->spliced_seq();
    my $cds_start_in_cdna = $transcript->cdna_coding_start();
    my $cds_end_in_cdna = $transcript->cdna_coding_end();
    
    unless( $mut_pos_in_cdna >=  $cds_start_in_cdna  and  $mut_pos_in_cdna  <=  $cds_end_in_cdna-3  ){
        return;
    }
    
    my $mut_pos_in_cds = $mut_pos_in_cdna - $cds_start_in_cdna + 1;
    my $mut_codon_pos;
    my $mut_aa_ps;
    if(($mut_pos_in_cds+2)/3 == int(($mut_pos_in_cds+2)/3)){
        $mut_codon_pos = 1;
        $mut_aa_ps = ($mut_pos_in_cds+2)/3;
    }elsif(($mut_pos_in_cds+1)/3 == int(($mut_pos_in_cds+1)/3)){
        $mut_codon_pos = 2;
        $mut_aa_ps = ($mut_pos_in_cds+1)/3;
    }else{
        $mut_codon_pos = 3;
        $mut_aa_ps = ($mut_pos_in_cds)/3;
    }
    
    my $translation = $transcript->translation();
    my $pfeatures = $translation->get_all_ProteinFeatures();
    my %dom;
    while ( my $pfeature = shift @{$pfeatures} ) {
        my $logic_name = $pfeature->analysis()->logic_name();
        my $start  = $pfeature->start();
        my $end = $pfeature->end(); 
        my $ac = $pfeature->interpro_ac() or next;
        my $desc = $pfeature->idesc() or next;
        if( $mut_aa_ps >= $start  and $mut_aa_ps <= $end){
            $dom{"$desc($ac)"}=1;
        }
    }   
    my $ret = join(",", sort keys %dom);
    return $ret;
}


