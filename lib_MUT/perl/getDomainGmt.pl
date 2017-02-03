#!/usr/local/bin/perl                                                                                                                            
                                                                                                                                                 
use strict;                                                                                                                                      
use warnings;                                                                                                                                    
$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require MUT;                                                                   
set_environment();
          
unshift(@INC, get_home_dir()."/perl/ensembl/modules");
require Bio::EnsEMBL::Registry; 


my $registry = 'Bio::EnsEMBL::Registry';                                                                                                         
$registry->load_all();                                                                                                                           
my $slice_adaptor = $registry->get_adaptor('Human' ,'Core','Slice') or die;                     
my $gene_adaptor = $registry->get_adaptor('Human' ,'Core','Gene') or die;     

my @genes =  @{$gene_adaptor->fetch_all()};  


my %geneset;

foreach  my $gene (@genes){ 
    my $id = $gene->external_name();
    my $transcript = $gene->canonical_transcript() or next;
    my $translation = $transcript->translation() or next;
    my $pfeatures = $translation->get_all_ProteinFeatures() or next;
    while ( my $pfeature = shift @{$pfeatures} ) {
	my $logic_name = $pfeature->analysis()->logic_name();
	my $start  = $pfeature->start(); 
	my $end = $pfeature->end(); 
        my $ac = $pfeature->interpro_ac() or next;
        my $desc = $pfeature->idesc() or next;
	$geneset{"$desc($ac)"}{$id} = 1;
    }
}

foreach my $gs (keys %geneset){
    my @genes = keys %{$geneset{$gs}};
    print join("\t", ($gs, "", @genes))."\n";
}
