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

foreach  my $gene (@genes){ 
    my $id = $gene->external_name();
    my $transcript = $gene->canonical_transcript() or next;
    my $start = $transcript->cdna_coding_start() or next;
    my $end = $transcript->cdna_coding_end() or next;
    my $length = $end - $start + 1;
    print $id."\t".$length."\n";
}
