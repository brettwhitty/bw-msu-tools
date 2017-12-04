#!/opt/rocks/bin/perl

## use bioperl to convert fgenesh/genscan output to gff
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;
use Carp;
use Getopt::Long;

use Bio::Tools::Fgenesh;
use Bio::Tools::Genscan;
use Bio::Tools::GFF;

my $method = undef;
my $input = undef;
my $output = undef;
my $gff_version = 3;
my $exons = 0;
my $utrs = 0;
my $noskip = 0;
my $spec_mode = 0;

my $result = GetOptions(
                            'method|m=s'    =>  \$method,
                            'input|i=s'     =>  \$input,
                            'output|o=s'    =>  \$output,
                            'version|v=i'   =>  \$gff_version,
                            'exons!'        =>  \$exons,
                            'utrs!'         =>  \$utrs,
                            'noskip!'       =>  \$noskip,
                            'spec!'         =>  \$spec_mode,
                       );

if ($spec_mode) {
    $utrs = 1;
    $noskip = 1;
    $exons = 1;
}    

## if adding other genefinders, create an empty hash ref here because
## I'm using this to check for supported types as well
my $term_map =  {
                    'Fgenesh'   =>  {
                                        #'Exon'          =>  'exon',
                                        'Promoter'      =>  'promoter',
                                        'SingletonExon'  =>  'CDS',
                                        'InitialExon'   =>  'CDS',
                                        'InternalExon'  =>  'CDS',
                                        'TerminalExon'  =>  'CDS',
                                        'Exon'          =>  'CDS',
                                        'Poly_A_site'   =>  'polyA_site',
                                    },

                    'Genscan'   =>  {
                        #                'Initial'       =>  'exon',
                        #                'Internal'      =>  'exon',
                        #                'Terminal'      =>  'exon',
                                        'Promoter'      =>  'promoter',
                                        'Exon'          =>  'CDS',
                                        'Initial'       =>  'CDS',
                                        'Internal'      =>  'CDS',
                                        'Terminal'      =>  'CDS',
                                        'Poly_A_site'   =>  'polyA_site',
                                    },                
                };

my $gene_features = {
                        'promoter'      =>  1,
                        'polyA_site'    =>  1,
                    };
                
my $mrna_features = {
                        'exon'  =>  1,
                        'CDS'   =>  1,
                    };


## this hash provides a way to skip features we don't want to output to GFF                    
my $skip_features = {
                        'promoter'      =>  1,
                        'polyA_site'    =>  1,
                    };
                    
unless ($method) {
    confess "Must provide value to --method flag, supported methods are: '".join("', '",keys(%{$term_map}))."'";
}                    
my $lcmethod = $method;
$lcmethod =~ tr/A-Z/a-z/;
unless (exists($term_map->{$method})) {
    confess "Method '$method' is unsupported, only supports '".join("', '",keys(%{$term_map}))."'";
}

unless ($input) {
    confess "Must specify an input file with the --input flag";
}
unless (-f $input) {
    confess "Input file '$input' does not exist";
}

my $gff3 = new Bio::Tools::GFF( 
                                #  -file => $outfile, 
                                -gff_version => 3, 
                              );

my $class = 'Bio::Tools::'.$method;
                              
my $parser = new $class( -file => $input );

my $outfh = get_outfh($output);

# parse the results
# note: this class is-a Bio::Tools::AnalysisResult which implements
# Bio::SeqAnalysisParserI, i.e., $fgensh->next_feature() is the same
my $gene_counter = 0;
my %feature_counter = (); 
while(my $gene = $parser->next_prediction()) {
  
    ## adjust the boundaries of the predicted gene/mRNA feature to exclude
    ## transcription start sites, poly-a signals, and any other UTRs
    ## unless the utrs flag is enabled
    unless ($utrs) {
        my @exons = $gene->exons_ordered();

        my ($first_exon, $last_exon) = ($exons[0], $exons[$#exons]);
        if ($exons[0]->location->strand == -1) {
            ($last_exon, $first_exon) = ($first_exon, $last_exon);
        }
        
        my $new_start = ($first_exon->location->start < $first_exon->location->end) ? $first_exon->location->start : $first_exon->location->end;
        my $new_end = ($last_exon->location->end > $last_exon->location->start) ? $last_exon->location->end : $last_exon->location->start;
        
        $gene->location->start($new_start);
        $gene->location->end($new_end);
    }
    
    my $seq_id = $gene->seq_id;
    my $gene_id = "$seq_id.$lcmethod.gene.".++$gene_counter;
    $gene->primary_tag($gene_id);
    
    $gene->source_tag($lcmethod);
    
    $gene->primary_tag('gene');
    $gene->add_tag_value("ID", $gene_id);
    $gene->add_tag_value("Name", $gene_id);
    
    print $outfh $gff3->gff_string($gene)."\n";
    
    $gene->primary_tag('mRNA');
    my $mrna_id = "$seq_id.$lcmethod.mRNA.$gene_counter";
    $gene->remove_tag('ID');
    $gene->add_tag_value("ID", $mrna_id);
    $gene->remove_tag('Name');
    $gene->add_tag_value("Name", $mrna_id);
    $gene->add_tag_value("Parent", $gene_id);
    
    print $outfh $gff3->gff_string($gene)."\n";
    
    my @features = $gene->features_ordered();
    
    foreach my $feature(@features) {
        $feature->source_tag($lcmethod);
        $feature->seq_id($seq_id);

        ## map feature terms if we need to
        if (exists($term_map->{$method}->{$feature->primary_tag})) {
            $feature->primary_tag($term_map->{$method}->{$feature->primary_tag});
        }            
        
        ## skip outputting this feature type if it's in the skip_feature hash
        if (! $noskip && $skip_features->{$feature->primary_tag}) {
            next;
        }
        
      
        my $feature_id = ''; 
        
             if ($feature->primary_tag eq 'CDS') {#
            $feature_id = "$seq_id.$lcmethod.CDS.$gene_counter";#
        } else {#
            $feature_id = "$seq_id.$lcmethod.".$feature->primary_tag.".".++$feature_counter{$feature->primary_tag};
            }#
        
        if ($gene_features->{$feature->primary_tag}) {
            $feature->add_tag_value("Parent", $gene_id);
        } elsif ($mrna_features->{$feature->primary_tag}) {
            $feature->add_tag_value("Parent", $mrna_id);
        }
        $feature->add_tag_value("ID", $feature_id);
        
        print $outfh $gff3->gff_string($feature)."\n";
        if ($feature->primary_tag eq 'CDS' && $exons) {
            $feature->frame('.');
            $feature->primary_tag('exon');
            $feature_id = "$seq_id.$lcmethod.".$feature->primary_tag.".".++$feature_counter{$feature->primary_tag};
            $feature->remove_tag('ID');
            $feature->add_tag_value('ID', $feature_id);
            print $outfh $gff3->gff_string($feature)."\n";
        }
    }
}
$parser->close();

## returns an output file handle
sub get_outfh {
    my ($file) = @_;
   
    if ($file) {
        open (my $fh, ">$file") || confess "Failed opening '$file' for writing: $!";
        return $fh;
    } else {
        return \*STDOUT;
    }
}
