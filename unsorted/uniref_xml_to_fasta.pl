#!/usr/bin/perl

## quick script to convert UniProt UniRef XML into FASTA
## 
## also optionally removes non-standard AAs


use strict;
use warnings;

use lib "/home/whitty/SVN/lib";
use NCBITaxonomyTree;

use XML::Twig;

my $xml_file = shift @ARGV || die "Provide a UniRef XML file as input";
my $clean_aa = shift @ARGV;
my $taxon_filter = shift @ARGV; ## provide the name of a taxonomy node to filter on eg: 'Solanaceae' or 'Solanum tuberosum'

my $tree = new NCBITaxonomyTree();

my $infh;
if ($xml_file =~ /\.gz$/) {
    open $infh, '<:gzip', $xml_file or die $!;
} else {
    open $infh, "<$xml_file" or die $!;
}

my $twig = new XML::Twig(
                            twig_handlers => {
                                                entry => \&entry_handler,
                                             }
                        );

$twig->parse($infh);

sub entry_handler {
    my ($twig, $entry) = @_;
    
    my $id = $entry->att('id');
    my $name = $entry->first_child('name')->text;
    my $sequence = $entry->first_child('representativeMember')->first_child('sequence')->text;
    
    my $taxon_id = '';
    foreach my $child($entry->first_child('representativeMember')->first_child('dbReference')->children) {
        if ($child->att('type') eq 'NCBI taxonomy') {
            $taxon_id = $child->att('value');
        }
    }

    ## allow filtering on taxonomy
    if (! $taxon_filter || ($tree->is_valid_taxon_id($taxon_id) && $tree->id_has_ancestor_name($taxon_id, $taxon_filter))) {

        $id =~ s/^\s+|\s+$//g;
        $sequence =~ s/^\s+|\s+$//g;
        $name =~ s/^\s+|\s+$//g;

        ## optionally remove exotic AAs that wu-blast (among others) can't deal with very well
        if ($clean_aa) {
            $sequence =~ tr/jJoOuU/iIkKcC/;
        }

        $name =~ s/Cluster: //;

        print ">$id $name\n$sequence\n";
    
    }

    $twig->purge();
}    
