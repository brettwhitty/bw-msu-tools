#!/usr/bin/perl

use lib "/home/whitty/SVN/lib"; 
use NCBITaxonomyTree;

my $taxon_id = shift @ARGV;

my $tree = new NCBITaxonomyTree();

if ($tree->is_valid_taxon_id($taxon_id)) {
    print $tree->get_common_name($taxon_id) || 'N/A';
} else {
    die "Invalid taxon ID '$taxon_id'!";
}
