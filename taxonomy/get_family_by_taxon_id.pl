#!/usr/bin/perl

use lib "/home/whitty/SVN/lib"; 
use NCBITaxonomyTree;

my $taxon_id = shift @ARGV;
my $rank = shift @ARGV || 'family';
my $flag = shift @ARGV || '';

my $tree = new NCBITaxonomyTree();

if ($tree->is_valid_taxon_id($taxon_id)) {
    if ($flag) {
        print $tree->get_taxon_id_for_rank_by_taxon_id($taxon_id, $rank) || 'N/A';
    } else {
        print $tree->get_name_for_rank_by_taxon_id($taxon_id, $rank) || 'N/A';
    }
} else {
    die "Invalid taxon ID '$taxon_id'!";
}
