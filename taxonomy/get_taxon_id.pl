#!/usr/bin/perl

use lib "/home/whitty/SVN/lib"; 
use NCBITaxonomyTree; 

my $species = shift @ARGV || die "Must provide a species name";

$species =~ s/subsp_/subsp\./;
$species =~ s/_/ /g;

my $tree = new NCBITaxonomyTree(path => "/home/whitty/db/NCBI/");

my $taxon_id = $tree->get_taxon_id($species) || die "Taxon name not found";

print $taxon_id;
