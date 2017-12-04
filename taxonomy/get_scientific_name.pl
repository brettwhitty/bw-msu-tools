#!/usr/bin/perl

use lib "/home/whitty/SVN/lib"; 
use NCBITaxonomyTree; 

my $tax_id = shift @ARGV || die "Must provide a taxon ID";

my $tree = new NCBITaxonomyTree(path => "/projects/whitty_home/db/NCBI/Taxonomy_dump/");

my $scientific_name = $tree->get_scientific_name($tax_id) || die "Taxon ID not found";

print $scientific_name;
