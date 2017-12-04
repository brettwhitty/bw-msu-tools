#!/usr/bin/perl

use lib "/home/whitty/SVN/lib"; 
use NCBITaxonomyTree; 

my $dir = shift @ARGV;

$dir =~ /^(\d+)\..*/ || die "Failed to match tax_id in dir name";

$tax_id = $1;

my $tree = new NCBITaxonomyTree(path => "/home/whitty/db/NCBI/"); 
if ($tree->id_has_ancestor_name($tax_id, "Solanaceae")) { 
    print "true";
} else { 
    print "false";
}
