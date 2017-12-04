#!/usr/bin/perl

## converts swissprot flat files to fasta
## with optional taxonomy filter

use lib '/home/whitty/SVN/lib';
use NCBITaxonomyTree;

use Bio::SeqIO;

use strict;
use warnings;

my $filter_taxon = shift @ARGV;

my $tree;
if ($filter_taxon) {
#    $tree = new NCBITaxonomyTree(path => '/projects/whitty_home/db/NCBI/Taxonomy_dump');
    $tree = new NCBITaxonomyTree();
}

my $in  = Bio::SeqIO->new('-fh' => \*STDIN , '-format' => 'swiss');
my $out = Bio::SeqIO->new('-fh' => \*STDOUT, '-format' => 'fasta');
# note: we quote -format to keep older perl's from complaining.

my $counter=0;
print STDERR "Converting";
while ( my $seq = $in->next_seq() ) {
    my $species = $seq->species;

    if ($species->id == 325569 
        || $species->id == 9862 
        || $species->id == 50413 
        || $species->id == 119024) {
        next;
    }
    if ($filter_taxon) {
        if ($tree->id_has_ancestor_name($species->id, $filter_taxon)) {
            $out->write_seq($seq);
        }
    } else {
        $out->write_seq($seq);
    }
    $counter++;
    if ($counter % 100 == 0) {
        print STDERR "...$counter";
    }
}
print STDERR "...done.\n";
