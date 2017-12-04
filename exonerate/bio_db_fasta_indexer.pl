#!/usr/bin/perl

##
## Forces indexing of fasta databases for Bio::DB::Fasta access
## and sets them to read-only
##
## For use with scripts that use Bio::DB::Fasta on the grid to
## prevent them from a causing collisions during re-indexing

use strict;
use warnings;
use Bio::DB::Fasta;

my $file = shift @ARGV || die "Please provide path to fasta file";

unless (-f $file) {
    die "Please provide a path to a fasta file";
}

my $db = Bio::DB::Fasta->new($file, -reindex => 1);
my $indexname = $db->index_name($file);
chmod(0444, $indexname);
