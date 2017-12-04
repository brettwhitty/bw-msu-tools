#!/usr/bin/perl

$| = 1;

use lib "/home/whitty/SVN/lib";
use Sol::SeqDB;

my $accession = shift @ARGV || die "Provide a genbank accession";

my $db = new Sol::SeqDB();

my $gi = $db->get_gb_gi($accession, 1);
unless ($gi) {
    die "No record found for '$gi'";
}
my $genbank = $db->get_genbank($gi) or die "Failed fetching genbank file for '$gi'";
print $genbank;
