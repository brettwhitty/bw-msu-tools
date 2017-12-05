#!/usr/bin/perl

$| = 1;

use lib "/home/whitty/SVN/lib";
use Sol::SeqDB;

my $accession = shift @ARGV || die "Provide a genbank accession";

my $db = new Sol::SeqDB();

my $gi = $db->get_gb_gi($accession);
my $fasta = $db->get_gb_fasta($gi) or die "Failed fetching fasta for '$gi'";
$fasta =~ s/>gi\|\d+\|gb\|([^.]+)[^\|]+\|/>$1/g;
print $fasta;
