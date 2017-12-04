#!/usr/bin/perl

$| = 1;

use lib "/home/whitty/SVN/lib";
use Sol::SeqDB;

my $gi = shift @ARGV || die "Provide a genbank GI";

my $db = new Sol::SeqDB();

my $genbank = $db->get_genbank($gi) or die "Failed fetching genbank file for '$gi'";
print $genbank;
