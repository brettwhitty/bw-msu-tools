#!/usr/bin/perl

$| = 1;

use lib "/home/whitty/SVN/lib";
use Sol::SeqDB;

my $db = new Sol::SeqDB();

while (<>) {
    chomp;

    my @t = split("\t", $_);
    my $gi = $t[4];

    my $genbank = $db->get_genbank($gi) or die "Failed fetching genbank file for '$gi'";
    if ($genbank =~ /Buell/) {
        print $_."\n";
    }
}
