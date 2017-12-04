#!/usr/bin/perl

$| = 1;

use strict;
use warnings;

while (<>) {
    chomp;

    /^(ORT[^(]+)\([^)]+\):\s+(.*)/ || die;

    my ($clid, $rest) = ($1, $2);
    my @t = split(" ", $rest);

    foreach my $member(@t) {
        $member =~ /^([^(]+)\(([^)]+)\)$/;
        my ($member_id, $source) = ($1, $2);
        print join("\t", ($member_id, $clid))."\n";
    }
}
