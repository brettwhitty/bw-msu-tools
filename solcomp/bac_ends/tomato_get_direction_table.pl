#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    chomp;

    if (/^>(\S+) (\S+).*clone (\S+) ([53]+)/) {

        my ($acc, $id, $clone, $end) = ($1, $2, $3, $4);
        print join("\t", (
                $clone,
                $end,
                $id,
                $acc,
                         )
                  )."\n";
    } elsif (/^>(\S+) (\S+).*(clone|similar to) (\S+)([fFrR]+)/) {
        my ($acc, $id, $clone, $end) = ($1, $2, $4, $5);
        $end = ($end =~ /^f$/i) ? '5' : '3';
        print join("\t", (
                $clone,
                $end,
                $id,
                $acc,
                         )
                  )."\n";
    } elsif (/^>/) {
        print STDERR $_."\n";
    }
}
