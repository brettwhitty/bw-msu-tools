#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    chomp;

    if (/^>(\S+) (\S+([fFrR]{1})) .*clone (\S+),/) {

        my ($acc, $id, $end, $clone) = ($1, $2, $3, $4);
        $end = ($end =~ /^f$/i) ? '5' : '3';
        print join("\t", (
                $clone,
                $end,
                $id,
                $acc,
                         )
                  )."\n";

   } elsif (/^>(\S+) (\S+([fFrR]{1}[A-Z]{1})) .*clone (\S+),/) {

        my ($acc, $id, $end, $clone) = ($1, $2, $3, $4);
        $end = ($end =~ /^f$/i) ? '5' : '3';
        print join("\t", (
                $clone,
                $end,
                $id,
                $acc,
                         )
                  )."\n";

   } elsif (/^>(\S+) (\S+) .*clone (\S+),/) {

        my ($acc, $id, $end, $clone) = ($1, $2, '?', $3);
#        $end = ($end =~ /^f$/i) ? '5' : '3';
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
