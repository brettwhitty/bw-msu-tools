#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    chomp;

    my @t = split(/\t/, $_);

    my @ids = split(/\s+/, $t[4]);

    foreach my $id(@ids) {
        if ($id =~ /^AT/i) {
            print join("\t", (
                    $id,
                    $t[0],
            ))."\n";
        }
    }
}

