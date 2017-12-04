#!/usr/bin/perl

use strict;
use warnings;

use DBM::Deep;

my $db = new DBM::Deep('.puts.db');

my $d = $db->export();

foreach my $k(keys(%{$d})) {
    print join("\t", (
            $k,
            'PlantGDB PUTs',
            $d->{$k}->{'version'},
            $d->{$k}->{'build_date'},
            $d->{$k}->{'build_ests'},
            $d->{$k}->{'puts_count'},
            $d->{$k}->{'fasta_url'},
        ))."\n";
}
