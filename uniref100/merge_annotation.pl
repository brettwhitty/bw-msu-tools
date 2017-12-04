#!/usr/bin/perl

use strict;
use warnings;

my $annot_table = shift @ARGV;
my $t2g_map = shift @ARGV;

my $annot = {};
open (IN, $annot_table) || die $!;
while (<IN>) {
    chomp;

    my @t = split("\t", $_);

    $annot->{$t[0]} = $t[1];
}

my $t2g = {};
open (IN, $t2g_map) || die $!;
while (<IN>) {
    chomp;

    my @t = split("\t", $_);

    $t2g->{$t[0]} = $t[1];
}

my $gannot = {};
foreach my $id(keys(%{$t2g})) {
    my $func = $annot->{$id} || 'N/A';
    $gannot->{$t2g->{$id}}->{$func} = 1;
}

foreach my $id(keys(%{$gannot})) {
    if (keys %{$gannot->{$id}} > 1 && $gannot->{$id}->{'N/A'}) {
        delete $gannot->{$id}->{'N/A'};
    }
    if (keys %{$gannot->{$id}} > 1 && $gannot->{$id}->{'Gene of unknown function'}) {
        delete $gannot->{$id}->{'Gene of unknown function'};
    }
    if (keys %{$gannot->{$id}} > 1 && $gannot->{$id}->{'Conserved gene of unknown function'}) {
        delete $gannot->{$id}->{'Conserved gene of unknown function'};
    }
    $annot->{$id} = join("; ", keys(%{$gannot->{$id}}));
}

foreach my $id(sort keys %{$gannot}) {
    print join("\t", (
        $id,
        $annot->{$id},
    ))."\n";
}
