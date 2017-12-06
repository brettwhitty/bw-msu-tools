#!/usr/bin/perl

use strict;
use warnings;
use Set::Scalar;

my $in1 = 'v3_arabidopsis_cluster_membership.txt';
my $in2 = 'v3_2_arabidopsis_cluster_membership.txt';

my $map1;
my $map2;

open my $infh, '<', $in1 || die;
while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    $map1->{$t[0]} = $t[1];
}

open $infh, '<', $in2 || die;
while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    $map2->{$t[0]} = $t[1];
}

my $set = new Set::Scalar(keys(%{$map1}), keys(%{$map2}));

print "## arabidopsis_accession\tv3_cluster_id\tv3.2_cluster_id\n";

foreach my $at_id(sort($set->members())) {
    print join("\t", (
        $at_id,
        $map1->{$at_id} || '-',
        $map2->{$at_id} || '-',
    ))."\n";
}

