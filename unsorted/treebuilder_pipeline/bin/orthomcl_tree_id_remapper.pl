#!/usr/bin/perl

$| = 1;

use strict;
use warnings;

use File::Basename;

#my $msa = shift @ARGV || die;
my $map_file = shift @ARGV || die;
my $newick = shift @ARGV || die;

my $dir = dirname($newick);
my $base = basename($newick, '.newick');

#$newick =~ /^([^.]+)/;
#my $pre = $1;
my $pre = $base;

my $out = "$dir/$pre.final.newick";
my $map = "$dir/$pre.codes.txt";

open (my $infh, '<', $map_file) || die "$!";

my $ids = {};
my $count = 0;
my $clid = '';
OUTER: while (<$infh>) {
    chomp;
    if (/^#(.*)/) {
        $clid = $1;
        if ($base ne $clid) {
            next;
        }
        while (<$infh>) {
            chomp;
            if (/^#/) {
                last OUTER;
            }
            my ($i, $g) = split("\t", $_);
            $ids->{$clid}->{$i} = $g;
        }
    }
}

open ($infh, '<', $newick) || die "$!";
open (my $outfh, '>', $out) || die "$!";
open (my $mapfh, '>', $map) || die "$!";
    foreach my $key(sort keys(%{$ids->{$pre}})) {
        my $id = $ids->{$pre}->{$key};
        print $mapfh "$key\t$id\n";
    }

while (<$infh>) {
    chomp;
    foreach my $key(keys(%{$ids->{$pre}})) {
        my $id = $ids->{$pre}->{$key};
        s/,$key:/,$id:/;
        s/\($key:/\($id:/;
        s/^$key:/$id:/;
    }
    print $outfh $_."\n";
}
