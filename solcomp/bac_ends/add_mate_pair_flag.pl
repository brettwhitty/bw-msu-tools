#!/usr/bin/perl

##
## This will add a flag to the BAC_ends GFF3 if mate pairs have been mapped across BACs
##

$| = 1;

use strict;
use warnings;

use File::Temp qw{ tempfile };
use File::Copy qw{ move };
use Cwd qw{ abs_path };

my $gff_file = shift @ARGV || die "Please provide a GFF3 file as an argument";

$gff_file = abs_path($gff_file);

my $mates = {};

my $infh;

## first pass to check end pair mappings
open ($infh, $gff_file) || die "Failed to open '$gff_file': $!";
while (<$infh>) {
    ## don't allow script to be run twice
    if (/match_mate/) {
        die "Already has match_mate flags!";
    }
    /clone=([^;]+)/ || die "Failed to match on clone:\n$_";
    my $clone = $1;
    /end=([53?]+)/ || die "Failed to match on end:\n$_";
    my $end = $1;
    
    $mates->{$clone}->{$end} = 1;
}

my ($outfh, $outfile) = tempfile(UNLINK => 0);

## second pass to update the gff
open ($infh, $gff_file) || die "Failed to open '$gff_file': $!";
while (<$infh>) {
    /clone=([^;]+)/ || die "Failed to match on clone:\n$_";
    my $clone = $1;

    ## both ends are mapped
    if (scalar(keys(%{$mates->{$clone}})) > 1) {
        s/(clone=[^;]+);/$1;match_mate=1;/;
    } else {
        s/(clone=[^;]+);/$1;match_mate=0;/;
    }

    print $outfh $_;
}
close $outfh;

move($outfile, $gff_file);
