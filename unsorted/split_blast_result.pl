#!/opt/rocks/bin/perl

use strict;
use warnings;

use Carp;
use POSIX;
use Getopt::Long;

## quick script to split wu-blast records, probably also works for NCBI BLAST

## wu-blast record separator
$/ = "";

## flush
$| = 1;

my ($input, $output, $records);

GetOptions(
    'input|i=s'             =>  \$input,
    'output_prefix|o=s'     =>  \$output,
    'records|r=i'           =>  \$records,
);

if (! defined($input) && ! -e $input) {
    confess "Need raw BLAST input file with --input flag";
}

$records ||= 10000;

open my $infh, '<', $input or confess "Failed to open '$input' for reading: $!";

my $output_prefix = '';
if (defined($output)) {
    $output_prefix = $output;
} else {
    $output_prefix = $input;
}

my $filen = 0;
my $counter = 0;
my $outfh;
while (<$infh>) {
    $counter++;

    my $curn = floor(($counter - 1) / $records) + 1;

    if ($curn != $filen) {
        $filen = $curn;
        print STDERR "Writing to file '$output_prefix.$filen'...\n";
        open $outfh, '>', $output_prefix.'.'.$filen
            or confess "Failed to open '$output_prefix.$filen' for writing: $!";
    }

    print $outfh $_;
}
