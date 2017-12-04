#!/usr/bin/perl

use strict;
use warnings;
use Carp;

#Chr1    genoscope   chromosome      1       30432563        .       .       .       ID=Chr1;Name=Chr1

my $infile = shift @ARGV;
my $source = shift @ARGV || 'genoscope';
my $feature = shift @ARGV || 'chromosome';

open (my $seqstat_fh, "seqstat -a $infile |") || confess "Failed to open pipe to seqstat on '$infile': $!";

my $lens = {};
while (<$seqstat_fh>) {
    chomp;
    if (/^\*\s+(\S+)\s+(\d+)/) {
        my ($seq_id, $len) = ($1, $2);
        my $chromosome;
        if ($seq_id =~ /^chr(\d+)/) {
            $chromosome = $1;
#            $seq_id = ucfirst($seq_id);
        }
        print join("\t", (
                            $seq_id,
                            $source,
                            $feature,
                            1,
                            $len,
                            '.',
                            '.',
                            '.',
                            'ID='.$seq_id.';Name='.$seq_id,
                         ))."\n";
    }
}
