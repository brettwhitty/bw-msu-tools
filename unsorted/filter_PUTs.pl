#!/usr/bin/perl

## filter PUTs files based on length and # of Ns

use strict;
use warnings;

use Bio::DB::Fasta;
use Cwd qw{ abs_path };
use File::Basename qw{ basename dirname}; 

my $min_length = 250;
my $max_n = 10;

my $file = shift @ARGV || die "Provide a fasta file";

$file = abs_path($file);

unless (-e $file) { die "File '$file' doesn't exist"; }

my $outdir = dirname($file);
my $base = basename($file, ".fasta");
my $outfile = "$outdir/$base.filter.fasta";

my $db = Bio::DB::Fasta->new($file, -reindex => 1);

open (my $outfh, ">$outfile") || die "failed to open '$outfile' for writing: $!";

my @ids = $db->ids();
foreach my $id(@ids) {
    my $len = $db->length($id);

    if ($len >= $min_length) {
        my $header = ">".$db->header($id)."\n";
        my $seq = $db->seq($id);
        my $n_count = 0;
        while ($seq =~ /n/ig) {
            $n_count++;
        }
        if ($n_count > $max_n) {
            print STDERR "$id\tN=$n_count\n";
        } else {
            $seq =~ s/(.{1,60})/$1\n/g;
            print $outfh $header.$seq;
        }
    } else {
        print STDERR "$id\tlen=$len\n";
    }
}
