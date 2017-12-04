#!/usr/bin/perl

##
## This script converts the output from a nucmer run aligning
## BAC end sequences to BAC sequences into GFF3
##
## eg:
## 
## nucmer -p 4081.Solanum_lycopersicum.bacs.vs.4081.Solanum_lycopersicum.bac_ends 4081.Solanum_lycopersicum.bacs.fna 4081.Solanum_lycopersicum.bac_ends.fna
##

use strict;
use warnings;

#use File::Basename;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;

my $gff = new Bio::Tools::GFF(
                                -gff_version => '3',
                              );

my $delta_file = shift @ARGV || die "Please provide a delta file as an argument";

$delta_file =~ /(([^.]+)\.[^.]+\.bac_ends)\.delta/ || die "Unexpected delta file name";
my $base = $1;
my $taxon = $2;

#my $base = basename($delta_file, ".bacs_ends.delta");
my $ends_table = "$base.table";

my $clones = {};
my $aliases = {};
my $bac_end = {};

open (my $infh, $ends_table) || die "$ends_table: ".$!;
while (<$infh>) {
    chomp;

    my @t = split("\t", $_);
    
    $clones->{$t[3]} = $t[0];
    $aliases->{$t[3]} = $t[2];
    $bac_end->{$t[3]} = $t[1];

}

open (my $show_coords, "show-coords -I 95 -qcTH $delta_file |") || die $!;

my $ids = {};

while (<$show_coords>) {
    chomp;

    my @t = split("\t", $_);

    ## coverage
    if ($t[8] < 90) {
        next;
    }
    ## identity (This is handled by the -I 95 so it shouldn't matter)
    if ($t[6] < 95) {
        next;    
    }
    
    my $revcomp = 0;
    if ($t[2] > $t[3]) {
        $revcomp = 1;
    }

    my ($start, $end, $bac_acc, $end_acc, $coverage, $identity) 
        = ($t[0], $t[1], $t[9], $t[10], $t[8], $t[6]);

    my $strand = ($revcomp) ? '-' : '+';

    my $tags_ref = {
                        ID           => "match-$bac_acc-$end_acc-"
                                        . ++$ids->{$bac_acc}->{$end_acc},
                        Name         => $end_acc,
                        Alias        => $aliases->{$end_acc},
                        Clone        => $clones->{$end_acc},
                        End          => $bac_end->{$end_acc},
                        Coverage     => $coverage,
                        Identity     => $identity,
                        Dbxref       => "taxon:$taxon",
                   };

    my $match_feature = Bio::SeqFeature::Generic->new(
            -seq_id       => $bac_acc,
            -start        => $start,
            -end          => $end,
            -strand       => $strand,
            -primary      => 'match', # -primary_tag is a synonym
            -source_tag   => 'BAC_ends',
            -display_name => 'Display name',
#            -score        => '.',
            -tag          => $tags_ref,
                                                         );

        print $gff->gff_string($match_feature)."\n";

}
