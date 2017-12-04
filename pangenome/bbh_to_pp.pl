#!/usr/bin/perl

use strict;
use warnings;
use Carp;

## Quick script to take a OrthoMCL BLAST best-hit file (.bbh)
## and generate a phylogenetic profile table file
## suitable as input to the pangenome analysis pipeline
##
## Brett Whitty
## whitty@msu.edu

my $hits = {};
while (<>) {
    chomp $_;

    my @t = split("\t", $_);
    
    my $qdb = get_db($t[0]);
    my $sdb = get_db($t[1]);
    $hits->{$qdb}->{$t[0]}->{$sdb} = 1;
    $hits->{$sdb}->{$t[1]}->{$qdb} = 1;
}
my @dbs = sort {$a cmp $b} keys(%{$hits});
my $profile_size = scalar(@dbs);

print "# GENOME\tGENE\t".join("\t", @dbs)."\n";
foreach my $query_db(@dbs) {
    foreach my $gene(sort {$a cmp $b} keys(%{$hits->{$query_db}})) {
        my @profile = (0) x $profile_size;
        for (my $i = 0; $i < $profile_size; $i++) {
            my $subject_db = $dbs[$i];
            if ($subject_db eq $query_db) {
                $profile[$i] = 1;
                next;
            }
            if (defined($hits->{$query_db}->{$gene}->{$subject_db})) {
                $profile[$i] = 1;
            }
        }
        print join("\t", (
                $query_db,
                $gene,
                join("\t", @profile),
        ))."\n";
    }
}

## matches gene identifiers to return species names
sub get_db {
    my ($id) = @_;

    if ($id =~ /^jgi\|Sorbi1/) {
        return 'Sorghum_bicolor';
    } elsif ($id =~ /^jgi\|Phypa1_1/) {
        return 'Physcomitrella_patens';
    } elsif ($id =~ /^jgi\|Chlre4/) {
        return 'Chlamydomonas_reinhardtii';
    } elsif ($id =~ /^jgi\|Selmo1/) {
        return 'Selaginella_moellendorffii';
    } elsif ($id =~ /^Bradi\d+/) {
        return 'Brachypodium_distachyon';
    } elsif ($id =~ /^evm/) {
        return 'Carica_papaya';
    } elsif ($id =~ /^jgi\|Poptr1_1/) {
        return 'Populus_trichocarpa';
    } elsif ($id =~ /^GSVIVP/) {
        return 'Vitis_vinifera';
    } elsif ($id =~ /^AT[\dMC]G\d+/) {
        return 'Arabidopsis_thaliana';
    } elsif ($id =~ /^PGSC\d+DM/) {
        return 'Solanum_phureja';
    } elsif ($id =~ /^LOC_Os/) {
        return 'Oryza_sativa';
    } else {
        confess 'Unknown Species: '.$id;
    }

}
