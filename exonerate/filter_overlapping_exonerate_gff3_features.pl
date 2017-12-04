#!/usr/bin/perl

##
## Quick script to filter exonerate GFF3 output for overlapping features
##
## coding2genome model produces alignments on both strands, so this will
## filter the alignments with the lower scores.
##
## This will also filter out duplicate features produced in the overlapping
## regions when the 'fastaoverlap' tool is used to split the target sequence.
##
## No measure for % overlap is used, because I only really want the best 
## alignment for a given query sequence within a given region anyway
##
## WARNING/NOTE: This script is meant only to remove duplicates of features 
##  (to self     for alignments of the *same* query sequence *aligned more 
##     and       than once in the same region*. This will not give the top
##   others)     hit amongst different aligned query sequences in the same 
##               region.
##
## Brett Whitty
## whitty@msu.edu
##

use strict;
use warnings;
use Carp;

#use Bio::Tools::GFF;
use Getopt::Long;

use lib "/home/whitty/SVN/lib";
use MyIO;

my ($input, $output);

my $result = GetOptions(
                           'input|i=s'      =>  \$input,
                           'output|o=s'     =>  \$output,
                       );

my $infh = get_infh($input);
my $outfh = get_outfh($output);

my $features = {};

while (<$infh>) {
    chomp;

    my @t = split("\t");

    my $type = $t[2];
    
    unless ($type eq 'cDNA_match' || $type eq 'gene') {
        next;
    }

    $t[8] =~ /Name=([^\;]+);/ || die "Failed parsing query name";
    my $query_acc = $1;
    $t[8] =~ /ID=([^\;]+);/ || die "Failed parsing ID";
    my $id = $1;
    
    my ($target_acc, $start, $end, $score) = ($t[0], $t[3], $t[4], $t[5]);
    
    push(@{$features->{$type}->{$target_acc}->{$query_acc}}, [ $start, $end, $score, $id ]);
}


my $delete_features = ();
foreach my $type(keys(%{$features})) {
    foreach my $target(keys(%{$features->{$type}})) {
        foreach my $query(keys(%{$features->{$type}->{$target}})) {
            my @coords = sort {$a->[0] <=> $b->[0]} @{$features->{$type}->{$target}->{$query}};
OUTER:          while (my $feat1 = shift @coords) {
                foreach my $feat2(@coords) {
                    if (features_overlap($feat1, $feat2)) {
                        if ($feat1->[2] > $feat2->[2]) {
                            $delete_features->{$target}->{$feat2->[3]} = 1;
                        } else {
                            $delete_features->{$target}->{$feat1->[3]} = 1;
                            next OUTER;
                        }
                    }
                }
            }
        }
    }
}

## open again to produce the output
$infh = get_infh($input);

while (<$infh>) {
    chomp;

    my @t = split("\t");

    my $target_acc = $t[0];
    
    my $id;
    if ($t[8] =~ /ID=([^\;]+);/) {
        $id = $1;
    }
    
    my $parent = '';
    if ($t[8] =~ /Parent=([^\;]+)(;)?/) {
        $parent = $1;
    }

    if (defined($id) && $delete_features->{$target_acc}->{$id}) {
        print STDERR "Skipping ID='$id'\n";
    } elsif (defined($parent) && $delete_features->{$target_acc}->{$parent}) {
        print STDERR "Skipping feat with parent='$parent'\n";
    } else {
        print $outfh $_."\n";
    }
}
    

sub features_overlap {
    my ($i1, $i2) = @_;
    ##
    ## |---i1---|
    ##            |---i2---|
    ##
    if ($i2->[0] > $i1->[1]) {
        return 0;
    ##
    ## |---i2---|
    ##            |---i1---|
    ##
    } elsif ($i2->[1] < $i1->[0]) {
        return 0;
    } else {
        print STDERR "***OVERLAP\n";
        return 1;
    }
}   
