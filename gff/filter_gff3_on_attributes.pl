#!/usr/bin/perl

=head1 NAME

filter_gff3_on_attributes - Remove features from the GFF3 based on the values of their attributes

=head1 SYNOPSIS

Usage:
  filter_gff3_on_attributes -i input.gff3 -o output.gff3 -f feature_type

=head1 DESCRIPTION

This script will filter features in a GFF3 file based on their attribute values.

The options are:
   
   -f/--feature     - parent feature type(s) to filter on (eg: 'match' for bp_search2gff 
                      or 'gene,cDNA_match' for exonerate raw output), can be more than one delimited
                      with a comma
   -i/--input       - input GFF3 file
   -o/--output      - output GFF3 file

   --skip           - optionally skip certain feature types and don't output 
                      (can be a comma delimited list)

   --id_cutoff      - percent identity cutoff (eg: 70) using 'identity' attribute
   --sim_cutoff     - percent similarity cutoff (eg: 70) using 'similarity' attribute
   --cov_cutoff     - percent coverage cutoff (eg: 50) using 'coverage' attribute
   --sig_cutoff     - significance cutoff (eg: 1e-5) using 'significance' attribute
   --rank_cutoff    - rank cutoff (eg: 1) using 'rank' attribute
   --len_cutoff     - length cutoff, using start/end coordinates of feature
   --top            - keep only the top X features (eg: X = 1) by score
   --skip_self      - skip features where Name attribute value equals the sequence id
   
=head1 AUTHOR

Brett Whitty, whitty@msu.edu

=cut

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;

use lib "/home/whitty/SVN/lib";
use MyIO;

my $id_cutoff;
my $sim_cutoff;
my $cov_cutoff;
my $sig_cutoff;
my $rank_cutoff;
my $feature;
my $input;
my $output;
my $top_n; ## top n hits
my $len_cutoff;
my $help;
my $skip;
my $skip_self;

my $result = GetOptions(
                            'id_cutoff=i'   =>  \$id_cutoff,
                            'sim_cutoff=i'  =>  \$sim_cutoff,
                            'cov_cutoff=i'  =>  \$cov_cutoff,
                            'sig_cutoff=s'  =>  \$sig_cutoff,
                            'rank_cutoff=i' =>  \$rank_cutoff,
                            'len_cutoff=i'  =>  \$len_cutoff,
                            'top=i'         =>  \$top_n,
                            'skip_self=i'   =>  \$skip_self,
                            'feature|f=s'   =>  \$feature,
                            'skip|s=s'      =>  \$skip,
                            'input|i=s'     =>  \$input,
                            'output|o=s'    =>  \$output,
                            'help|h!'       =>  \$help,
                            'man|m!'        =>  \$help,
                      );

if ($help) {
    pod2usage(-verbose => 2);
}
                      
$id_cutoff  = (defined($id_cutoff))  ?  $id_cutoff  : 0;
$sim_cutoff = (defined($sim_cutoff)) ?  $sim_cutoff : 0;
$cov_cutoff = (defined($cov_cutoff)) ?  $cov_cutoff : 0;
$sig_cutoff = (defined($sig_cutoff)) ?  $sig_cutoff : '10';
$feature    = (defined($feature))    ?  $feature    : 'match';
$skip       = (defined($skip))       ?  $skip       : '';
$top_n      = (defined($top_n))      ?  $top_n      : 0;
$skip_self  = (defined($skip_self))  ?  $skip_self  : 0;
$len_cutoff = (defined($len_cutoff)) ?  $len_cutoff : 0;
                      
if ($input && $output && $input eq $output) {
    pod2usage();
    confess "Input file '$input' and output file '$output' are the same file";
}
    
my $parent_features = {};
foreach my $feat_type(split(",", $feature)) {
    $parent_features->{$feat_type} = 1;
} 
my $skip_features = {};
foreach my $skip_type(split(",", $skip)) {
    $skip_features->{$skip_type} = 1;
} 

my $delete_feature = {};

my $last_score  = {};
my $top_counter = {};
my $feat_scores = {};


my $infh = get_infh($input);
while (<$infh>) {
    if (/^#/) {
        next;
    }
    chomp;
    
    my @cols = split("\t", $_);
   
    if (scalar(@cols) != 9) { 
        last; 
    }
    
    my $feature_type = $cols[2];
 
    ## skip features that aren't the parent feature types specified by the 'feature' option
    unless ($parent_features->{$feature_type}) { 
        next;
    }
   
    my ($id, $parent, $coverage, $identity, $similarity, $significance, $rank, $name);
    
    my $seq_id = $cols[0];
    
    my $score = $cols[5];

    my $start = $cols[3];
    my $end   = $cols[4];
    my $len = $end - $start + 1; 
    
    if ($cols[8] =~ /ID=([^;]+)/) {
        $id = $1;
    }

    unless ($id) {
        confess "Script requires that an ID attribute be present on the feature to filter, but none found.";
    }
    if ($cols[8] =~ /Parent=([^;]+)/i) {
        $parent = $1;
    }
    if ($cols[8] =~ /identity=([^;]+)/i) {
        $identity = $1;
    }
    if ($cols[8] =~ /similarity=([^;]+)/i) {
        $similarity = $1;
    }
    if ($cols[8] =~ /coverage=([^;]+)/i) {
        $coverage = $1;
    }
    if ($cols[8] =~ /significance=([^;]+)/i) {
        $significance = $1;
    }
    if ($cols[8] =~ /rank=([^;]+)/i) {
        $rank = $1;
    }
    if ($cols[8] =~ /Name=([^;]+)/i) {
        $name = $1;
    }

    ## optionally skip features that have the same Name and sequence_id 
    if ($skip_self && $name && $name eq $cols[0]) {
        $delete_feature->{$id} = 1;
        next;
    }
   
    ## this code assumes match features are sorted in GFF3 from highest to lowest scores
    ## this is true for output from bp_search2gff.pl
#    if (!defined($last_score->{$id})) {
#        $last_score->{$id} = 0;
#        $top_counter->{$id} = 1;
#    }
#    if ($score > $last_score->{$id}) {
#        if ($last_score->{$id} > 0) {
#            $top_counter->{$id}++;
#        }
#        $last_score->{$id} = $score;
#    }

#    if ($top_n && $top_counter->{$id} > $top_n) {
#        $delete_feature->{$id} = 1;
#    }
    
    if ($cov_cutoff && $coverage && $coverage < $cov_cutoff) {
        $delete_feature->{$id} = 1;
    }
    if ($id_cutoff && $identity && $identity < $id_cutoff) {
        $delete_feature->{$id} = 1;
    }
    if ($sim_cutoff && $similarity && $similarity < $sim_cutoff) {
        $delete_feature->{$id} = 1;
    }
    if ($sig_cutoff && $significance && $significance > $sig_cutoff) {
        $delete_feature->{$id} = 1;
    }
    if ($rank_cutoff && $rank && $rank > $rank_cutoff) {
        $delete_feature->{$id} = 1;
    }
    if ($len_cutoff && $len < $len_cutoff) {
        $delete_feature->{$id} = 1;
    } 
   
    unless($delete_feature->{$id}) {
        $feat_scores->{$seq_id}->{$id}->{'score'} = $score;
        $feat_scores->{$seq_id}->{$id}->{'evalue'} = $significance;
    }
    
}

## do deleting for > top n features
if ($top_n) {
    foreach my $seq_id(keys(%{$feat_scores})) {
        my $last_score = 0;
        my $top_counter = 0;
    
        my @hits = sort {
            $feat_scores->{$seq_id}->{$b}->{'score'} <=> $feat_scores->{$seq_id}->{$a}->{'score'}
                        } keys %{$feat_scores->{$seq_id}};
    
        foreach my $hit(@hits) {
            my $score = $feat_scores->{$seq_id}->{$hit}->{'score'};
            if ($score > $last_score) {
                $top_counter++;
            }
            if ($top_n && $top_counter > $top_n) {
                $delete_feature->{$hit} = 1;
            }
        }
    }
}

my $outfh = get_outfh($output);

$infh = get_infh($input);
while (<$infh>) {
    chomp;
    
    my @cols = split("\t", $_);
    
    if (scalar(@cols) == 9) { 
    
        my ($id, $parent);

        my $feature_type = $cols[2];
        ## optionally skip certain feature types
        if ($skip_features->{$feature_type}) {
            next;
        }
 
        if ($cols[8] =~ /ID=([^;]+)/) {
            $id = $1;
        }
        if ($cols[8] =~ /Parent=([^;]+)/) {
            $parent = $1;
        }   

        if ($id) {
            if ($delete_feature->{$id}) {
                next;
            }
        }
        if ($parent) {
            if ($delete_feature->{$parent}) {
                next;
            }
        }
    }
    print $outfh $_."\n";    
}


