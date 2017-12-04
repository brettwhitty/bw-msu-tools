#!/usr/bin/perl

##
## Quick script to filter GFF to n tiers of features
##
## To be used with exonerate GFF, or other cDNA_match/gene feature type GFF
##
## Brett Whitty
## whitty@msu.edu
##

use strict;
use warnings;
use Carp;

use Bio::Tools::GFF;
use Getopt::Long;

use sort 'stable';

use lib "/home/whitty/SVN/lib";
use MyIO;
use TierFeatures;

my ($input, $output, $tiers);

my $result = GetOptions(
                           'input|i=s'      =>  \$input,
                           'output|o=s'     =>  \$output,
                           'tiers|t=i'      =>  \$tiers,
                       );
                      
$tiers = $tiers || 5; ## default to 5 tiers

my $infh = get_infh($input);
my $outfh = get_outfh($output);

print STDERR "WARNING: this script only supports filtering where parent features are:\n'cDNA_match', 'protein_match', 'gene', 'match'\n";

my $features = {};
#my @features = ();

my $tiered_feature = {};

while (<$infh>) {
    chomp;

    if (/^#/) { next; } ## skip comments
    
    if (/^##FASTA/) { last; }

    my @t = split("\t");

    my $type = $t[2];
    
    unless (
        $type eq 'cDNA_match' 
        || $type eq 'protein_match' 
        || $type eq 'gene' 
        || $type eq 'match'
        ) {
        next;
    }

    $t[8] =~ /Name=([^\;]+)/ || die "Failed parsing query name:\n$_";
    my $query_acc = $1;
    $t[8] =~ /ID=([^\;]+)/ || die "Failed parsing ID\n$_";
    my $id = $1;
    
    my ($target_acc, $start, $end, $score) = ($t[0], $t[3], $t[4], $t[5]);
    
    my $feature = new TierFeatures::Feature($start, $end, $id);
    $feature->{'chain_score'} = $score;

    push (@{$features->{$type}->{$target_acc}}, $feature);
}

foreach my $feat_type(keys(%{$features})) {
    foreach my $target_acc(keys(%{$features->{$feat_type}})) {

        my $feature_tierer = new TierFeatures($tiers);
    
        @{$features->{$feat_type}->{$target_acc}} 
            = reverse sort {$a->{'chain_score'} <=> $b->{'chain_score'}} 
              @{$features->{$feat_type}->{$target_acc}};
        
        my @tiers = $feature_tierer->tier_features(@{$features->{$feat_type}->{$target_acc}});

        my @tiered_features;
        for (my $i = 1; $i <= $tiers; $i++) {
            my $tier = shift @tiers;
            if (ref $tier) {
                push (@tiered_features, @$tier);
            } else {
                # no tiers left.  all done.
                last;
            }
        }
        ## set an array of flags
        foreach my $feature (@tiered_features) {
            $tiered_feature->{$target_acc}->{$feature->{'feat_id'}} = 1;
        }

    }
}

## open again to produce the output
$infh = get_infh($input);

while (<$infh>) {
    chomp;
    
    my @t = split("\t");
   

    ## REALLY NEED TO DEAL WITH FASTA HERE
    ## --> SPLIT FILE AFTER ##FASTA
    ## --> CREATE NEW FASTA DATABASE
    ## --> BUT REMOVE ID'S NOT IN OUTPUT TIERS
    
    ## comment or sequence line
    if (scalar(@t) != 9) {
        print $outfh $_."\n";
        next;
    }

    my $target_acc = $t[0];
    my $type = $t[2];
    
    my $id;
    if ($t[8] =~ /ID=([^\;]+);/) {
        $id = $1;
    }
    
    my $parent = '';
    if ($t[8] =~ /Parent=([^\;]+)(;)?/) {
        $parent = $1;
    }

    if (defined($id) && $tiered_feature->{$target_acc}->{$id}) {
        print $outfh $_."\n";
    } elsif (defined($parent) && $tiered_feature->{$target_acc}->{$parent}) {
        print $outfh $_."\n";
    } else {
#        print STDERR "Skipping ID='$id'\n";
    }
}
