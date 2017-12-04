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

use Bio::DB::Fasta;
use Bio::Tools::GFF;
use Getopt::Long;

use sort 'stable';

use lib "/home/whitty/SVN/lib";
use MyIO;
use TierFeatures;

my ($input, $output, $tiers, $query_fasta, $id_cutoff, $cov_cutoff);

my $result = GetOptions(
                           'input|i=s'      =>  \$input,
                           'output|o=s'     =>  \$output,
                           'tiers|t=i'      =>  \$tiers,
                           'query-fasta|q=s'    =>  \$query_fasta,
                           'id_cutoff=i'    =>  \$id_cutoff,
                           'cov_cutoff=i'   =>  \$cov_cutoff,
                       );
                      
$tiers = $tiers || 5; ## default to 5 tiers

my $infh = get_infh($input);
my $outfh = get_outfh($output);


my $features = {};
#my @features = ();
my $child_mrna = {};

my $tiered_feature = {};

my $qfdb = undef;
if (defined($query_fasta)) {
    $qfdb = new Bio::DB::Fasta($query_fasta);
}

while (<$infh>) {
    chomp;

    if (/^#/) { next; } ## skip comments
    
    if (/^##FASTA/) { last; }

    my @t = split("\t");

    my $type = $t[2];
    
    #unless ($type eq 'cDNA_match' || $type eq 'gene') {
    unless ($type eq 'mRNA') {
        next;
    }

    $t[8] =~ /Name=([^;]+)/ || die "Failed parsing query name:\n$_";
    my $query_acc = $1;
    $t[8] =~ /ID=([^;]+)/ || die "Failed parsing ID\n$_";
    my $id = $1;
    $t[8] =~ /Parent=([^;]+)/ || die "Failed parsing ID\n$_";
    my $parent = $1;
    $t[8] =~ /Coverage=([^;]+)/i || die "Failed parsing Coverage:\n$_";
    my $coverage = $1;
    $t[8] =~ /Identity=([^;]+)/i || die "Failed parsing Identity:\n$_";
    my $identity = $1;

    if (defined($id_cutoff) && $identity < $id_cutoff)   { next; }
    if (defined($cov_cutoff) && $coverage < $cov_cutoff) { next; }

    my ($target_acc, $start, $end, $score) = ($t[0], $t[3], $t[4], $t[5]);
    
    $child_mrna->{$parent} = $id;

    my $feature = new TierFeatures::Feature($start, $end, $id);
    #$feature->{'chain_score'} = $score;
    my $feat_chain_score = (($coverage / 100) * $identity) / 100; ## perfect = 100%

    ## if a fasta file is provided, use the length of the sequence
    ## to calculate a feature score for tiering (best overall alignment
    ## as total # of identical aligned bases)
    if (defined($qfdb)) {
        $feat_chain_score *= $qfdb->length($query_acc);
    }

    $feature->{'chain_score'} = $feat_chain_score;

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

my $last = ''; 

while (<$infh>) {
    chomp;

    if (/^#/) {
        if ($last !~ /^###$/) {
            $last = $_;
            print $outfh $_."\n";
        }
        next;
    }
    if (/^##FASTA/) { last; }    

    my @t = split("\t");

    ## REALLY NEED TO DEAL WITH FASTA HERE
    ## --> SPLIT FILE AFTER ##FASTA
    ## --> CREATE NEW FASTA DATABASE
    ## --> BUT REMOVE ID'S NOT IN OUTPUT TIERS

    my $target_acc = $t[0];
    my $type = $t[2];
    
    my $id;
    if ($t[8] =~ /ID=([^\;]+)/) {
        $id = $1;
    }
    
    my $parent = '';
    if ($t[8] =~ /Parent=([^\;]+)/) {
        $parent = $1;
    }

    if (defined($id) && defined($tiered_feature->{$target_acc}->{$id})) {
        $last = $_;
        print $outfh $_."\n";
    } elsif ($type eq 'gene' && defined($id) && defined($tiered_feature->{$target_acc}->{$child_mrna->{$id}})) {
        $last = $_;
        print $outfh $_."\n";
    } elsif (defined($parent) && defined($tiered_feature->{$target_acc}->{$parent})) {
        $last = $_;
        print $outfh $_."\n";
    } else {
#        print STDERR "Skipping ID='$id'\n";
    }
}
