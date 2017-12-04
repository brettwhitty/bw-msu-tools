#!/usr/bin/perl

## use bioperl to convert AAT btab output to gff
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;
use Carp;

use lib "/home/whitty/SVN/lib";
use TierFeatures;

use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use Getopt::Long;

my $method = 'aat';
my $input = undef;
my $output = undef;
my $gff_version = 3;
my $infh = undef;
my $min_percent_id = 0;
my $min_percent_sim = 0;
my $min_percent_cov = 0;
my $tiers = 0;
my $add_note = 0;

my $db_path = '';

my $feature_match       = 'match';
my $feature_match_part  = 'match_part';

my $result = GetOptions(
                            'method|m=s'            =>  \$method,
                            'match=s'               =>  \$feature_match,
                            'match_part=s'          =>  \$feature_match_part,
                            'input|i=s'             =>  \$input,
                            'output|o=s'            =>  \$output,
                            'version|v=s'           =>  \$gff_version,
                            'min_percent_id|p=s',   =>  \$min_percent_id,
                            'min_percent_sim|s=s',  =>  \$min_percent_sim,
                            'min_percent_cov|c=s',  =>  \$min_percent_cov,
                            'tiers|t=i'             =>  \$tiers,
                            'note|n!'               =>  \$add_note,
                            'db_path|d=s'           =>  \$db_path,
                       );

unless ($input) {
    confess "Must specify input file with --input flag";
}
unless (-f $input) {
    confess "Must specify input file that exists with --input flag";
}
                       
my $lcmethod = $method;
$lcmethod =~ tr/A-Z/a-z/;

my $gff = new Bio::Tools::GFF( 
                                -gff_version => $gff_version, 
                              );

## convert aat btab strand text to Bio::SeqFeature                             
my %strand = (  
                'Plus'  => '+1',
                'Minus' => '-1'
             );
                              
## get input file handle
$infh = get_infh($input);

## do a pass through the file to get the ranges for the parent match features
## doing this to avoid having to keep a lot of things in memory
my $match = {};
my $percent_id = {};
my $percent_sim = {};
my $percent_cov = {};
my $tiered_feature = {};
my $subject_db = undef;
while (<$infh>) {
    my @btab = split(/\t/, $_);
   
    $subject_db = undef;

    ## clean up query sequence id
    my $seq_id = get_clean_accession($btab[0]);
   
    ## check if subject database can be found
    ## then we'll be using to to calculate the percent coverage
    if (defined($db_path) && -e $db_path) {
        $subject_db = new Bio::DB::Fasta( $db_path );
    } elsif (-e $btab[4]) {
        $subject_db = new Bio::DB::Fasta( $btab[4] );
    }

    ## subject accession
    my $subject_accession = $btab[5];

    ## chain index #
    my $chain_id = $btab[13];
    
    ## score
    my $hsp_score = $btab[9];
    
    ## percent id
    my $hsp_pct_id = $btab[10];

    ## percent similarity
    my $hsp_pct_sim = $btab[11];
    
    ## set the strand for the match
    unless(defined($match->{$seq_id}->[$chain_id]->{'strand'})) { 
        $match->{$seq_id}->[$chain_id]->{'strand'} = $strand{$btab[17]};
    }
    
    ## sum the hsp scores to give the match feature score
    $match->{$seq_id}->[$chain_id]->{'score'} += $hsp_score;

    ## start and end positions are flipped for Minus strand
    my ($start, $end) = ($btab[6], $btab[7]);
    if ($start > $end) {
        ($start, $end) = ($end, $start);
    }

    my $hsp_len = $end - $start + 1;

    ## calculate hsp lengths for subject to use for coverage calculations
    my ($subject_start, $subject_end) = ($btab[8], $btab[9]);
    if ($subject_start > $subject_end) {
        ($subject_start, $subject_end) = ($subject_end, $subject_start);
    }
    my $subject_hsp_len = $subject_end - $subject_start + 1;


    $percent_id->{$seq_id}->[$chain_id]->{'pctsum'} += $hsp_len * $hsp_pct_id; 
    $percent_sim->{$seq_id}->[$chain_id]->{'pctsum'} += $hsp_len * $hsp_pct_sim; 
    $percent_id->{$seq_id}->[$chain_id]->{'lensum'} += $hsp_len; 
    $percent_sim->{$seq_id}->[$chain_id]->{'lensum'} += $hsp_len; 
    
    ## for coverage calculations
    if (defined($subject_db)) {
        $percent_cov->{$seq_id}->[$chain_id]->{'subject_accession'} = $subject_accession;
        $percent_cov->{$seq_id}->[$chain_id]->{'lensum'} += $subject_hsp_len;
    }
    
    my $match_start = ( $match->{$seq_id}->[$chain_id]->{'match'}->[0] ) 
                      ? $match->{$seq_id}->[$chain_id]->{'match'}->[0] : 999999999;
                      
    my $match_end = ( $match->{$seq_id}->[$chain_id]->{'match'}->[0] ) 
                    ? $match->{$seq_id}->[$chain_id]->{'match'}->[0] : -1;
    
    ## 0 = start, 1 = end
    if ($start < $match_start) {
        $match->{$seq_id}->[$chain_id]->{'match'}->[0] = $start;
    }
    if ($end > $match_end) {
        $match->{$seq_id}->[$chain_id]->{'match'}->[1] = $end;
    }
}

## calculate percent identity for chains
foreach my $seq_id(keys(%{$percent_id})) {
    foreach my $chain_ref(@{$percent_id->{$seq_id}}) {
        unless (defined($chain_ref)) {
            next;
        }
        my $pctsum = $chain_ref->{'pctsum'};
        my $lensum = $chain_ref->{'lensum'};
        $chain_ref = $pctsum / $lensum;
    }
}    

## calculate percent similarity for chains
foreach my $seq_id(keys(%{$percent_sim})) {
    foreach my $chain_ref(@{$percent_sim->{$seq_id}}) {
        unless (defined($chain_ref)) {
            next;
        }
        my $pctsum = $chain_ref->{'pctsum'};
        my $lensum = $chain_ref->{'lensum'};
        $chain_ref = $pctsum / $lensum;
    }
}    

if (defined($subject_db)) {
    ## calculate percent coverage for chains
    foreach my $seq_id(keys(%{$percent_cov})) {
        foreach my $chain_ref(@{$percent_cov->{$seq_id}}) {
            unless (defined($chain_ref)) {
                next;
            }
            my $subject_accession = $chain_ref->{'subject_accession'};
            my $subject_len = $subject_db->length($subject_accession);
            my $lensum = $chain_ref->{'lensum'};
            $chain_ref = $lensum / $subject_len * 100;
        }
    }
}    

## if # of tiers has been specified then filter on tiers
if ($tiers) {
    foreach my $seq_id(keys(%{$match})) {
        my @features;
        my $chain_count = scalar(@{$match->{$seq_id}});
        for (my $chain_id = 0; $chain_id < $chain_count; $chain_id++) {

            my $chain_ref = $match->{$seq_id}->[$chain_id];
            
            my $match_start = $chain_ref->{'match'}->[0];
            my $match_end   = $chain_ref->{'match'}->[1];
            my $score       = $chain_ref->{'score'};

            ## skip empty records                                                            
            unless ($match_start && $match_end) {
                next;
            }
            ## filter on percent identity
            if (defined($percent_cov) && $percent_cov->{$seq_id}->[$chain_id] < $min_percent_cov) {
                next;
            }
            
            ## filter on percent identity
            if ($percent_id->{$seq_id}->[$chain_id] < $min_percent_id) {
                next;
            }
    
            ## filter on percent similarity 
            if ($percent_sim->{$seq_id}->[$chain_id] < $min_percent_sim) {
                next;
            }
 
            my $feature = new TierFeatures::Feature($match_start, $match_end, $chain_id);
        
            $feature->{chain_score} = $score;

            push (@features, $feature);

        }

        my $feature_tierer = new TierFeatures($tiers);
    
        ## sort features by score, desc, so best features are tiered first:
        @features = reverse sort {$a->{chain_score} <=> $b->{chain_score}} @features; 

        my @tiers = $feature_tierer->tier_features(@features);

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
            $tiered_feature->{$seq_id}->[$feature->{'feat_id'}] = 1;
        }
    }
}

## get input file handle
$infh = get_infh($input);

## get output file handle
my $outfh = get_outfh($output);

my $feature_count = {};
my $parent_id = {};

my $match_to_gff = {};
while (<$infh>) {
    chomp;
    my @btab = split(/\t/, $_);
  
    ## clean up query sequence id
    my $seq_id = get_clean_accession($btab[0]);
   
    ## score
    my $hsp_score = $btab[9];

    ## percent id
    my $hsp_pct_id = $btab[10];

    ## percent similarity
    my $hsp_pct_sim = $btab[11];
    
    ## chain index #
    my $chain_id = $btab[13];
    
    ## hsp index #
    my $hsp_id = $btab[14];
   
    ## subject accession
    my $subject_accession = $btab[5];

    ## subject header
    my $note = get_clean_note($btab[15]);
    
    ## filter on # of tiers if enabled
    if ($tiers && ! $tiered_feature->{$seq_id}->[$chain_id]) {
        next;
    }   
    
    ## filter output on percent identity
    if ($percent_id->{$seq_id}->[$chain_id] < $min_percent_id) {
        next;
    }
    
    ## filter output on percent similarity 
    if ($percent_sim->{$seq_id}->[$chain_id] < $min_percent_sim) {
        next;
    }
    
    ## filter on percent identity
    if (defined($percent_cov) && $percent_cov->{$seq_id}->[$chain_id] < $min_percent_cov) {
        next;
    }
 
    ## start and end positions are flipped for Minus strand
    my ($start, $end) = ($btab[6], $btab[7]);
    if ($start > $end) {
        ($start, $end) = ($end, $start);
    }
    
    ## if we haven't created a match feature yet
    unless ($match_to_gff->{$seq_id}->[$chain_id]) {

        my $match_feature_id = $seq_id.'.'.lc($feature_match).'_'.++$feature_count->{$seq_id}->{'match'};
        $parent_id->{$seq_id}->[$chain_id] = $match_feature_id; 
       
        ## setup feature tags hash
        my $tags_ref = {
                        ID           => $match_feature_id,
                        Name         => $subject_accession,
                        Identity     => sprintf("%.f", $percent_id->{$seq_id}->[$chain_id]),
                        Similarity   => sprintf("%.f", $percent_sim->{$seq_id}->[$chain_id]),
#                               Match    =>
                       }; 
        ## add a Note attribute if the flag was enabled                               
        if ($add_note) {
            $tags_ref->{'Note'} = $note;
        }       

        ## Add coverage if we've done the calculation
        if (defined($percent_cov->{$seq_id}->[$chain_id])) {
            $tags_ref->{'Coverage'} = sprintf("%.f", $percent_cov->{$seq_id}->[$chain_id]);
        }
        
        my $match_feature = Bio::SeqFeature::Generic->new( 
            -seq_id       => $seq_id,
            -start        => $match->{$seq_id}->[$chain_id]->{'match'}->[0], 
            -end          => $match->{$seq_id}->[$chain_id]->{'match'}->[1], 
            -strand       => $match->{$seq_id}->[$chain_id]->{'strand'}, 
            -primary      => $feature_match, # -primary_tag is a synonym
            -source_tag   => $method,
            -display_name => 'Display name',
            -score        => $match->{$seq_id}->[$chain_id]->{'score'},
            -tag          => $tags_ref,
                                                         );
        
        print $outfh $gff->gff_string($match_feature)."\n";

        $match_to_gff->{$seq_id}->[$chain_id] = 1;
    }
    
    my $hsp_feature_id = $parent_id->{$seq_id}->[$chain_id].'.'.lc($feature_match_part).'_'.++$feature_count->{$seq_id}->{'hsp'};
    
    my $hsp_feature = Bio::SeqFeature::Generic->new( 
            -seq_id       => $seq_id,
            -start        => $start, 
            -end          => $end,
            -strand       => $strand{$btab[17]}, 
            -primary      => $feature_match_part, # -primary_tag is a synonym
            -source_tag   => $method,
            -display_name => 'Display name',
            -score        => $hsp_score,
            -tag          => { 
                                ID          => $hsp_feature_id,
                                Parent      => $parent_id->{$seq_id}->[$chain_id],
                                Identity    => sprintf("%.f", $hsp_pct_id),
                                Similarity  => sprintf("%.f", $hsp_pct_sim),
                             } );

    print $outfh $gff->gff_string($hsp_feature)."\n";
} 

## opens a file for reading and returns the filehandle
sub get_infh {
    my ($file) = @_;

    open (my $fh, $file) || confess "Failed opening '$file' for reading: $!";

    return $fh;
}

## returns filehandle to a file for output or STDOUT if no file name is provided
sub get_outfh {
    my ($file) = @_;
   
    if ($file) {
        open (my $fh, ">$file") || confess "Failed opening '$file' for writing: $!";
        return $fh;
    } else {
        return \*STDOUT;
    }
}

## cleans up the query header line from the btab
## and returns a useful accession
sub get_clean_accession {
    my ($header) = @_;
    
    ## clean up query sequence id
    $header =~ /^(\S+)/;
    my $seq_id = $1;
    
    if ($seq_id =~ /gb\|([^\|^\.]+)(\.\d+)?\|/) {
        $seq_id = $1;
    }

    return $seq_id;
}

## do some cleanup on recognized subject header formats
sub get_clean_note {
    my ($note) = @_;

    if ($note =~ /^Cluster: /) {
        $note =~ /^Cluster: (.*?)(; |$)/ || confess "Failed to capture cluster function from uniref header";
        $note = $1;
    }
    
    return $note;
}
