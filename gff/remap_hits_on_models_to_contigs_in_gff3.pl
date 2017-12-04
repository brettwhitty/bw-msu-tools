#!/opt/rocks/bin/perl

$| = 1;

use strict;
use warnings;
use Carp;

use Getopt::Long;

use Bio::Tools::GFF;

use Data::Dumper;

my $genome_gff   = '';
my $search_gff   = '';
my $gff_version  = 3;
my $dna;
my $add_note;

my $result = GetOptions(
                            'genome_gff|g=s'    =>  \$genome_gff, ## contains coordinates of CDS features for models
                            'gff_version|v=s'   =>  \$gff_version,
                            'search_gff|s=s'    =>  \$search_gff, ## contains coordinates of model-level search alignments
                            'dna!'              =>  \$dna,
                            'add_note!'         =>  \$add_note,
                       );

unless (-f $genome_gff) {
    confess "Must provide genome sequence encoding GFF file that defines positions of models with --genome_gff flag";
}

unless (-f $search_gff) {
    confess "Must provide search result encoding GFF file with --search_gff flag";
}

my $feature;
my $gff;

## keep track of the models in the gff
my $search_models = {};
## keep track of what strand the models are on
my $strand = {};
## map models to their parent genomic sequence feature
my $genome_models = {};

## first do a scan through the search results gff to identify
## features that we need to pull the locations of from the
## genome gff file
$gff = new Bio::Tools::GFF(
                                -file => $search_gff,
                                -gff_version => $gff_version,
                          );

print STDERR "Processing search GFF '$genome_gff' to find model IDs\n";

while ($feature = $gff->next_feature()) {
    $search_models->{$feature->seq_id} = [];
}
$gff->close();

## now scan through the genome gff file containing the models
## and store the CDS feature coordinates for models that have
## search results
$gff = new Bio::Tools::GFF(
                                -file => $genome_gff,
                                -gff_version => $gff_version,
                          );

print STDERR "Processing genome GFF '$genome_gff'\n";

my $feat_counter = {};
my $counter = 0;
## keep track of the models in the gff
while ($feature = $gff->next_feature()) {
    ## only look at features with Parents
    unless ($feature->has_tag('Parent')) { next; }


    my ($id) = ($feature->has_tag('ID')) ? $feature->each_tag_value('ID') : ($feature->primary_tag.++$feat_counter->{$feature->primary_tag});

    my @parents = $feature->each_tag_value('Parent');

    ## some CDSs can have multiple parents (different mRNA isoforms)
    foreach my $parent(@parents) {

        ## we're looking for CDSs with parents that are in the search results gff
        if ($feature->primary_tag eq 'CDS' && defined($search_models->{$parent})) {
            $genome_models->{$parent} = $feature->seq_id;

            unless (defined($strand->{$parent})) {
                $strand->{$parent} = $feature->strand;
            }

            push(@{$search_models->{$parent}}, [$feature->location->start, $feature->location->end]);
        }

    }
}
$gff->close();

## finally scan through the search results gff and remap the feature coordinates
$gff = new Bio::Tools::GFF(
                                -file => $search_gff,
                                -gff_version => $gff_version,
                          );

print STDERR "Remapping search GFF '$search_gff'\n";

## keep track of the models in the gff
while ($feature = $gff->next_feature()) {
    if (scalar @{$search_models->{$feature->seq_id}} < 1) {
        #confess "No coordinates were found in genome gff file '$genome_gff' for '".$feature->seq_id."'";
        carp "No coordinates were found in genome gff file '$genome_gff' for '".$feature->seq_id."', skipping...";
        next;
    }

    my $target_value;
    my $target_id;
    if ($feature->has_tag('Target')) {
        ($target_value) = $feature->each_tag_value('Target');
        $target_value = clean_target_value($target_value);
        $target_value =~ /^(\S+)/;
        $target_id = $1;
        if ($add_note) {
            $feature->add_tag_value('Note', split(/\s+/, $target_value));
        }
        unless ($feature->has_tag('Name')) {
            $feature->add_tag_value('Name', $target_id);
        }
    }

    my $id = '';
## commenting out the code that reassigns IDs
##
#    if ($feature->has_tag('ID')) {
#        ($id) = $feature->each_tag_value('ID');
#        $feature->remove_tag('ID');
#    }
##
##
    my $hsp_id;
    if ($target_id) {
        $hsp_id = $genome_models->{$feature->seq_id}."-".$target_id."-".get_hsp_id($feature->primary_tag);
    } else {
        $hsp_id = $genome_models->{$feature->seq_id}."-".get_hsp_id($feature->primary_tag);
    }
##    $feature->add_tag_value('ID', $hsp_id);

    my $parent = '';
    if ($feature->has_tag('Parent')) {
        ($parent) = $feature->each_tag_value('Parent');
    }

    ## add strand information to feature
    if (defined($strand->{$feature->seq_id})) {
        $feature->strand($strand->{$feature->seq_id});
    }

    my $match_ref = [$feature->location->start, $feature->location->end];
    my $cds_arr_ref = $search_models->{$feature->seq_id};

    my @remapped_match;
    if ($strand->{$feature->seq_id} == 1) {
        @remapped_match = @{remap_forward_interval($match_ref, $cds_arr_ref)};
    } elsif ($strand->{$feature->seq_id} == -1) {
        @remapped_match = @{remap_reverse_interval($match_ref, $cds_arr_ref)};
    } else {
        confess "Unexpected strand information for '".$feature->seq_id."'";
    }

    $feature->seq_id($genome_models->{$feature->seq_id});

    if ($feature->primary_tag eq 'match') {

        ## this should really be the minimum start position
        ## and the maximum end position
        $feature->location->start($remapped_match[0]->[0]);
        $feature->location->end($remapped_match[$#remapped_match]->[1]);
#        $feature->add_tag_value('Target', $target_value);

        print $gff->gff_string($feature)."\n";
    } else {
        foreach my $pos_ref(@remapped_match) {

            $feature->location->start($pos_ref->[0]);
            $feature->location->end($pos_ref->[1]);
#            $feature->add_tag_value('Target', $target_value);

            print $gff->gff_string($feature)."\n";
        }
    }
}
$gff->close();


## remaps the match features and HSPs so they are relative to the assembly and not the model
sub remap_forward_interval {
    ## match_ref will be an array with start and end position for a match or match_part feature
    ## (that needs to be remapped)
    ## cds_arr_ref will be an array of arrays with start and end positions for each CDS interval
    my ($match_ref, $cds_arr_ref) = @_;

    ## sort the CDS segments by start from low to high
    @{$cds_arr_ref} = sort { $a->[0] <=> $b->[0] } @{$cds_arr_ref};

    my $remapped_arr_ref = [];

    my ($match_start, $match_end) = @{$match_ref};

    unless ($dna) {
        $match_start = 1 + ($match_start - 1) * 3;
        $match_end = 1 + ($match_end - 1) * 3 + 2;
    }

    my $cds_ref;
    foreach my $cds_ref(@{$cds_arr_ref}) {
        my ($cds_start, $cds_end) = @{$cds_ref};

        ## the offset we need to add to the match feature
        my $offset = $cds_start - 1;

        ## calculate the mapped start position of the hit
        my $new_start = $match_start + $offset;

        ## if the new start is outside the bounds of this cds, go to the next cds
        if ($new_start > $cds_end) {
            next;
        }

        my $new_end = $offset + $match_end;

        my $overhang_len = $new_end - $cds_end;

        if ($new_end > $cds_end) {
            $new_end = $cds_end;
        }

        push (@{$remapped_arr_ref}, [$new_start, $new_end]);

        if ($overhang_len <= 0) {
            last;
        }
        $match_start = 1;
        $match_end = $overhang_len;
    }

    return $remapped_arr_ref;
}

## remaps the match features and HSPs so they are relative to the assembly and not the model
sub remap_reverse_interval {
    ## match_ref will be an array with start and end position for a match or match_part feature
    ## (that needs to be remapped)
    ## cds_arr_ref will be an array of arrays with start and end positions for each CDS interval
    my ($match_ref, $cds_arr_ref) = @_;

    ## sort the CDS segments by start from low to high
    @{$cds_arr_ref} = sort { $b->[0] <=> $a->[0] } @{$cds_arr_ref};

    my $remapped_arr_ref = [];

    my ($match_start, $match_end) = @{$match_ref};

    unless ($dna) {
        $match_start = 1 + ($match_start - 1) * 3;
        $match_end = 1 + ($match_end - 1) * 3 + 2;
    }

    my $cds_ref;
    foreach my $cds_ref(@{$cds_arr_ref}) {
        my ($cds_start, $cds_end) = @{$cds_ref};

        ## calculate the mapped start position of the hit
        my $new_start = $cds_end - ($match_start - 1); #verified

        ## if the new start is outside the bounds of this cds, go to the next cds
        if ($new_start < $cds_start) {
            next;
        }

        my $new_end = $cds_end - ($match_end - 1); # verified

        my $overhang_len = $cds_start - $new_end;

        if ($new_end < $cds_start) {
            $new_end = $cds_start;
        }

        push (@{$remapped_arr_ref}, [$new_end, $new_start]);

        if ($overhang_len <= 0) {
            last;
        }
        $match_start = 1;
        $match_end = $overhang_len;
    }

    @{$remapped_arr_ref} = sort { $a->[0] <=> $b->[0] } @{$remapped_arr_ref};

    return $remapped_arr_ref;
}


sub clean_target_value {
    my ($target) = @_;

    if ($target =~ /^Sequence:(.*)/) {
        $target = $1;
    }

    return $target;
}

{
my $hsp_counter = 0;

    sub get_hsp_id {
        my ($type) = @_;
        return "$type-".$hsp_counter++;
    }

}
