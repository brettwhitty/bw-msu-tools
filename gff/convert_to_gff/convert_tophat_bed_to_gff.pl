#!/usr/bin/perl

##
## Conversion script for converting TopHat bed files into GFF3
## suitable for loading to gbrowse to visualize support for junctions
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;

use Carp;
use Getopt::Long;
use Bio::Tools::GFF;

my ($input, $output, $version, $method, $parent_type, $child_type, $zero_strip);

GetOptions(
    'i|input=s'         =>  \$input,
    'o|output=s'        =>  \$output,
    'v|version=s'       =>  \$version,
    'm|method=s'        =>  \$method,
    'p|parent_type=s'   =>  \$parent_type,
    'c|child_type=s'    =>  \$child_type,
    'z|zero_strip!'     =>  \$zero_strip,
);

## GFF version
$version        ||= 3;
## method defaults to tophat_junctions
$method         ||= 'tophat_junctions';
## parent GFF feature
$parent_type    ||= 'match';
## child GFF feature
$child_type     ||= 'match_part';

my $infh;
if (defined($input)) {
    open $infh, '<', $input || croak "Failed opening '$input' for reading: $!";
} else { 
    $infh = \*STDIN;
}
my $outfh;
if (defined($output)) {
    open $outfh, '>', $output || croak "Failed opening '$output' for writing: $!";
} else { 
    $outfh = \*STDOUT;
}


my $gff = new Bio::Tools::GFF(
                                -gff_version => $version,
                              );

my $juncs = {};

while (<$infh>) {
    chomp;

    ## skip header
    if (/^track name/) {
        next;
    }
    
    my @t = split("\t", $_);
    
    ## NOTE: numbering is interbase in BED files
    my ($start, $end) = ($t[1], $t[2]);

    ##      [------------]...........[--------------]
    ##  start                                       end
    ##         left_size                right_size
    ##      left_offset              right_offset
    ##      0                        0 + left_size + intron_size
    ##
    ## left_start = start
    ## left_end = left_start + left_len
    ## right_start =  start + right_offset
    ## right_end = right_start + right_len
    
    my ($left_len, $right_len) = split(/,/, $t[10]);
    my ($left_offset, $right_offset) = split(/,/, $t[11]);

    my $left_start = $start;
    my $left_end = $left_start + $left_len;

    my $right_start = $start + $right_offset;
    my $right_end = $right_start + $right_len;

    ## create a unique junction id using the coordinates
    ## seq_id.left_end.right_start.strand
    my $junc_id = join('.', ($t[0], $left_end, $right_start, $t[5]));

    ## keep a total of junction depth of coverage
    $juncs->{$junc_id}->{'coverage'} += $t[4];

    ## keep a counter of the number of libraries it's found in
    $juncs->{$junc_id}->{'libs'}++;

    ## set left bound
    $juncs->{$junc_id}->{'left_bound'} = ( 
            ! defined($juncs->{$junc_id}->{'left_bound'})
            || $left_start < $juncs->{$junc_id}->{'left_bound'}
        ) 
        ? $left_start : $juncs->{$junc_id}->{'left_bound'};
    
    ## set right bound
    $juncs->{$junc_id}->{'right_bound'} = ( 
            ! defined($juncs->{$junc_id}->{'right_bound'})
            || $right_end > $juncs->{$junc_id}->{'right_bound'}
        ) 
        ? $right_end : $juncs->{$junc_id}->{'right_bound'};
    
}

foreach my $junc_id(keys(%{$juncs})) {

    my ($seq_id, $left_end, $right_start, $strand) = split(/\./, $junc_id);

    my $junc_atts = $juncs->{$junc_id};

    my $parent_id = $junc_id;


    if ($zero_strip) {
        ## strip zero padding from PGSC IDs
        $parent_id =~ s/([A-Z]+)[0]+/$1/g;
    }
    
    my $left_len = $left_end - $junc_atts->{'left_bound'};
    my $right_len = $junc_atts->{'right_bound'} - $right_start;
    my $intron_len = $right_start - $left_end;

    ## create new SeqFeature for the parent feature
    my $parent_feature = Bio::SeqFeature::Generic->new(
            -seq_id       => $seq_id,
            -source_tag   => $method,
            -primary      => $parent_type, ## -primary_tag is a synonym
            -start        => $junc_atts->{'left_bound'} + 1, ## GFF is 1-based
            -end          => $junc_atts->{'right_bound'},
            -score        => $junc_atts->{'coverage'},
            -strand       => $strand,
            -display_name => 'Display name',
            -tag          => { 
                'ID'        =>  $parent_id,
                'Gap'       =>  "M$left_len D$intron_len M$right_len",
                'libs'      =>  $junc_atts->{'libs'},
            },
    );

    ## create new SeqFeature for left child
    my $child_feature_left = Bio::SeqFeature::Generic->new(
            -seq_id       => $seq_id,
            -source_tag   => $method,
            -primary      => $child_type, ## -primary_tag is a synonym
            -start        => $junc_atts->{'left_bound'} + 1, ## GFF is 1-based
            -end          => $left_end,
            -score        => $junc_atts->{'coverage'},
            -strand       => $strand,
            -display_name => 'Display name',
            -tag          => {
                'ID'        =>  $parent_id.'.0',
                'Parent'    =>  $parent_id,
                'libs'      =>  $junc_atts->{'libs'},
            },
    );

    ## create new SeqFeature for right child
    my $child_feature_right = Bio::SeqFeature::Generic->new(
            -seq_id       => $seq_id,
            -source_tag   => $method,
            -primary      => $child_type, ## -primary_tag is a synonym
            -start        => $right_start + 1, ## GFF is 1-based
            -end          => $junc_atts->{'right_bound'},
            -score        => $junc_atts->{'coverage'},
            -strand       => $strand,
            -display_name => 'Display name',
            -tag          => {
                'ID'        =>  $parent_id.'.1',
                'Parent'    =>  $parent_id,
                'libs'      =>  $junc_atts->{'libs'},
            },
    );

    print $outfh $gff->gff_string($parent_feature)."\n";
    print $outfh $gff->gff_string($child_feature_left)."\n";
    print $outfh $gff->gff_string($child_feature_right)."\n";

}
