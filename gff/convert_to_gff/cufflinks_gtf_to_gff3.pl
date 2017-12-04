#!/usr/bin/perl

use strict;
use warnings;

## Quick script to convert cufflinks GTF to GFF3
##
## Brett Whitty
## whitty@msu.edu

use lib "/home/whitty/SVN/lib";
use MyIO;

use Bio::Tools::GFF;
use Getopt::Long;
use Carp;

my ($input, $output, $method, $gff_version);

## defaults
$gff_version = '3';

my $result = GetOptions(
                            'method|m=s'            =>  \$method,
                            'input|i=s'             =>  \$input,
                            'output|o=s'            =>  \$output,
                            'version|v=s'           =>  \$gff_version,
);

my $gff = new Bio::Tools::GFF(
                                -gff_version => $gff_version,
);

my $infh = get_infh($input);
my $outfh = get_outfh($output);

while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    ## collect attributes for making a SeqFeature::Generic object
    my $feat_atts = {
        'seq_id'        =>  $t[0],
        'source_tag'    =>  $method || $t[1],   ## method
        'primary'       =>  $t[2],              ## feature type
        'start'         =>  $t[3],
        'end'           =>  $t[4],
        'score'         =>  $t[5],
        'strand'        =>  $t[6],
    };

    ## parse the attributes column of the GTF file
    my @attribs = split(/;/, $t[8]);

    my $atts = {};
    foreach my $attrib(@attribs) {
        $attrib =~ s/^\s+|\s+$//g;
        $attrib =~ /^(\S+) "([^"]+)"/ || croak "Failed to match attribute string: $attrib";
        my ($att, $value) = ($1, $2);
        $atts->{$att} = $value;
    }

    ## create appropriate IDs for the feature types
    if ($feat_atts->{'primary'} eq 'transcript') {
        $atts->{'ID'} = $atts->{'gene_id'};
    } elsif ($feat_atts->{'primary'} eq 'exon') {
        $atts->{'ID'} = $atts->{'gene_id'}.'.'.$atts->{'exon_number'};
        $atts->{'Parent'} = $atts->{'gene_id'};
    } else {
        croak "Unexpected feature type '".$feat_atts->{'primary'}."' in GTF";
    }

    ## create new SeqFeature
    my $feature = Bio::SeqFeature::Generic->new(
            -seq_id       => $feat_atts->{'seq_id'},
            -source_tag   => $feat_atts->{'source_tag'},
            -primary      => $feat_atts->{'primary'}, ## -primary_tag is a synonym
            -start        => $feat_atts->{'start'},
            -end          => $feat_atts->{'end'},
            -score        => $feat_atts->{'score'},
            -strand       => $feat_atts->{'strand'},
            -display_name => 'Display name',
            -tag          => $atts,
    );

    ## write GFF line
    print $outfh $gff->gff_string($feature)."\n";

}
