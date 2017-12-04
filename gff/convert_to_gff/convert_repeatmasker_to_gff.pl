#!/usr/bin/perl

## use bioperl modules to convert repeatmasker output to GFF
## for eventual loading to a gbrowse database
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;
use Carp;
use Getopt::Long;

use Bio::Tools::RepeatMasker;
use Bio::Tools::GFF;

use Data::Dumper;

my $input;
my $output;

my $result = GetOptions(
                            'input|i=s'     =>  \$input,
                            'output|o=s'    =>  \$output,
                       );

unless ($input && -f $input) {
    confess "Must provide correct path to an input file with --input flag";
}    

my $parser = new Bio::Tools::RepeatMasker(-file => $input);

my $gff = new Bio::Tools::GFF(-gff_version => 3);

my $counter = 0;
while (my $result = $parser->next_result()) {
   
    my $query_feature = $result->feature1;
    
    ## clean up the genbank accession if it didn't happen before the repeatmasker run
    if ($query_feature->seq_id =~ /^gi\|\d+\|gb\|([^\.^\|]+)/) {
        $query_feature->seq_id($1);
    }
    my $id = $query_feature->seq_id;
    ## strip leading zeroes from PGSC IDs for a shorter feature ID
    if ($id =~ /^PGSC/) {
        $id =~ s/([A-Z]+)[0]+/$1/g;
    }
    $query_feature->add_tag_value('ID', $id.'.repeat.'.++$counter);
    $query_feature->add_tag_value('Name', $query_feature->primary_tag);
    $query_feature->primary_tag('repeat_region'); ## correct SO type
    
    my $match_feature = $result->feature2;
    
    print $gff->gff_string($query_feature)."\n";
}

