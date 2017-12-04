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

use Bio::Tools::tRNAscanSE;
use Bio::Tools::GFF;

my $input;
my $output;
my $gff_version = 3;

my $result = GetOptions(
                            'input|i=s'         =>  \$input,
                            'output|o=s'        =>  \$output,
                            'gff_version|v=s'   =>  \$gff_version,
                       );

unless ($input && -f $input) {
    confess "Must provide correct path to an input file with --input flag";
}    

my $parser = new Bio::Tools::tRNAscanSE(-file => $input);

my $gff = new Bio::Tools::GFF(-gff_version => $gff_version);

my $counter = 0;
while (my $result = $parser->next_prediction()) {
   
#    $result->add_tag_value('Name', $query_feature->primary_tag);    

    $result->primary_tag('tRNA'); ## SOFA term

    my $gff_string = $gff->gff_string($result);

    my @t = split("\t", $gff_string);

    $t[8] =~ /ID=([^;]+);aminoAcid=([^;]+);codon=(.*)/ || confess "Failed to match attributes";
    
    my ($name, $aa, $codon) = ($1, $2, $3);

    $name =~ tr/:/-/;
    $aa =~ s/-\d+$//;
    
    $t[8] = "ID=".$t[0]."-trna".++$counter.";Name=tRNA-$aa;Note=$codon";
    
    print join("\t", @t)."\n";
}

