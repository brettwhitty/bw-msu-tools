#!/usr/bin/perl

use strict;
use warnings;
use Carp;

use Bio::Tools::GFF;
use Getopt::Long;

use lib "/home/whitty/SVN/lib";
use MyIO;

my $input;
my $output;
my $version;

my $result = GetOptions(
                           'input|i=s'      =>  \$input,
                           'output|o=s'     =>  \$output, 
                           'version=s'      =>  \$version,
                       );

my $infh = get_infh($input);
my $outfh = get_outfh($output);

my $in_gff = new Bio::Tools::GFF(
                        -fh =>  $infh,
                                );

my $out_gff = new Bio::Tools::GFF(
                        -fh             =>  $outfh,
                        -gff_version    => $version,
                                 );

while (my $feature = $in_gff->next_feature) {
    print $outfh $out_gff->gff_string($feature)."\n";
}
