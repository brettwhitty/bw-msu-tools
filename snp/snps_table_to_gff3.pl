#!/opt/rocks/bin/perl

## Converts a predicted SNPs table into GFF3

#use strict;
use warnings;

use lib "/home/whitty/SVN/lib";
use MyIO;

use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Getopt::Long;

my ($input, $output, $method, $gff_version);

## defaults
$method = 'predicted_SNP';
$gff_version = '3';

my $result = GetOptions(
                            'method|m=s'            =>  \$method,
                            'input|i=s'             =>  \$input,
                            'output|o=s'            =>  \$output,
                            'version|v=s'           =>  \$gff_version,
                       );

my $infh = get_infh($input);
my $outfh = get_outfh($output);

my $gff = new Bio::Tools::GFF(
                                -gff_version => $gff_version,
                              );

## read the header line from the SNP file
my $header = <$infh>;
chomp $header;
$header =~ s/^# (.*)/$1/ || die "Unexpected header format";
## header line will be used to reference the columns by name
my @col_names = split(/\t/, $header);

my $feat_counter = {};
while (<$infh>) {
    chomp;

    my @cols = split(/\t/, $_);
    
    my $num_cols = scalar(@cols);
    for (my $i = 0; $i < $num_cols; $i++) {
        ${$col_names[$i]} = $cols[$i];
    }

    my $feature_id = ${'put_id'}.'-SNP-'.${'snp_loc'}.'-'.${'ref_base'}.'-'.${'snp_base'};

    my $ref_support = ${
                            {
                                split(/:/, ${'base_freq_a'}),
                                split(/:/, ${'base_freq_c'}),
                                split(/:/, ${'base_freq_g'}),
                                split(/:/, ${'base_freq_t'}),
                            }
                          }{${'ref_base'}};

    my $snp_support = ${
                            {
                                split(/:/, ${'base_freq_a'}),
                                split(/:/, ${'base_freq_c'}),
                                split(/:/, ${'base_freq_g'}),
                                split(/:/, ${'base_freq_t'}),
                            }
                          }{${'snp_base'}};

    my $snp_feature = Bio::SeqFeature::Generic->new(
            -seq_id       => ${'put_id'},
            -start        => ${'snp_loc'},
            -end          => ${'snp_loc'},
            -strand       => '+',
            -primary      => 'SNP', # -primary_tag is a synonym
            -source_tag   => $method,
            -display_name => 'Display name',
            -score        => undef,
            -tag          => {
                                ID              => $feature_id,
                                ref_base        => ${'ref_base'},
                                ref_support     => $ref_support,
                                alt_base        => ${'snp_base'},
                                alt_support     => $snp_support,
                                depth           => ${'snp_cov'},
                                base_support    => [ 
                                                        ${'base_freq_a'},
                                                        ${'base_freq_c'},
                                                        ${'base_freq_g'},
                                                        ${'base_freq_t'},
                                                   ],
                                Dbxref      => 'taxon:'.${'taxon_id'},

                             } );

    print $outfh $gff->gff_string($snp_feature)."\n";
};
