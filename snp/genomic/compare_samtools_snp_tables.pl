#!/opt/rocks/bin/perl

use strict;
use warnings;

use Carp;
use Getopt::Long;

## provided a list of SNP tables output from samtools consensus pileup output
## will generate a table of SNP coordinates with a vector representing the
## presence/absence of the SNPs
##
## Brett Whitty, whitty@msu.edu

my ($snp_file_list, $present, $absent);

GetOptions(
    'snp_files|s=s'     =>  \$snp_file_list,
    'present|p=s'       =>  \$present,
    'absent|a=s'        =>  \$absent,
);

if (! defined($snp_file_list)) {
    confess "Must provide a comma-delimited list of SNP files with --snp_files";
}

## using y/n here for easier grep'ing
$present ||= 'y';
$absent ||= 'n';

## only one character allowed
$present = substr($present, 0, 1);
$absent = substr($absent, 0, 1);

## SNP files are to be provided as a list delimited by ','
my @snp_files = split(/,/, $snp_file_list);

if (scalar(@snp_files) < 2) {
    confess "Can't compare <2 datasets!";
}

#my @snp_files = (
#    '/run/whitty/proposal/snps/final_dataset/snps_for_browser/67.STP_AA.filtered.consensus.snp',
#    '/run/whitty/proposal/snps/final_dataset/snps_for_browser/rnaseq/67_leaf.rnaseq.snp',
#    '/run/whitty/proposal/snps/final_dataset/snps_for_browser/90.STP_AB.filtered.consensus.snp',
#    '/run/whitty/proposal/snps/final_dataset/snps_for_browser/rnaseq/90_leaf.rnaseq.snp',
#);

my $bits = {};

my $file_count = scalar(@snp_files);

for (my $i = 0; $i < $file_count; $i++) {

    open my $infh, '<', $snp_files[$i] or die "failed to open '".$snp_files[$i]."' for reading:$!";

    while (<$infh>) {
        my @t = split("\t", $_);

        ## generate a unique key using the sequence molecule ID, position and base
        my $key = $t[0]."\t".$t[1]."\t".$t[3];

        if (! defined($bits->{$key})) {
            @{$bits->{$key}} = split(//, $absent x $file_count);
        }
        $bits->{$key}[$i] = $present;
    }
}

foreach my $key (sort keys %{$bits}) {
    print $key."\t".join('', @{$bits->{$key}})."\n";
}
