#!/usr/bin/perl

use strict;
use warnings;

use lib '/home/whitty/SVN/lib';
use NCBITaxonomyTree;
use Net::FTP;
use File::Temp qw{ :POSIX };
use Date::Calc qw{ Decode_Month };

my $tree = new NCBITaxonomyTree;

my $ftp_root = 'ftp://ftp.ncbi.nih.gov/repository/UniGene/';
my $host = 'ftp.ncbi.nih.gov';
my $path = '/repository/UniGene';

my $ftp = new Net::FTP($host, Debug => 0);
$ftp->login('anonymous', 'me@here.edu');
$ftp->cwd($path);
my @files = $ftp->ls();

foreach my $file(@files) {
    my $name = $file;
    $name =~ s/_/ /;

    if ($tree->valid_scientific_name($name) && $tree->name_has_ancestor_name($name, 'Solanaceae')) {
        my @info = $ftp->ls($file."/*.info");

        my $taxon_id = $tree->get_taxon_id($name);
        my $url = $ftp_root.$file.'/';

        my $tmpfile = tmpnam();
        $ftp->get($info[0], $tmpfile);

        open my $infh, '<', $tmpfile || die "Failed to open fetched file '$tmpfile': $!";

        my ($version, $date, $seq_count, $cluster_count);
        while (<$infh>) {
            chomp;

            if (/^UniGene Build \D+(\d+)/) {
                $version = $1;
            }
            if (/(\d+)\s+total sequences in clusters\s*$/) {
                $seq_count = $1;
            }
            if (/(\d+)\s+sets total\s*$/) {
                $cluster_count = $1;
            }
            if (/ESTs are from dbEST through (\d+)\s+(\S+)\s+(\d+)/) {
                my ($d, $m, $y) = ($1, $2, $3);

                $m = sprintf("%02d", Decode_Month($m));

                $date = "$y-$m-$d";
            }
        }
        print join("\t", (
                $taxon_id,
                'NCBI UniGene',
                $version,
                $date,
                $seq_count,
                $cluster_count,
                $url,
        ))."\n";
    }
}
