#!/opt/rocks/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;

my $db_list = shift @ARGV;

open (my $infh, '<', $db_list) || die "Provide tab delimited list of database name and fasta file: $!";

while (<$infh>) {
    chomp;

    my ($name, $fasta) = split("\t", $_);

    my $db = new Bio::DB::Fasta($fasta);

    my @ids = $db->ids();

    print join(':', (
            $name,
            join(' ', @ids),
    ))."\n";
}
