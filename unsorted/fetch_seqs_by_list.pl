#!/usr/bin/perl

##
## Fetch a set of sequences by identifier using a list
##

use strict;
use warnings;

use Bio::DB::Fasta;
use Getopt::Long;

use lib '/home/whitty/SVN/lib';
use MyIO;

my ($db_file, $input, $output);

GetOptions(
    'db|d=s'        =>  \$db_file,
    'input|i=s'     =>  \$input,
    'output|o=s'    =>  \$output,
);

unless (-f $db_file) {
    die "Please provide a path to a fasta file";
}

my $db = Bio::DB::Fasta->new($db_file, -reindex => 1);

my $outfh = get_outfh($output);

my $infh = get_infh($input);
while (<$infh>) {
    chomp;

    my $id = $_;

    my $header = ">".$db->header($id)."\n";
    my $seq = $db->seq($id);
    $seq =~ s/(.{1,60})/$1\n/g;
    print $outfh $header.$seq;
}
