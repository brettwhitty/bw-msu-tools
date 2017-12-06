#! /usr/bin/perl

# convert_and_decode_phyi_file.pl

# 6 August 2009
# Kevin Childs

# This script decodes the sequence names in a phylip interleaved multiple sequence file.  
# The file is also converted to a multiple sequence alignment fasta file.

use Getopt::Std;
use Bio::AlignIO;
use Bio::SeqIO;

use strict;

my $usage = "\n$0 -i phylip_file -o msa_file -d decoder_file\n\n";

our ( $opt_i, $opt_o, $opt_d, $opt_h );
getopts("i:o:d:h") or die usage();

if ($opt_h) {
    print $usage;
    exit;
}

my $input_file = $opt_i;
my $output_file = $opt_o;
my $decoder_file = $opt_d;

if (   !defined($input_file)
       || !( -e $input_file )
       || !defined($output_file)
       || ( -e $output_file )
       || !defined($decoder_file)
       || !(-e $decoder_file)) {
    die "\nMissing or invalid input values.\n$input_file\n$output_file\n$decoder_file\n" . $usage;
}

my $short_phylip_file_name;
if ($input_file =~ /([\w\.]+)\.phyi$/) {
    $short_phylip_file_name = $1;
}
else {
    die "\nUnable to parse cluster name out of file name.\nX$input_file" . "X\n\n";
}

my %decoder;

open DECODER, "$decoder_file" || die "\nUnable to open decoder file, $decoder_file, for reading.\n\n";
while (my $line = <DECODER>) {
    chomp $line;
    my $cluster_found = 0;
    if ($line =~ /#$short_phylip_file_name$/) {
    # This is our cluster.
    while ($line = <DECODER>) {
        chomp $line;
        my @elems = split "\t", $line;
        if ($elems[0] =~ /\D/) {
        # This is something beside a number.  We're finished.
        $cluster_found = 1;
        last;
        }

        $decoder{$elems[0]} = $elems[1];
    }
    }
    if ($cluster_found == 1) {
    last;
    }
}

my $phylip = Bio::AlignIO->new(-file => $input_file,   -format=>'phylip');
my $msf = Bio::SeqIO->new(-file => ">$output_file",   -format=>'fasta');

# Do the conversion.
while(my $align = $phylip->next_aln) { 
    foreach my $seq ($align->each_seq()) {
    my $coded_id = $seq->display_id();
    $seq->display_id($decoder{$coded_id});
    print "$coded_id\n$decoder{$coded_id}\n";
    my $sequence = $seq->seq();
    print "$sequence\n";
    $msf->write_seq($seq);
    }
}

exit;

