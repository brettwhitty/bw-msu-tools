#!/usr/bin/perl

$| = 1;

use strict;
use warnings;

use Cwd qw(abs_path);
use Bio::SeqIO;

my $infile = shift @ARGV || die "Must specify and input genbank file";

$infile = abs_path($infile);

$infile =~ /^(.*)\.(gbk|gbff)$/ || die "Input file must be a genbank file with the file extension '.gbk' or '.gbff'";

my $outfile = $1.".fsa";

#$in  = Bio::SeqIO->new('-fh' => \*STDIN , '-format' => 'genbank');
#$out = Bio::SeqIO->new('-fh' => \*STDOUT, '-format' => 'fasta');

open(my $infh, $infile) || die "Failed to open input file '$infile' for reading: $!";
open(my $outfh, ">$outfile") || die "Failed to open output file '$outfile' for writing: $!";

my $in  = Bio::SeqIO->new('-fh' => $infh,  '-format' => 'genbank');
#my $out = Bio::SeqIO->new('-fh' => $outfh, '-format' => 'fasta');

while ( my $seq = $in->next_seq() ) {
    my $head = '>gi|'.$seq->primary_id().'|gb|'.$seq->id().'.'.$seq->version().'| '.$seq->desc()
."\n";
    my $seq  = $seq->seq();
    $seq =~ s/(\S{1,60})/$1\n/g;

    print $outfh $head.$seq;
}
