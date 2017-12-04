#!/usr/bin/perl

## 
## Quick script to generate a pseudomolecule given a fasta database
## and an ordered list of identifiers --- will link with a 100bp
## string of Nn's
##
## Brett Whitty
## whitty@msu.edu
##

use strict;
use warnings;

use Getopt::Long;
use Bio::DB::Fasta;
use Carp;

my ($db_fasta, $id_list, $out_prefix, $defline);

## the linker for joining the seqs
my $linker = 'Nn' x 50;
## can change the linker, but leave this calculation alone
my $linker_len = length($linker);

GetOptions(
    'db|d=s'        =>  \$db_fasta,
    'input|i=s'     =>  \$id_list,
    'output|o=s'    =>  \$out_prefix,
    'defline|n=s'   =>  \$defline,
);

unless (-f $db_fasta) {
    confess "Provide a fasta db with --db flag";
}
unless (-f $id_list) {
    confess "Provide identifier list for scaffolding order with --input flag";
}
unless (defined($out_prefix)) {
    confess "Must provide an output prefix with --output flag";
}
unless (defined($defline)) {
    confess "Must provide definition line for output pseudomolecule with --defline flag";
}

#my $db = new Bio::DB::Fasta($db_fasta, -reindex => 1);
my $db = new Bio::DB::Fasta($db_fasta);

open (my $infh, '<', $id_list) || confess "Failed to open '$id_list' for reading: $!";

my $out_fasta = $out_prefix.'.fsa';
my $out_offset = $out_prefix.'.offsets';

open (my $outfh_fsa, '>', $out_fasta) || confess "Failed to open '$out_fasta' for writing: $!";
open (my $outfh_off, '>', $out_offset) || confess "Failed to open '$out_offset' for writing: $!";

print $outfh_fsa ">$defline\n";

my @ids = ();
while (<$infh>) {
    chomp;
    push (@ids, $_);
}
my $id_count = scalar(@ids);

my $offset = 0;
my $count = 0;
foreach (my $i = 0; $i < $id_count; $i++) {
    
    ## add a linker
    if ($i > 0) {
        $offset = write_char($outfh_fsa, $offset, \$linker);
    }
    
    ## write line to offset file
    print $outfh_off join("\t", (
            $ids[$i],
            $offset
    ))."\n";
    
    ## fetch sequence
    my $seq_obj = $db->get_Seq_by_id($ids[$i]);
    
    ## write sequence to file
    my $seq = $seq_obj->seq();
    $offset = write_char($outfh_fsa, $offset, \$seq);

}
## add a newline to last line if not already done
if ($offset % 60 != 0) {
    print $outfh_fsa "\n";
}

## write string to output file character by character
## add a line feed after every X characters determined
## by '$wrap' --- returns an updated offset value
sub write_char {
    my ($fh, $offset, $string_ref) = @_;

    ## line wrap width
    my $wrap = 60;

    foreach my $char (split("", ${$string_ref})) {
        print $fh $char;
        $offset++;
        if ($offset % $wrap == 0) {
            print $outfh_fsa "\n";
        }
    }

    return $offset;
}
