#!/usr/bin/perl

use strict;
use warnings;
use Carp;

use lib "/home/whitty/SVN/lib";
use MyIO;

use Getopt::Long;
use Cwd;
use Datastore::MD5;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $mcl_result;
my $all_fasta;
my $out_dir;

my $result = GetOptions(
                            'fasta_db|f=s'      =>  \$all_fasta,
                            'input_file|i=s'    =>  \$mcl_result,
                            'output_dir|o=s'    =>  \$out_dir,
                       );

unless (-e $all_fasta) {
    confess "Input fasta database '$all_fasta' does not exist";
}
unless (-e $mcl_result) {
    confess "Specified orthoMCL output file '$mcl_result' does not exist";
}

my $db = Bio::DB::Fasta->new($all_fasta);

unless ($out_dir) {
    $out_dir = getcwd()."/clusters/";
}
                       
my $ds = new Datastore::MD5( root => $out_dir, depth => 2 );

my $infh = get_infh($mcl_result);

while (<$infh>) {
    chomp;

    /^(ORTHOMCL\d+)\((\d+) genes,(\d+) taxa\):\s+(.*)$/;
    my ($id, $genes, $taxa, $rest) = ($1, $2, $3, $4);

    my @rest = split(/\s+/, $rest);
    
    if (scalar(@rest) != $genes) {
        confess "Mismatch of split operation contents and gene count in cluster --- this is a problem";
    }
    
#    print "$id\t$genes\t$taxa\n";

    $ds->chdir($id);
    
    my $outfh = get_outfh("$id.fsa");
    my $fasta_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $outfh );

    foreach my $cluster_member( @rest ) {
        $cluster_member =~ /^([^(]+)\((.*)\)$/ || confess "Failed to match cluster member name";
        my ($seq_id, $src) = ($1, $2);

        my $seq_obj = $db->get_Seq_by_id($seq_id);    
        $fasta_out->write_seq($seq_obj);
    }
}
