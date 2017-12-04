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
                            'input_file|i=s'    =>  \$mcl_result,
                            'output_dir|o=s'    =>  \$out_dir,
                       );

unless ($mcl_result) {
    confess "Must specify an orthomcl mcl output file with --input_file";
}

unless ($out_dir) {
    confess "Must specify a Datastore::MD5 root diretory containing the clusters";
}
                       
unless (-e $mcl_result) {
    confess "Specified orthoMCL output file '$mcl_result' does not exist";
}
                       
my $ds = new Datastore::MD5( root => $out_dir, depth => 2 );

my $infh = get_infh($mcl_result);

while (<$infh>) {
    chomp;

    /^(ORTHOMCL\d+)\((\d+) genes,(\d+) taxa\):\s+(.*)$/;
    my ($id, $genes, $taxa, $rest) = ($1, $2, $3, $4);

    print STDERR "Entering '$id' directory '".$ds->id_to_dir($id)."'...\n";
    
    if (-e $ds->id_to_dir($id)."/$id.uniref50.blastp.raw") {
        print STDERR "Output file exists for '$id', skipping BLAST run...\n";
    } else {
        $ds->system($id, "nice -10 blastp /home/whitty/projects/db/UniRef/uniref50.fasta $id.fsa -cpus 8 -e 1e-5 -b 1 -v 1 >$id.uniref50.blastp.raw");
    }
    
    $ds->system($id, "/home/whitty/SVN/uniref50_blast_output_to_table.pl $id.uniref50.blastp.raw >$id.uniref50.blastp.txt");
}
