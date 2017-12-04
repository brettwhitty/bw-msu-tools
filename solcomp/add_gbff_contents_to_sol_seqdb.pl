#!/usr/bin/env perl

## Take a GenBank flat file containing one or more records
## and split into individual files that correspond to GI numbers
## that have already been loaded to a Sol::SeqDB database
## and update the database information for those records

use lib '/home/whitty/work/testseqdb/';
#use lib '/home/whitty/SVN/lib';

use strict;
use warnings;
use Carp;

use File::Path;
#use DateTime;
#use Roman;
use Time::Piece;

use Bio::SeqIO;

## my modules
#use NCBITaxonomyTree;
use Sol::SeqDB;      

#my $current_date = DateTime->now()->ymd('');
my $lt = localtime();
my $current_date = $lt->ymd('');

my $months = {
                'JAN' => '01',
                'FEB' => '02',
                'MAR' => '03',
                'APR' => '04',
                'MAY' => '05',
                'JUN' => '06',
                'JUL' => '07',
                'AUG' => '08',
                'SEP' => '09',
                'OCT' => '10',
                'NOV' => 11,
                'DEC' => 12
             };

my $input_file = shift @ARGV || die "Must provide a GenBank flat file";

#my $tree = new NCBITaxonomyTree(                               
#                                'path'  =>  '/projects/whitty_home/db/NCBI/Taxonomy_dump',
#                               ); 

my $db = new Sol::SeqDB();

my $outfile;
my $outfh;
my $acc;
my $gi;
my $acc_ver;
my $org;
my $taxon_id;
my $def;
my $record_date;
my $division;
my $mol_type;
my $topology;
my $seq_len;
my $chrom = 'unknown';

my %stats = ();


my $gbk_out_string = '';

open (my $in_fh, $input_file) || croak "Failed to open '$input_file' for reading: $!";

my $count = 0;
while (<$in_fh>) {
    if (/^$/) { next; }
   
    $gbk_out_string .= $_;
    
    ## match the LOCUS line and parse out accession
    if (/^LOCUS\s+(\S+)\s+(\d+)\sbp\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
        ($acc, $seq_len, $mol_type, $topology, $division, $record_date) = ($1, $2, $3, $4, $5, $6);   

        my ($d, $m, $y) = split("-", $record_date);
        $m = $months->{$m};
        $record_date = "$y-$m-$d";

        ## don't allow corrupted records to slip through
        if (++$count > 1) {
            croak "Fatal: Encountered a second LOCUS line without hitting end of record";
        }
    } 
    
    if (/^VERSION\s+(\S+)\s+GI:(\d+)/) {
        ($acc_ver, $gi) = ($1, $2);
    }

    ## match the end of the record
    if (/^\/\/$/) {
        unless ($gi) {
            confess "Got to end of record without parsing gi";
        }
    
        my $out_gbk = "$gi.gbk";

        ## write the genbank flat file to the datastore
        my $gbk_fh = $db->{'file'}->create_fh($out_gbk);
        print $gbk_fh $gbk_out_string; 
        $gbk_out_string = '';
        close $gbk_fh;

        ## add the record to the database
        my $gb_atts = $db->get_gb_atts($gi) or confess "Failed to get atts for GI '$gi'";
        $db->add_gi($gi, $gb_atts);

        ## reset all variables
        $count = 0;
        $acc = '';
        $gi = '';
        $acc_ver = '';
        $org = '';
        $taxon_id = '';
        $def = '';
        $record_date = '';
        $division = '';
        $mol_type = '';
        $topology = '';
        $seq_len = '';
        $chrom = 'unknown';
       
    }
}
## update taxon table if there's new taxa added maybe
$db->add_taxon();
## set obsolete anything that's not actually sol
$db->set_gb_obsolete_not_sol();
## set obsolete any old BACs if we've loaded new versions
$db->set_gb_obsolete_not_current();
