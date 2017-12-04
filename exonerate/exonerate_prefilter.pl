#!/usr/bin/perl

##
## This script can be used to pre-filter an exonerate query database
## and exclude query sequences with no significant sequence similarity to the target
## to reduce the size of the query database for the exonerate run and reduce overall run-time
##
## The BLAST parameters used also support using softmasked target sequence
##
## CAUTION: 
## The first time this script is run, Bio::SeqIO will index the query database.
## Concurrent execution of the script before indexing is complete can cause corruption
## of the index file, and the script won't be able to properly fetch sequences from the
## query database when generating the output database.

use warnings;
use strict;

use Getopt::Long;

use Bio::SearchIO;
use File::Basename;
use Cwd qw/ getcwd abs_path /;
use File::Temp qw/ :POSIX tempfile tempdir /;
use File::Path;

use Bio::DB::Fasta;
use Bio::SeqIO;

my $cpus = 1;
my $cutoff = 0.01; ## this shouldn't be decreased from 0.01 which was found to be lossless in test runs

my $db_file;
my $query_file;
my $type = 'n';
my $output;
my $work_dir = getcwd();

my $result = GetOptions(
                           'target|t=s'     =>  \$db_file,
                           'query|q=s'      =>  \$query_file,
                           'query_type=s'   =>  \$type,
                           'output|o=s'     =>  \$output,
                           'work_dir|w=s'   =>  \$work_dir,
                       );

my $prog = {
            'n' =>  'tblastx',
            'p' =>  'tblastn',
        };

$work_dir .= "/tmp";        

mkpath($work_dir);

my $temp_dir = tempdir( DIR => $work_dir, CLEANUP => 1);                    
my (undef, $tmp_out) = tempfile(DIR => $temp_dir, OPEN => 0);

$tmp_out = abs_path($tmp_out);

$query_file = abs_path($query_file);

system("blastall -p $prog->{$type} -d $db_file -i $query_file -U T -F \"m S\" -a $cpus -e $cutoff -o $tmp_out");

    my $fasta_db = new Bio::DB::Fasta($query_file);
    
    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -file   =>  $tmp_out,
                                        );
    
                                    
    my ($out_fh, $out_name) = tempfile( DIR => $work_dir, CLEANUP => 0 );
    
    my $fasta_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $out_fh );
    
    while( my $result = $blast_report->next_result ) {

        my $query_name = $result->query_name;

        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;

            my $seq_obj = $fasta_db->get_Seq_by_id($query_name);

            $fasta_out->write_seq($seq_obj);
            
        }

    }   

print $out_name."\n";
unlink($tmp_out); 
