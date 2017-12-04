#!/usr/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use Getopt::Long;

use Data::Dumper;
use Bio::SearchIO;
use File::Basename;
use Cwd;
use File::Temp qw/ :POSIX tempfile tempdir /;
use File::Path;

use Bio::DB::Fasta;
use Bio::SeqIO;

my $cpus = 4;
my $cutoff = 0.01; ## this shouldn't be changed from 0.01 which was found to be lossless in test runs

my $db_file;
my $query_file;
my $type = 'n';
my $output;
my $work_dir = getcwd();

my $result = GetOptions(
                           'target|t=s'     =>  \$db_file,
                           'query|q=s'      =>  \$query_file,
                           'query_type=s'         =>  \$type,
                           'output|o=s'     =>  \$output,
                           'work_dir|w=s'   =>  \$work_dir,
                       );

my $prog = {
            'n' =>  'tblastx',
            'p' =>  'tblastn',
        };

$work_dir .= "/tmp";        

mkpath($work_dir) || die "Failed to make work dir '$work_dir': $!";

my $temp_dir = tempdir( DIR => $work_dir, CLEANUP => 1);                    
my (undef, $tmp_out) = tempfile(DIR => $temp_dir, OPEN => 0);
        
system("blastall -p $prog->{$type} -d $db_file -i $query_file -F \"m S\" -a $cpus -e 0.01 -o $tmp_out");

    my $fasta_db = new Bio::DB::Fasta($query_file);
    
    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -file   =>  $tmp_out,
                                        );
    
                                    
    my ($out_fh, $out_name) = tempfile( DIR => $work_dir, CLEANUP => 0 );
    
    my $fasta_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $out_fh );
    
    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
#        my $hit_count = 0;

        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
#            my $coverage = sprintf("%.2f", $hit->frac_aligned_query());
#            my $percent_coverage = $coverage * 100;

#            my $total_length = 0;
#            my $sum = 0;
#            while( my $hsp = $hit->next_hsp ) {
#                ## calculate average % identity across all hsps
#                $total_length += $hsp->length('query');
#                $sum += $hsp->frac_identical * $hsp->length('query');
#            }
#            my $identity = sprintf("%.2f", $sum/$total_length);
#            my $percent_identity = $identity * 100;
        
#            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff) {
                
#                $hit_count++;
                
#                print join("\t", (
#                                    $query_name,
#                                    $hit_name,
#                                    $rank,
#                                    $score,
#                                    $sig,
#                                    $identity,
#                                    $coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
#                                 ))."\n";

            my $seq_obj = $fasta_db->get_Seq_by_id($query_name);

            $fasta_out->write_seq($seq_obj);
            
#            }
            
        }

        #   if ($hit_count) {
            #print STDOUT $query_name."\n"; 
        #   } else {
            #print STDERR $query_name."\n";
        #       print $query_name."\t-\n";
        #   }
        
    }   

print $out_name."\n";
unlink($tmp_out); 
