#!/usr/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::DB::Fasta;

use File::Basename;
use Getopt::Long;

use Data::Dumper;


my ($query_file, $db_file, $identity_cutoff, $similarity_cutoff, $coverage_cutoff, 
    $evalue_cutoff, $b, $v, $cpu, $program, $len_cutoff);

$b = 500;
$v = 500;
$cpu = 4;

my $result = GetOptions(
                'query|q=s'     =>  \$query_file,
                'db|d=s'        =>  \$db_file,
                'identity|i=s'  =>  \$identity_cutoff,
                'similarity|s=s'  =>  \$similarity_cutoff,
                'coverage|c=s'  =>  \$coverage_cutoff,
                'evalue|e=s'    =>  \$evalue_cutoff,
                'cpu=i'         =>  \$cpu,
                'b=i'           =>  \$b,
                'v=i'           =>  \$v,
                'len|l=i'       =>  \$len_cutoff,
                       );

unless ($identity_cutoff) { $identity_cutoff = 30;}
unless ($similarity_cutoff) { $similarity_cutoff = 30;}
unless ($coverage_cutoff) { $coverage_cutoff = 20;}
unless ($evalue_cutoff) { $evalue_cutoff = '1e-5';}
unless ($program) { $program = 'blastp'; }
unless ($len_cutoff) { $len_cutoff = 75; }

my $factory= Bio::Tools::Run::StandAloneBlast->new(
                                                    'program'   =>  $program,
                                                    'database'  =>  "$db_file",
                                                    'e'         =>  $evalue_cutoff,
                                                    'a'         =>  $cpu,  # 4 processors
                                                    'b'         =>  $b,
                                                    'v'         =>  $v,
                                                    'F'         =>  '"m S"',
                                                    's'         =>  'T',
                                                  ); 

my $fasta_file = Bio::SeqIO->new(-file => $query_file, -format => 'Fasta', );

while (my $seq = $fasta_file->next_seq()) {

    ## len cutoff was 40 here
    if ($seq->length <= $len_cutoff) { next; }

    my $blast_report = $factory->blastall($seq);
 
    process_blast($blast_report);
    
}

sub process_blast {
    my ($blast_report) = @_;

    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
        my $hit_count = 0;
        
        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            
            ## skip self hits
            if ($hit_name eq $query_name) { next; }
            
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
            my $coverage = sprintf("%.3f", $hit->frac_aligned_query());
            my $percent_coverage = $coverage * 100;

            my $total_length = 0;
            my $id_sum = 0;
            my $sim_sum = 0;
            while( my $hsp = $hit->next_hsp ) {
                ## calculate average % identity across all hsps
                $total_length += $hsp->length('query');
                $id_sum += $hsp->frac_identical * $hsp->length('query');
                $sim_sum += $hsp->frac_conserved * $hsp->length('query');
            }
            my $identity = sprintf("%.3f", $id_sum/$total_length);
            my $similarity = sprintf("%.3f", $sim_sum/$total_length);
            my $percent_identity    = $identity * 100;
            my $percent_similarity  = $similarity * 100;
        
            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff && $percent_similarity >= $similarity_cutoff && $total_length >= $len_cutoff) {
                
                $hit_count++;
              
                print join("\t",
                                 (
                                    $query_name,
                                    $hit_name,
                                    #    $rank,
                                    $score,
                                    $sig,
                                    #$identity,
                                    #$similarity,
                                    #$coverage,
                                    $percent_identity,
                                    $percent_similarity,
                                    $percent_coverage,
                                 )
                          )."\n";
            }
            
        }

    }    
}
