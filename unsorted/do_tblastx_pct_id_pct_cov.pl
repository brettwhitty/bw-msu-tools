#!/usr/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use Data::Dumper;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;
use File::Basename;

my $query_file = shift @ARGV;
my $db_file = shift @ARGV;
my $identity_cutoff = shift @ARGV;
my $coverage_cutoff = shift @ARGV;
my $evalue_cutoff = shift @ARGV;

unless ($identity_cutoff) { $identity_cutoff = 95;}
unless ($coverage_cutoff) { $coverage_cutoff = 95;}
unless ($evalue_cutoff) { $evalue_cutoff = '1e-25';}

my $hits_list = basename($query_file,".fsa")."--vs--".basename($db_file,".fasta").".hits.list";
my $nohits_list = basename($query_file,".fsa")."--vs--".basename($db_file,".fasta").".nohits.list";

my $factory= Bio::Tools::Run::StandAloneBlast->new(
                                                    'program'   =>  'tblastx',
                                                    'database'  =>  "$db_file",
                                                    'e'         =>  $evalue_cutoff,
                                                    'a'         =>  4,  # 4 processors
                                                    'b'         =>  1,
                                                    'v'         =>  1,
                                                  ); 

my $fasta_file = Bio::SeqIO->new(-file => $query_file, -format => 'Fasta', );

open (my $ofh_hits, ">".$hits_list) || die "Couldn't open '$hits_list' for writing: $!";
open (my $ofh_nohits, ">".$nohits_list) || die "Couldn't open '$nohits_list' for writing: $!";

while (my $seq = $fasta_file->next_seq()) {
    my $blast_report = $factory->blastall($seq);
    
    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
        my $hit_count = 0;

        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
            my $coverage = sprintf("%.2f", $hit->frac_aligned_query());
            my $percent_coverage = $coverage * 100;

            my $total_length = 0;
            my $sum = 0;
            while( my $hsp = $hit->next_hsp ) {
                ## calculate average % identity across all hsps
                $total_length += $hsp->length('query');
                $sum += $hsp->frac_identical * $hsp->length('query');
            }
            my $identity = sprintf("%.2f", $sum/$total_length);
            my $percent_identity = $identity * 100;
        
            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff) {
                
                $hit_count++;
                
                print join("\t", (
                                    $query_name,
                                    $hit_name,
                                    $rank,
                                    $score,
                                    $sig,
                                    $identity,
                                    $coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
                                 ))."\n";
            }
            
        }

        if ($hit_count) {
            print $ofh_hits $query_name."\n"; 
        } else {
            print $ofh_nohits $query_name."\n";
        }
        
    }    
}
