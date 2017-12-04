#!/opt/rocks/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use Data::Dumper;
use Bio::SearchIO;
use File::Basename;

my $result_file = shift @ARGV;
my $identity_cutoff = shift @ARGV;
my $coverage_cutoff = shift @ARGV;
my $evalue_cutoff = shift @ARGV;

unless ($identity_cutoff) { $identity_cutoff = 1;}
unless ($coverage_cutoff) { $coverage_cutoff = 1;}
unless ($evalue_cutoff) { $evalue_cutoff = '10';}


    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -file   =>  $result_file,
                                        );
    
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
            my $simsum = 0;
            while( my $hsp = $hit->next_hsp ) {
                ## calculate average % identity across all hsps
                $total_length += $hsp->length('query');
                $sum += $hsp->frac_identical * $hsp->length('query');
                $simsum += $hsp->frac_conserved * $hsp->length('query');
            }
            my $identity = sprintf("%.2f", $sum/$total_length);
            my $similarity = sprintf("%.2f", $simsum/$total_length);
            my $percent_identity = $identity * 100;
            my $percent_similarity = $similarity * 100;
        
            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff) {
                
                $hit_count++;
                
                print join("\t", (
                                    $query_name,
                                    $hit_name,
                                    $rank,
                                    $score,
                                    $sig,
                                    $identity,
                                    $similarity,
                                    $coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
                                 ))."\n";
            }
            
        }

        if ($hit_count) {
            #print STDOUT $query_name."\n"; 
        } else {
            #print STDERR $query_name."\n";
            print $query_name."\t-\n";
        }
        
    }    
