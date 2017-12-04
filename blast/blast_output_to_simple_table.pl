#!/usr/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use lib "/home/whitty/SVN/lib";

use MyIO;

use Data::Dumper;
use Bio::SearchIO;
use File::Basename;

my $result_file = shift @ARGV;
my $evalue_cutoff = shift @ARGV;

unless ($evalue_cutoff) { $evalue_cutoff = '1e-5';}

    my $infh = get_infh($result_file);

    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -fh   =>    $infh,
                                        );
    
    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
        my $hit_count = 0;

        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
                
                print join("\t", (
                                    $query_name,
                                    $hit_name,
                                    $rank,
                                    $score,
                                    $sig,
                                 ))."\n";
            }
            
    }

