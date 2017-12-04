#!/usr/bin/perl

## quick script to do tblastx with cutoffs for percent identity and percent coverage
## of the query sequence for all HSPs against a given subject sequence

use warnings;
use strict;

use lib "/opt/rocks/lib/perl5/site_perl";

use Data::Dumper;
use Bio::SearchIO;
use File::Basename;
use Getopt::Long;

my ($result_file, $identity_cutoff, $coverage_cutoff, $evalue_cutoff);

GetOptions(
    'input|i=s' =>  \$result_file,
    'id=i'      =>  \$identity_cutoff,
    'cov=i'     =>  \$coverage_cutoff,
    'e=s'       =>  \$evalue_cutoff,
);

$identity_cutoff ||= 1;
$coverage_cutoff ||= 1;
$evalue_cutoff   ||= '10';

my $infh;

if (! defined($result_file)) {
    $infh = \*STDIN;
} elsif ($result_file =~ /\.gz/) {
    open $infh, '<:gzip', $result_file || die "Failed to open '$result_file': $!";
} else {
    open $infh, '<', $result_file || die "Failed to open '$result_file': $!";
}

    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -fh   =>  $infh,
                                        );
    

    my $query_len = 0;
    my $aln_index = 0;                                        
    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
        my $query_len = $result->query_length;
        my $hit_count = 0;

        while( my $hit = $result->next_hit ) {
        
            $aln_index++;

            my $hit_name = $hit->name;
            my $hit_len = $hit->length;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
            my $coverage = sprintf("%.2f", $hit->frac_aligned_query());
            my $percent_coverage = $coverage * 100;

            my $total_length = 0;
            my $sum = 0;
            my $simsum = 0;
            my $hsp_index = 0;
            my $hsp_string = '';
            while( my $hsp = $hit->next_hsp ) {

                $hsp_index++;

                my $q_start = $hsp->start('query');
                my $q_end   = $hsp->end('query');
#                my $q_len   = $hsp->length('query');
                my $s_start = $hsp->start('hit');
                my $s_end   = $hsp->end('hit');
#                my $s_len   = $hsp->length('hit');

                $hsp_string .= "$hsp_index:$q_start-$q_end:$s_start-$s_end.";

                ## calculate average % identity across all hsps
                $total_length += $hsp->length('query');
                $sum += $hsp->frac_identical * $hsp->length('query');
                $simsum += $hsp->frac_conserved * $hsp->length('query');
            }
            my $identity = sprintf("%.2f", $sum/$total_length);
            my $similarity = sprintf("%.2f", $simsum/$total_length);
            my $percent_identity = $identity * 100;
            my $percent_similarity = $similarity * 100;
        
#            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff) {
                
                $hit_count++;
                
                print join(";", (
                                    $aln_index,
                                    $query_name,
                                    $query_len,
                                    $hit_name,
                                    $hit_len,
                                    #$score,
                                    $sig,
                                    $identity,
                                    #$similarity,
                                    #$coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
                                    $hsp_string,
                                 ))."\n";
#            }
            
        }

#        if ($hit_count) {
            #print STDOUT $query_name."\n"; 
#        } else {
            #print STDERR $query_name."\n";
#            print $query_name."\t-\n";
#        }
        
    }    
