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


my ($query_file, $db_file, $identity_cutoff, $similarity_cutoff, $coverage_cutoff, $evalue_cutoff, $b, $v, $cpu);

$b = 1;
$v = 1;
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
                       );

unless ($identity_cutoff) { $identity_cutoff = 30;}
unless ($similarity_cutoff) { $similarity_cutoff = 30;}
unless ($coverage_cutoff) { $coverage_cutoff = 20;}
unless ($evalue_cutoff) { $evalue_cutoff = '1e-5';}

## will be used to store the RBHs
my $rbh = {};

my $factory= Bio::Tools::Run::StandAloneBlast->new(
                                                    'program'   =>  'blastp',
                                                    'database'  =>  "$db_file",
                                                    'e'         =>  $evalue_cutoff,
                                                    'a'         =>  $cpu,  # 4 processors
                                                    'b'         =>  $b,
                                                    'v'         =>  $v,
                                                    'F'         =>  '"m S"',
                                                    's'         =>  'T',
                                                  ); 

my $s_factory = Bio::Tools::Run::StandAloneBlast->new(
                                                    'program'   =>  'blastp',
                                                    'database'  =>  "$query_file",
                                                    'e'         =>  $evalue_cutoff,
                                                    'a'         =>  $cpu,  # 4 processors
                                                    'b'         =>  $b,
                                                    'v'         =>  $v,
                                                    'F'         =>  '"m S"',
                                                    's'         =>  'T',
                                                  ); 
                                                  
my $fasta_file = Bio::SeqIO->new(-file => $query_file, -format => 'Fasta', );

my $search_db = new Bio::DB::Fasta($db_file);

#tie my %search_db, 'Bio::DB::Fasta', $db_file;

while (my $seq = $fasta_file->next_seq()) {
    if ($seq->length <= 40) { next; }

    my $blast_report = $factory->blastall($seq);
 
    $rbh = {};
    
    process_blast($blast_report);
    
    print_rbh();
}

sub process_blast {
    my ($blast_report, $forward_query_name) = @_;

    while( my $result = $blast_report->next_result ) {
        
        my $query_name = $result->query_name;
        my $hit_count = 0;
        my $prev_score = -1;
        my $prev_sig = -1;
        
        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;
            my $coverage = sprintf("%.2f", $hit->frac_aligned_query());
            my $percent_coverage = $coverage * 100;

            if ($prev_sig != -1) {
                unless ($score eq $prev_score && $sig eq $prev_sig) {
                    last;
                }
            }
            $prev_score = $score;
            $prev_sig = $sig;
       
            my $total_length = 0;
            my $id_sum = 0;
            my $sim_sum = 0;
            while( my $hsp = $hit->next_hsp ) {
                ## calculate average % identity across all hsps
                $total_length += $hsp->length('query');
                $id_sum += $hsp->frac_identical * $hsp->length('query');
                $sim_sum += $hsp->frac_conserved * $hsp->length('query');
            }
            my $identity = sprintf("%.2f", $id_sum/$total_length);
            my $similarity = sprintf("%.2f", $sim_sum/$total_length);
            my $percent_identity    = $identity * 100;
            my $percent_similarity  = $similarity * 100;
        
            if ($percent_coverage >= $coverage_cutoff && $percent_identity >= $identity_cutoff && $percent_similarity >= $similarity_cutoff) {
                
                $hit_count++;
              
                unless ($forward_query_name) {
                    ## fetch hit seq
                    my $hit_seq = $search_db->get_Seq_by_id($hit_name);
        
                    my $s_blast_report = $s_factory->blastall($hit_seq);
                    
                    process_blast($s_blast_report, $query_name);
                }
    
                if ($forward_query_name) {
                
                    if ($hit_name eq $forward_query_name) {
                    
                        $rbh->{$forward_query_name}->{'reverse'}->{$query_name} =  
                                 [
                                 #   $hit_name,
                                    #    $rank,
                                    $score,
                                    $sig,
                                    $identity,
                                    $similarity,
                                    $coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
                                 ];
                    }
                } else {
                    $rbh->{$query_name}->{'forward'}->{$hit_name} =  
                                 [
                                    $hit_name,
                                    #   $rank,
                                    $score,
                                    $sig,
                                    $identity,
                                    $similarity,
                                    $coverage,
                                    #$percent_identity.'%',
                                    #$percent_coverage.'%',
                                 ];
         
                }
            }
            
        }

    }    
}

sub print_rbh {
#print Dumper $rbh;
#return;
    foreach my $q(keys(%{$rbh})) {
        foreach my $hit(keys(%{$rbh->{$q}->{'forward'}})) {
            print $q."\t";
            print join("\t", @{$rbh->{$q}->{'forward'}->{$hit}})."\t";
            if (defined($rbh->{$q}->{'reverse'}->{$hit})) {
                print join("\t", @{$rbh->{$q}->{'reverse'}->{$hit}})."\n";
            } else {
                print "\t" x 4 . "\n";
            }
        } 
    }
}
