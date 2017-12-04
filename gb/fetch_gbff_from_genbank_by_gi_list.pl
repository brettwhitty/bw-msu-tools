#!/usr/bin/env perl

$| = 1;

use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use Bio::DB::EUtilities;
use Bio::SeqIO;
#use DateTime;
#use Date::Format;
use Time::Piece;
use Carp;
use Set::Scalar;

my $genbank_db = ''; # = 'nucest';
my $alphabet   = 'dna';
my $group_size = 100;
my $input_list = '';
my $output = '';
my $log = '';
my $email = 'nobody@nosuch.com';

my $debug = \*STDERR;

my %valid_db = ( 
                    gene => 1,
                    genome => 1,
                    nucleotide => 1,
                    nuccore => 1,
                    nucest => 1,
                    nucgss => 1,
                    protein => 1,
                    popset => 1,
                    snp => 1,
                    sequences => 1,
                );
                
my $result = GetOptions(
                        'email=s'           =>  \$email,
                        'input_list|i=s'    =>  \$input_list,
                        'db=s'              =>  \$genbank_db,
                        'alphabet=s'        =>  \$alphabet,
                        'group_size|g=i'    =>  \$group_size,
                        'output|o=s'        =>  \$output,
                        'log|l=s'           =>  \$log,
                       );

print $debug "NOTICE: This script has been modified to use GB accession numbers instead of now obsoleted GIs.\n";                       

unless ($input_list) {
    confess "Must provide an input list of GI numbers with --input_list flag";
}
unless (-f $input_list) {
    confess "Specified input list '$input_list' doesn't exist!";
}                      
unless ($genbank_db) {
    confess "Must provide a valid genbank db to search using --db flag ('".join("', '", keys(%valid_db))."')";
}
unless ($valid_db{$genbank_db}) {
    confess "Must provide a valid genbank db to search using --db flag ('".join("', '", keys(%valid_db))."')";
}
                    
if ($log) {
    open ($debug, ">>$log") || confess "Couldn't open log file '$log' for writing: $!";
}

my $prefix = basename($input_list);
my $output_file;
if ($output) {
    $output_file = $output;
} else {
    $output_file = $prefix.".gbk";
}
unlink($output_file);

my @output_segments = (); 

my $begin_time = localtime();
#print $debug time2str("BEGIN: %a %b %e %T %Y\n", localtime(time()));;
print $debug "BEGIN: ".$begin_time->datetime."\n";
print $debug "Fetch by list file: $input_list\n";
## read in the list of GIs
my @ids = read_gi_list($input_list);

my $count = scalar(@ids);
print $debug "# of Sequences: $count\n";

for (my $i = 0; $i < $count; $i += $group_size) {
    my $j = $i + $group_size - 1;

    if ($j > $count - 1) { $j = $count - 1;}
    
    my $segment_file_name = $output_file.".$i..$j";
    
    my @group = @ids[$i .. $j];

## DEBUG   
#    open(OUT, ">$segment_file_name.list");
#    print OUT join("\n", @group)."\n";
#    close OUT;
## DEBUG

    if (-e $segment_file_name && valid_fasta($segment_file_name, \@group)) {
        push(@output_segments, $segment_file_name);
        next;
    } else {
        unlink $segment_file_name;
    }
   
    my $factory = new Bio::DB::EUtilities(
                                            -email  => $email,
                                            -eutil  => 'epost',
                                            -db     => $genbank_db,
                                            -id     => \@group,
                                            -keep_histories => 1,
                                         );

                                         
    my $history = $factory->next_History();
    
    print "Posted successfully\n";
    print "WebEnv    : ", $history->get_webenv, "\n";
    print "Query_key : ", $history->get_query_key, "\n";

    $factory->set_parameters(
                                -email      => $email,
                                -eutil      => 'efetch',
                                -rettype    => 'gb',
                                -retmode    => 'text',
                                -history    => $history,
                            );

    my $retry = 0;
    my $flag = 0;
    RETRIEVE_SEQS: while ($flag == 0) {
        $factory->set_parameters(
                                -email      => $email,
                                -retmax     => $group_size,
                                -retstart   => 0,
                            );
        eval{
            $factory->get_Response(
                                -file   =>  ">$segment_file_name",
                              );
        };
        if ($@) {
            confess "Server error: $@.  Try again later" if $retry == 5;
            print $debug "Server error, redo #".$retry++."\n";
            redo RETRIEVE_SEQS;
        }
        if (valid_fasta($segment_file_name, \@group)) {
            push(@output_segments, $segment_file_name);
        } else {
            confess "Server error: $@.  Try again later" if $retry == 5;
            print $debug "Server error, fasta file did not validate, redo #".$retry++."\n";
            redo RETRIEVE_SEQS;
        }
        $flag = 1;
    };
}

## merge the output files
cat_output($output_file, @output_segments);
my $end_time = localtime();
print $debug "END: ".$end_time->datetime."\n";

## read the list of GIs from the list file
sub read_gi_list {
    my ($input_list) = @_;
    
    open (IN, $input_list) || confess "Failed opening input list '$input_list' for reading: $!";

    my @ids = ();
    while (<IN>) {
        chomp;
 
        push (@ids, $_);   
    }

    return @ids;

}


## uses Bio::SeqIO to do a simple validation test on the fasta file
sub valid_fasta {
    my ($file, $group_ref) = @_;

    my $expected_sequence_count = scalar(@{$group_ref});

    print $debug "Validating '$file'...";

    open (my $infh, $file) || die "Failed to open file '$file' for reading: $!";

    my $lens = {};
    my $record_count = 0;
    my @found_gi = ();
    while (<$infh>) {
        if (/^LOCUS\s+(\S+)\s+(\d+) (bp|aa|rc)/) {
            my $acc = $1;
            my $len = $2;
            my $unit = $3;
            if ($unit eq 'bp') {
                $lens->{$acc} = $len;
            } else {
                $lens->{$acc} = 0;
            }
            $record_count++;
        } elsif (/^VERSION\s+([^.]+)\S+\s+GI:(\d+)/) {
            push(@found_gi, $2);
        } elsif (/^VERSION\s+(\S+)/) { ## hack fix for GB not having GI #s anymore
            push(@found_gi, $1);
        }
    }

    my $query_set = new Set::Scalar(@{$group_ref});
    my $fetch_set = new Set::Scalar(@found_gi);

    my $missing_set = $query_set - $fetch_set;
    my $bonus_set = $fetch_set - $query_set;

    my $missing_count = $missing_set->size();
    my $bonus_count = $bonus_set->size();

    if ($missing_count > 0) {
        print $debug "FAIL!\n";
        carp "Found '$missing_count' missing GI #s in retrieved set: "
             .join(",", $missing_set->members);
        if ($bonus_set > 0) {
            carp "Found '$bonus_count' unexpected GI #s in retrieved set: "
             .join(",", $bonus_set->members);
        } else {
            return 0;
        }
    }

#    ## check that record count matches what we expect
#    if ($record_count != $expected_sequence_count) {
#        carp "Found '$record_count' records when '$expected_sequence_count' in '$file': $@";
#        return 0;
#    }

    my $fasta = new Bio::SeqIO(
                                -file       => $file, 
                                -format     => 'genbank',
                                -alphabet   => 'dna',
                            );

#    my $record_count = scalar(keys(%{$lens}));
    for (my $i = 0; $i < $record_count; $i++) {
        my $seq;
        eval{
            $seq = $fasta->next_seq();
        };
        if ($@) {
            carp "Encountered an error reading a sequence from the fasta file '$file': $@";
            return 0;
        }

        my $seq_len = (defined($seq->seq())) ? length($seq->seq()) : 0;
        
        ## check that sequence length matches what we've read from the header
        if ($lens->{$seq->id} != $seq_len) {
            carp "Sequence '".$seq->id."' length different from expected in '$file': $@";
            return 0;
        }
    }
    print $debug "OK\n";
    return 1;    
}

sub cat_output {
    my ($outfile, @segments) = @_;

    print $debug "Merging output files...";
    
    open (OUT, ">$outfile") || confess "Failed opening '$outfile' for write: $!";

    foreach my $infile(@segments) {
        open (IN, "<$infile") || confess "Failed opening '$infile' for read: $!";

        while (<IN>) {
            print OUT $_;
        }
    }
    foreach my $infile(@segments) {
        unlink($infile);
    }
    print $debug "done.\n";
}
