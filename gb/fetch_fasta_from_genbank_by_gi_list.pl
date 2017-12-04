#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Carp;

my $genbank_db = ''; # = 'nucest';
my $alphabet   = 'dna';
my $group_size = 500;
my $input_list = '';
my $output = '';

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
                        'input_list|i=s'    =>  \$input_list,
                        'db=s'              =>  \$genbank_db,
                        'alphabet=s'        =>  \$alphabet,
                        'group_size|g=i'    =>  \$group_size,
#                        'output|o=s'        =>  \$output,
#                        'log|l=s'           =>  \$log,
                       );


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
                       
my $prefix = basename($input_list);
my $output_file = $prefix.".fsa";
unlink($output_file);

my @output_segments = (); 

## read in the list of GIs
my @ids = read_gi_list($input_list);

my $count = scalar(@ids);
print $debug "# of Sequences in list: $count\n";

for (my $i = 0; $i < $count; $i += $group_size) {
    my $j = $i + $group_size - 1;

    if ($j > $count - 1) { $j = $count - 1;}
    
    my $segment_file_name = $output_file.".$i..$j";
    
    my @group = @ids[$i .. $j];
    
    if (-e $segment_file_name && valid_fasta($segment_file_name, scalar(@group))) {
        push(@output_segments, $segment_file_name);
        next;
    } else {
        unlink $segment_file_name;
    }
   
    my $factory = new Bio::DB::EUtilities(
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
                                -eutil      => 'efetch',
                                -rettype    => 'fasta',
                                -history    => $history,
                            );

    my $retry = 0;
    RETRIEVE_SEQS:
    $factory->set_parameters(
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
    if (valid_fasta($segment_file_name, scalar(@group))) {
        push(@output_segments, $segment_file_name);
    } else {
        confess "Server error: $@.  Try again later" if $retry == 5;
        print $debug "Server error, fasta file did not validate, redo #".$retry++."\n";
        redo RETRIEVE_SEQS;
    }
        
}

## merge the output files
cat_output($output_file, @output_segments);

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
    my ($file, $expected_sequence_count) = @_;
    
    print $debug "Validating '$file'...";

    my $fasta = new Bio::SeqIO(
                                -file       => $file, 
                                -format     => 'fasta',
                                -alphabet   => 'dna',
                              );
    for (my $i = 0; $i < $expected_sequence_count; $i++) {
        eval{
            my $seq = $fasta->next_seq();
        };
        if ($@) {
            carp "Encountered an error reading a sequence from the fasta file '$file': $@";
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
