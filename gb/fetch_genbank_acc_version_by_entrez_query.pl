#!/usr/bin/env perl

## uses Bio::DB::EUtilities to fetch a list of genbank GI numbers by Entrez query
##
## Brett Whitty 2008
## whitty@msu.edu

use lib "/home/whitty/lib/perl5/site_perl/5.8.8";

use warnings;
use strict;
use Carp;

use Getopt::Long;
use Bio::DB::EUtilities;
use Time::Piece;

my $log = '';
my $output = '';
my $genbank_db = ''; # = 'nucest';
my $search_terms = ''; # = 'txid4070[Organism:exp]';
my $email = 'nobody@nosuch.com';

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
                        'email=s'       =>  \$email,
                        'db=s'          =>  \$genbank_db,
                        'query=s'       =>  \$search_terms,
                        'output|o=s'    =>  \$output,
                        'log|l=s'       =>  \$log,
                       );

unless ($genbank_db) {
    confess "Must provide a valid genbank db to search using --db flag ('".join("', '", keys(%valid_db))."')";
}
unless ($valid_db{$genbank_db}) {
    confess "Must provide a valid genbank db to search using --db flag ('".join("', '", keys(%valid_db))."')";
}
unless ($search_terms) {
    confess "Must provide search string with --query flag";
}

my $debug = \*STDERR;
if ($log) {
    open ($debug, ">".$log) || confess "Failed opening log file '$log' for writing: $!";
}

my $out = \*STDOUT;
if ($output) {
    open ($out, ">".$output) || confess "Failed opening output file '$output' for writing: $!";
    print $debug "Will write output to file '$output'.\n";
}
    
my $begin_time = localtime();
print $debug "BEGIN: ".$begin_time->datetime()."\n";
print $debug "Query string is: \"$search_terms\"\n";
print $debug "Database to search is: $genbank_db\n";

print $debug 'Querying Entrez...';
my $factory = Bio::DB::EUtilities->new(
                                        -email      => $email,
                                        -eutil      => 'esearch',
                                        -db         => $genbank_db,
                                        -term       => $search_terms,
                                        -usehistory => 'yes',
                                      );
print $debug "done.\n";
                                      
my $count = $factory->get_count();
print $debug "$count records found.\n";

my @accessions = ();
my $retstart = 0;
my $remaining = $count;
my $fetch_count = 0;
if ($count > 0) {
    my $history = $factory->next_History() or confess 'No history data returned';

    while ($remaining > 0 ) {

        print $debug "Doing fetch for $count records [retstart=$retstart, retmax=10000] ... ";
        $factory->set_parameters(
                                        -email      => $email,
                                        -eutil      => 'efetch',
                                        -db         => $genbank_db,
                                        -rettype    => 'acc',
                                        -retmode    => 'text',
                                        -history    => $history,
                                        -retstart   => $retstart,
                                        -retmax     => 10000,
                                      );

        my $result = $factory->get_Response->content;
        chomp $result;
        
        my @chunk_accessions = split(/\n/, $factory->get_Response->content);

        $fetch_count = scalar(@chunk_accessions);

        print $debug "results=$fetch_count ... ";

        if ($fetch_count != $remaining && $fetch_count != 10000) {
            confess "Fetch [retstart=$retstart] returned an unexpected number of records!\n";
        }
        
        push(@accessions, @chunk_accessions);
        
        $remaining -= $fetch_count;

        if ($remaining == 0) {
            print $debug "success. Fetch complete.\n";
        } elsif ( $remaining > 0 && $fetch_count == 10000) {
            print $debug "success. Remaining=$remaining.\n";
        } else {
            confess "This should never happen! [retstart=$retstart]\n";
        }

        ## next chunk
        $retstart += 10000;
    }
    if (scalar(@accessions) != $count) {
        confess "Search returned $count GIs, but only ".scalar(@accessions)." were retrieved\n";
    }
} else {
    print $debug "Nothing to fetch...done.\n";
}
my $end_time = localtime();
print $debug "END: ".$end_time->datetime()."\n";

foreach my $accession(@accessions) {
    print $out $accession."\n";
}
