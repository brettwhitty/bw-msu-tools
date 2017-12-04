#!/opt/rocks/bin/perl

=head1  NAME 

ngap.pl - find N-gap positions in nucleotide sequences

=head1 SYNOPSIS

USAGE: ngap.pl 
        --input=/path/to/fasta_file.fsa
        --output=/path/to/output.gff

=head1 OPTIONS

B<--input,-i> 
    Input FASTA or multi-FASTA format file.

B<--output,-o> 
    Output GFF file (will be created, must not exist)

B<--min_len,-l> 
    Minimum length of N-gap regions to output
    
B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to identify N-gap regions in FASTA format nucleotide sequences.

=head1 INPUT

The input should be an individual or multiple nucleotide sequences in one FASTA format file.
Sequence identifiers are recognized as the first string of non-whitespace characters occurring 
after the '>' character.

=head1 OUTPUT

The output of this script is a GFF3 format file.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'min_len|l=i',
              'help|h') || pod2usage();

if (scalar keys(%options) < 1) {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

my $next_id = 1;

## open the input file for parsing
open (IN, $options{'input'}) || die("can't open input file for reading");

my $seq = '';
my $seq_id = '';
my %ngaps = ();

## Read in FASTA sequences and find Ngaps
while (<IN>) {
    chomp;
    ## skip comments
    if (/^#/) { next; }
    ## skip blank lines
    if (/^\s+$/) { next; }
    
    if (/^>([^\s]+)/) {
        if ($seq ne '') {
            $ngaps{$seq_id} = &find_ngaps(\$seq);
        }
        $seq = '';
        $seq_id = $1;
    } else {
        $seq .= $_;
    }
}
if ($seq ne '') {
    $ngaps{$seq_id} = &find_ngaps(\$seq);
}
close IN;

## Write GFF output file

foreach $seq_id(keys(%ngaps)) {
   
    foreach my $ngap_ref(@{$ngaps{$seq_id}}) {

        if ($options{'min_len'} && ($ngap_ref->[1] - $ngap_ref->[0]) < $options{'min_len'}) {
            next;
        }
        
        print join("\t", (
                        $seq_id,
                        'ngap',
                        'gap', ## SO term
                        ($ngap_ref->[0] + 1), ## convert from interbase
                        $ngap_ref->[1],
                        '.',
                        '.',
                        '.',
                        'ID='.$seq_id.'-ngap'.$next_id++.';Name=N-gap;Length='.($ngap_ref->[1] - $ngap_ref->[0]),
                         ))."\n";
    
    }
}

exit(0);

sub find_ngaps {
    my ($seq_ref) = shift @_;
    my @ngaps = ();
    
    while (${$seq_ref} =~ /[N]{1,}/ig) {
        push (@ngaps, [$-[0],$+[0]]);
    }
    return \@ngaps;
}

sub check_parameters {
    ## check required options
    my @required = qw( input );
    for ( @required ) {
        unless ( defined $options{$_} ) {
                die("--$_ is a required parameter");
        }
    }

    ## make sure input file exists
    if (! -e $options{'input'}) { die("input file $options{'input'} does not exist") }
    
    return 1;
}
