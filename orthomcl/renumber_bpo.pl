#!/opt/rocks/bin/perl 

=head1  NAME 

renumber_bpo.pl - Renumbers the lines in a bpo file

=head1 SYNOPSIS

USAGE: renumber_bpo.pl 
        --input=/path/to/input
        --output=/path/to/output

=head1 OPTIONS

B<--input,-i> 
    Description of input.

B<--output,-o> 
    Description of output.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

For renumbering the lines in an OrthoMCL bpo file that has been prepared in chunks and merged.

=head1 INPUT

Merged bpo file.

=head1 OUTPUT

Renumbered bpo file.

=head1 CONTACT

Brett Whitty
whitty@msu.edu

=cut

use strict;
use warnings;

use Carp;
use Pod::Usage;
use Getopt::Long;

my ($input, $output, $help);

GetOptions(
    'input|i=s'     => \$input,
    'output|o=s'    => \$output,
    'help|h!'       => \$help,
);

if ($help || ! $input) {
    pod2usage(verbose => 2);
}

my $infh;
open $infh, '<', $input or die "$!";

my $outfh;
if ($output) {
    open $outfh, '>', $output or die "$!";
} else {
    $outfh = \*STDOUT;
}

my $count = 0;
while (<$infh>) {
    $count++;
    s/^\d+;/$count;/;
    print $_;
}
