#!/opt/rocks/bin/perl 

=head1  NAME 

to_snp_table.pl - Short description

=head1 SYNOPSIS

USAGE: to_snp_table.pl 
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

Description of script.

=head1 INPUT

Description of input.

=head1 OUTPUT

Description of output.

=head1 CONTACT

Brett Whitty
whitty@msu.edu

=cut

use strict;
use warnings;

use Carp;
use Pod::Usage;
use Getopt::Long;

my ($help,$input, $org, $ref_genome, $snp_genome);

GetOptions(
    'input|i=s'         => \$input,
    'org|o=s'           => \$org,
    'ref_genome|r=s'    => \$ref_genome,
    'snp_genome|s=s'    => \$snp_genome,
    'help|h!'           => \$help,
);

if ($help) {
    pod2usage(verbose => 2);
}

open my $infh, '<', $input or croak "Failed to open '$input' for reading: $!";
while (<$infh>) { 
    chomp;

    my @t = split("\t", $_);

    print join("\t", (
            $org,
            $ref_genome,
            $snp_genome,
            $t[0],
            $t[1],
            $t[2],
            $t[3],
        ))."\n";

}
