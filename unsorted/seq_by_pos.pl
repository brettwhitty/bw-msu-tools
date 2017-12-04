#!/opt/rocks/bin/perl 

=head1  NAME 

seq_pos.pl - Short description

=head1 SYNOPSIS

USAGE: seq_pos.pl 
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
use Bio::DB::Fasta;

my ($id, $output, $help, $fasta, $start, $end, $revcomp);

GetOptions(
#    'input|i=s'     => \$input,
#    'output|o=s'    => \$output,
    'id|i=s'        =>  \$id,
    'fasta|f=s'        =>  \$fasta,
    'help|h!'       => \$help,
    'start|s=i'     =>  \$start,
    'end|e=i'       =>  \$end,
    'revcomp|r!'    =>  \$revcomp,
);

if ($help || ! defined($fasta)) {
    pod2usage(verbose => 2);
}

my $db = new Bio::DB::Fasta($fasta);

my $seq = $db->seq($id, $start, $end);
if ($revcomp) {
    $seq = reverse_complement_dna($seq);
}
print $seq;

sub reverse_complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
    $r_seq = reverse($r_seq);

    return $r_seq;
}
