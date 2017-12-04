#!/opt/rocks/bin/perl 

=head1  NAME 

make_chrUn.pl - Short description

=head1 SYNOPSIS

USAGE: make_chrUn.pl 
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

my ($fasta, $spacer, $list, $input, $output, $help);

GetOptions(
    'fasta|f=s'     => \$fasta,
    'list|l=s'      => \$list,
#    'agp|s=s'       => \$agp,
    'spacer|s=i'    =>  \$spacer,
#    'output|o=s'    => \$output,
    'help|h!'       => \$help,
);

$spacer ||= 100;

if ($help) {
    pod2usage(verbose => 2);
}

open my $infh, '<', $list or confess $!;

my $db = new Bio::DB::Fasta($fasta);

my $lens = ();
my @ids = ();
while (<$infh>) {
    chomp;

    push(@ids, $_);
    $lens->{$_} = $db->length($_) or confess "'$_' not in sequence database";
}

my $outseq;
my $seq_count = scalar(@ids);

my $offset = 0;
my $agp_counter = 0;
for (my $i = 0; $i < $seq_count; $i++) {

    print STDERR join("\t", (
            'chrUn',
            $offset + 1,
            $offset + $lens->{$ids[$i]},
            $agp_counter++,
            'W',
            $ids[$i],
            1,
            $lens->{$ids[$i]},
            '+',
        ))."\n";
    $offset += $lens->{$ids[$i]};

    $outseq .= $db->seq($ids[$i]);

    if ($i != $seq_count - 1) {
        print STDERR join("\t", (
            'chrUn',
            $offset + 1,
            $offset + $spacer,
            $agp_counter++,
            'N',
            $spacer,
            'contig',
            'no',
            '',
        ))."\n";
        $offset += $spacer;
        $outseq .= 'N' x $spacer;
    }
}

    $outseq =~ s/(.{1,60})/$1\n/g;
    print STDOUT ">chrUn\n".$outseq;

