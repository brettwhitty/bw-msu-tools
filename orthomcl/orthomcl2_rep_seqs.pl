#!/opt/rocks/bin/perl 

=head1  NAME 

orthomcl2_rep_seqs.pl - Short description

=head1 SYNOPSIS

USAGE: orthomcl2_rep_seqs.pl 
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
use Cwd qw{abs_path};
use File::Basename;
use Bio::DB::Fasta;

my ($input, $output, $fastadir, $help);

GetOptions(
    'input|i=s'     => \$input,
    'output|o=s'    => \$output,
    'fastadir|f=s'  => \$fastadir,
    'help|h!'       => \$help,
);

if ($help) {
    pod2usage(verbose => 2);
}

my @fastafiles;
if (-d $fastadir) {
    $fastadir = abs_path($fastadir);
    @fastafiles = <$fastadir/*.fasta>;
} else {
    croak "Must provide path to the orthomcl fasta database files with --fastadir flag";
}

my $db = {};
foreach my $file(@fastafiles) {
    my $base = basename($file, '.fasta');
    $db->{$base} = new Bio::DB::Fasta($file);
}

my $infh;

open $infh, '<', $input or confess "Failed to open input file '$input' for reading: $!";

while (<$infh>) {
    chomp;

    my $rep_len = 0;
    my $rep_id = '';

    my ($cluster_id, $member_string) = split(/: /, $_, 2);
    my @members = split(/ /, $member_string);

    foreach my $member(@members) {
        my ($member_db, $member_acc) = split(/\|/, $member, 2);

        my $member_len = $db->{$member_db}->length($member);
        if ($member_len > $rep_len) {
            $rep_len = $member_len;
            $rep_id = $member;
        }
    }
    my ($rep_db, undef) = split(/\|/, $rep_id, 2);
    my $rep_seq = $db->{$rep_db}->seq($rep_id);
    $rep_seq =~ s/(.{1,60})/$1\n/g;
    print '>'.$cluster_id.' '.$rep_id."\n";
    print $rep_seq;
}
