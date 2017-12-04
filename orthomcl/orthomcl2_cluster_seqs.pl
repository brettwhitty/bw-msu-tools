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

$| = 1;

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

if (! -d $output) {
    croak "Provide a directory for output fasta files with --output flag";
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

    my ($cluster_id, $member_string) = split(/: /, $_, 2);
    my @members = split(/ /, $member_string);
    
    my $output_file = $output.'/'.$cluster_id.'.fasta';

    open my $outfh, '>', $output_file or confess "Failed to open file '$output_file' for writing: $!";

    foreach my $member(@members) {
        my ($member_db, $member_acc) = split(/\|/, $member, 2);
   

        my $member_seq = $db->{$member_db}->seq($member);
        if (! defined($member_seq)) {
            croak "Unable to retrieve sequence for '$member' from DB '$member_db'";
        }
        $member_seq =~ s/(.{1,60})/$1\n/g;
        print $outfh '>'.$member_acc."\n";
        print $outfh $member_seq;
    }
}
