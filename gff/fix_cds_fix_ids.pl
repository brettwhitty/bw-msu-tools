#!/opt/rocks/bin/perl 

=head1  NAME 

fix_cds.pl - Fixes CDS ordering in GFF files

=head1 SYNOPSIS

USAGE: fix_cds.pl 
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

my ($input, $output, $help);

GetOptions(
    'input|i=s'     => \$input,
    'output|o=s'    => \$output,
    'help|h!'       => \$help,
);

if ($help) {
    pod2usage(verbose => 2);
}

#my $infh;
#if (defined($input)) {
#    open $infh, '<', $input
#        or die $!;
#} else {
#    $infh = \*STDIN;
#}


my $outfh;
if (defined($output)) {
    open $outfh, '>', $output
        or die $!;

} else {
    $outfh = \*STDOUT;
}

my $infh;

open $infh, '<', $input
    or die $!;

my $cdses = {};
while (<$infh>) {
    chomp;
    my @t = split("\t", $_);
    if (scalar(@t) != 9) {
        next;
    }

    if ($t[2] ne 'CDS') {
        next;
    }

    $t[8] =~ /ID=([^;]+)/;
    my $id = $1;

    push(@{$cdses->{$t[0].'.'.$id}}, [ @t ]);
}

foreach my $id(keys %{$cdses}) {
    @{$cdses->{$id}} = sort {$a->[3] <=> $b->[3]} @{$cdses->{$id}};
}
my $seen = {};

open $infh, '<', $input
    or die $!;

my $ids = {};
while (<$infh>) {
    chomp;
    my @t = split("\t", $_);
    if (scalar(@t) != 9) {
        next;
    }

    $t[8] =~ /ID=([^;]+)/;

    ## if not a CDS feature, then make sure we don't have duplicate IDs in the GFF file
    my $id = $1;
    if ($t[2] ne 'CDS') {
        $ids->{$id}++;
        if ($ids->{$id} > 1) {
            my $new_id = $id.'.'.$ids->{$id};
            $t[8] =~ s/(ID=)[^;]+/$1$new_id/;
        }
    }

    if (defined($cdses->{$t[0].'.'.$id})) {

        if (! defined($seen->{$id})) {
            $seen->{$id} = 1;
        
            foreach my $row(@{$cdses->{$t[0].'.'.$id}}) {
                print $outfh join("\t", @{$row})."\n";
            }
        }
    } else {
        print $outfh join("\t", @t)."\n";
    }
}
