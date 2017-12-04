#!/opt/rocks/bin/perl 

=head1  NAME 

recall_samtools_pileup_snps.pl - Replaces ambiguous SNP base calls with true SNP base where possible

=head1 SYNOPSIS

USAGE: recall_samtools_pileup_snps.pl 
        --input=/path/to/input
        --output=/path/to/output

=head1 OPTIONS

B<--input,-i> 
    SNP pileup output from samtools.

B<--output,-o> 
    SNP pileup output from samtools with modified base call.

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

my ($input, $output, $help, $min_depth, $consensus_pct);

GetOptions(
    'input|i=s'     => \$input,
    'output|o=s'    => \$output,
    'min_depth|d=i' => \$min_depth,
    'consensus_percent|p=i' => \$consensus_pct,
    'help|h!'       => \$help,
);

if (! defined($min_depth)) {
    $min_depth = 0;
}
if (! defined($consensus_pct)) {
    $consensus_pct = 60;
}

if ($help) {
    pod2usage(verbose => 2);
}

my $infh;
if (defined($input)) {
    open $infh, '<', $input or confess "Failed to open '$input' for reading: $!";
} else {
    $infh = \*STDIN;    
}

my $outfh;

if (defined($output)) {
    open $outfh, '>', $output or confess "Failed to open '$output' for writing: $!";
} else {
    $outfh = \*STDOUT;
}

while (<$infh>) {
    chomp;
    
    my $skip = 0; 

    my @t = split("\t", $_);

    my $base_pileup = $t[8];
    $base_pileup =~ s/[\$.,n]//ig;
    $base_pileup =~ s/\^.//g;
    $base_pileup = uc($base_pileup);

    my $snp_base_total = length($base_pileup);

    my $base_count = {};
    foreach my $base(split("", $base_pileup)) {
        $base_count->{$base}++;
    }

    my @base_freqs = ();
    foreach my $base(keys(%{$base_count})) {
        push(@base_freqs, { 'base' => $base, 'count' => $base_count->{$base} });
    }

    @base_freqs = sort{ $b->{'count'} <=> $a->{'count'} } @base_freqs;

    my $most_frequent_base = $base_freqs[0];

    ## skip this SNP if the most frequent base occurs < our min_depth cutoff    
    if ($min_depth && $most_frequent_base->{'count'} < $min_depth) {
       $skip = 1; 
    }

    my $pct_frequent_base = int($most_frequent_base->{'count'} / $snp_base_total * 100);

    ## if the most frequent base represents > consensus_pct % of the SNP bases in the pileup,
    ## we'll ignore the other bases
    if ($consensus_pct && $pct_frequent_base >= $consensus_pct) {
        ## eliminate other bases from consideration
        $base_count = { $most_frequent_base->{'base'} => $base_count->{$most_frequent_base->{'base'}} };
    }

    my @bases = keys(%{$base_count});   
    ## only one SNP base, so replace ambiguous SNP base
    if (scalar(@bases) == 1) {
        $t[3] = $bases[0];
    } 
    
    unless ($skip) {
        print $outfh join("\t", @t)."\n";
    }
}


