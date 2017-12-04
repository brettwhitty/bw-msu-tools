#!/opt/rocks/bin/perl

use strict;
use warnings;

use Cwd qw{ abs_path };
use File::Basename qw{ dirname basename };

my $window = 12;

my $infile = shift @ARGV || die "Please provide SNP table file";
my $snp_support = shift @ARGV;

unless ($infile =~ /\.snp\.txt$/) {
    die "Unexpected input filename '$infile'";
}

unless ($snp_support) {
    $snp_support = 2;
}

$infile = abs_path($infile);
my $outdir = dirname($infile);
my $outbase = basename($infile, ".txt");
my $outfile = "$outdir/$outbase.hc.txt";

## for logging
my $logfh = *STDERR;

open (my $outfh, ">$outfile") || die "Failed to open '$outfile' for writing: $!";

my $pos = {};

print $logfh "Processing '$infile' to filter SNPs\n";
my $infh;
open ($infh, $infile) || die "Failed opening file '$infile': $!";
while (<$infh>) {
    chomp;

    if (/^#/) { next; };    

    my (
        $taxon_id,
        $scientific_name,
        $put_id,
        $member_count,
        $seq_len,
        $snp_loc,
        $snp_cov,
        $ref_base,
        $alt_base,
        $a_count,
        $c_count,
        $g_count,
        $t_count,
       )
        = split("\t");


    my $alt_base_depth = ${ 
                            {
                                split(/:/, $a_count),
                                split(/:/, $c_count),
                                split(/:/, $g_count),
                                split(/:/, $t_count),
                            }
                          }{$alt_base};

    ## keep the SNP only if support for the alternative base is equal to
    ## or greater than the cutoff that's been set
    if ($alt_base_depth >= $snp_support) {
        push(@{$pos->{$put_id}}, $snp_loc);
    } else {
        print $logfh "Rejected SNP on '$put_id' at '$snp_loc', has depth '$alt_base_depth' < '$snp_support'\n";
    }
}

## filter SNPs on window
my $filter = {};
foreach my $id(sort keys %{$pos}) {
    my @out = ();

    ## get and sort the positions
    my @in = @{$pos->{$id}};
    @in = sort {$a <=> $b} @in;

    ## left and right distance will be the calculated distance
    ## between SNP positions
    my $left_dist = $window + 1;
    my $right_dist = $window + 1;
    while (my $this = shift @in) {
        my $next = shift @in;

        ## set right distance for this
        if (defined $next) {
            $right_dist = $next - $this;
            unshift(@in, $next);
        } else {
            $right_dist = $window + 1;
        }

        ## store SNP that passed window criteria
        if ($left_dist > $window && $right_dist > $window) {
            push(@out, $this);
        } else {
            print $logfh "Rejected SNP on '$id' at '$this', has adjacent SNP within '$window' bases\n";
        }
        ## set left distance for next
        if (defined $next) {
            $left_dist = $next - $this;
        } else {
            $left_dist = $window + 1;
        }
    }
    ## set filter hash
    foreach my $snp_loc(@out) {
        $filter->{$id}->{$snp_loc} = 1;
    }
}
$pos = {};

## second pass through the file to write the output
open ($infh, $infile) || die "Failed opening file '$infile': $!";
while (<$infh>) {
    my (
        $taxon_id,
        $scientific_name,
        $put_id,
        $member_count,
        $seq_len,
        $snp_loc,
        $snp_cov,
        $ref_base,
        $alt_base,
        $a_count,
        $c_count,
        $g_count,
        $t_count,
       )
        = split("\t");

    ## keep the header lines
    if (/^#/) {
        print $outfh $_;
    ## keep SNPs that have passed the filters
    } elsif ($filter->{$put_id}->{$snp_loc}) {
        print $outfh $_;
    }
}
