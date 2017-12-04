#!/usr/bin/perl

##
## Quick script to fix the deflines that estscan produces
##
## fixes the problem of reused accessions for multiple predicted ORFs per input sequence
##

use File::Temp qw(tempfile);
use Cwd;
use File::Copy;

my $idc = {};
my $cwd = getcwd();

my $infile = shift @ARGV || die "Must provide an input file";

$infile =~ /\.faa$|\.fna$/ || die "Unrecognized file extension";

my $type = ($infile =~ /\.faa$/) ? 'polypeptide' : 'CDS';

open(my $infh, $infile) || die "Can't open file '$infile': $!";

my ($tempfh, $tempfile) = tempfile(DIR => $cwd);

while (<$infh>) {
    if (/^>(\S+)/) {
        my $acc = $1;
        $acc .= '-'.$type.'-'.++$idc->{$acc};
        s/^>\S+/>$acc/;
        print $tempfh $_;
    } else {
        print $tempfh $_;
    }
}
close $tempfh;

move($tempfile, $infile) || die "Couldn't move file '$tempfile' to '$infile': $!";
