#!/usr/bin/perl


## This is a wrapper for bp_genbank2gff3.pl to help it support large multi-record GenBank flat files
##
## This is a workaround for the memory leak which I can't track down that causes bp_genbank2gff3.pl
## to consume GB++ of memory when converting a large multi-record flat file
##
##
## usage:
##
## bp_genbank2gff3_wrapper.pl "normal command line" input_file.gbk
##

$| = 1;

use FindBin qw{ $RealBin };

my $bp_genbank2gff3 = $RealBin.'/bp_genbank2gff3.pl';

my $cmd_line = shift @ARGV || die "Must provide command line";
my $infile = shift @ARGV || die "Must provide input file";

open (my $in_fh, $infile) || die "Failed to open '$infile' for reading: $!";

my $out_fh;
my $out_flag = 0;
while (<$in_fh>) {
    unless($out_flag) {
        open($out_fh, "| $bp_genbank2gff3 $cmd_line --input -") || die "Failed to open pipe to script";
        $out_flag = 1;
    }

    print $out_fh $_;
   
    if (/^\/\/\s*$/) {
        close $out_fh;
        $out_flag = 0;
    }
}

if ($out_flag) {
    close $out_fh;
}
