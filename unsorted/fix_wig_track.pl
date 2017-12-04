#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

use lib '/home/whitty/SVN/lib';
use MyIO;

my ($name, $desc, $input, $output);

GetOptions(
	'input|i=s'	=>	\$input,
	'output|o=s'	=>	\$output,
	'name|n=s'	=>	\$name,
	'desc|d=s'	=>	\$desc,
);

my $infh = get_infh($input);
my $outfh = get_outfh($output);

while (<$infh>) {
	if (/^track type=bedGraph/) {
		print $outfh qq{track type=wiggle_0 name="$name" description="$desc" trim=stdev1 transform=logsquared\n};
	} else {
		print $outfh $_;
	}
}
