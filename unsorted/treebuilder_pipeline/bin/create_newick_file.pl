#! /usr/bin/perl

# create_newick_file.pl

# 10 August 2009
# Kevin Childs

# This script will parse a protdist distance file and will create a newick file.
# This script is only designed to work with two node trees that cannot be created
# by proml.

use Getopt::Std;

use strict;

my $usage = "\n$0 -i input_file  -o outtree\n\n";

our ( $opt_i, $opt_o, $opt_h );
getopts("i:o:h") or die usage();

if ($opt_h) {
    print $usage;
    exit;
}

my $input_file = $opt_i;
my $output_file = $opt_o;

if (   !defined($input_file)
       || !( -e $input_file )
       || !defined($output_file)
       || (-e $output_file)) {
    die "\nMissing or invalid input values.\n" . $usage;
}

open IN, $input_file || die "\nUnable to open $input_file for reading.\n\n";
open OUT, ">$output_file" || die "\nUnable to open $output_file for writing.\n\n";
my $line = <IN>;
$line = <IN>;
chomp $line;
$line =~ s/\s+/ /g;
my @elems = split " ", $line;
my $dist = $elems[2];
my $half_dist = sprintf("%.5f", ($dist / 2));
print OUT "(1:$half_dist,2:$half_dist);";
close IN;
close OUT;

exit;

