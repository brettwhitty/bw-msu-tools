#!/usr/bin/perl

## creates a DB_File file of a tab-delimited table of key/value pairs

use strict;
use warnings;

use DB_File;
use File::Basename;
use Cwd qw( abs_path );

my $src_file = shift @ARGV || die "Expects key/value table source file as argument";

unless (-e $src_file) {
    die "Must provide a source file that exists";
}

unless ($src_file =~ /\.txt$/) {
    die "Expect that source file is a tab delimited key/value table ending with the extension '.txt'";
}

$src_file = abs_path($src_file);

my $base = basename($src_file, '.txt');
my $dir = dirname($src_file);

my $db_filename = "$dir/$base.txt.db";

## remove the database if it exists
unlink($db_filename);

## tie the db file
my %db;
tie %db, "DB_File", $db_filename
    or die "Can't tie the db file '$db_filename': $!\n";

open (IN, $src_file) || die "Failed to open file '$src_file' for reading: $!";

my $counter = 0;
while (<IN>) {
    chomp;
    
    $counter++;
    
    my @t = split(/\t/, $_);

    ## don't expect this to happen, but we'll ignore it if it does
    if (scalar(@t) != 2) {
        print STDERR "WARNING: skipping line $counter\n";
        next;
    }

    $db{$t[0]} = $t[1];
}
untie %db;
