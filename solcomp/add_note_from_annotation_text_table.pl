#!/usr/bin/perl

use strict;
use warnings;

## this script will add a Note to a GFF3 file based on the value of the Name attribute
## and given a table of functional text that has been added to a DB_File database using:
## ~/SVN/init_uniref_db_file_db.pl tab_delimited_key_value_table.txt
##
## it's a temporary script that needs to be rewritten
##
## usage:
##
## add_note_from_annotation_text_table.pl DB_File-filename.db GFF3_file
##
## Brett Whitty, whitty@msu.edu

use lib '/home/whitty/SVN/lib';
use GFFTextEncoder;

use DB_File;
use Carp;

my $db_file = shift @ARGV;
my $gff_file = shift @ARGV;

tie(my %db, 'DB_File', $db_file, O_RDONLY, 0600) || croak "Failed to tie db file '$db_file'";
open(my $infh, '<', $gff_file) || croak "Failed to open file '$gff_file': $!";

while (<$infh>) {
    chomp;

    my $name;
    if (/Name=([^;]+)/) {
        $name = $1;
        if (defined($db{$name})) {
            if (/Note=/) {
                croak "Note attribute already defined in GFF";
            }
            my $note = encode($db{$name});
            $_ =~ s/(Name=[^;]+)/$1;Note=$note/;
        }
        print $_."\n";
    } else {
        print $_."\n";
    }

}
