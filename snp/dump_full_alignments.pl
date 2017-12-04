#!/usr/bin/perl

use strict;
use warnings;

use DB_File;

my $prefix = shift @ARGV || die;
my $put_id = shift @ARGV || '';

tie my %aln, 'DB_File', "$prefix.full_alignment.db", O_RDONLY;

my @put_ids = keys %aln;

#print scalar(@put_ids)."\n";

if ($put_id) {

    print $aln{$put_id};
} else {
    foreach $put_id(@put_ids) {
        print "## $put_id\n";
        print $aln{$put_id};
    }
}
