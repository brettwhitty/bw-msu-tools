#!/usr/bin/perl

## quick script to parse GenBank dbEST est.full.YYYYMMDD.bcp.x dump files to give GI -> Library ID mappings
##
## see ftp://ftp.ncbi.nih.gov/repository/dbEST/bcp/

use strict;
use warnings;

while (<>) {

    chomp;

    my ($dbest_id, $gi, $_id2, $_id3, $date, $library_id, $_id4, $_id5, $est_name, $gb_acc, $__and_the_rest)
        = split(/\t/);

    print "$gi\t$library_id\n";
}
