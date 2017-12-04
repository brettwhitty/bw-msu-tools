#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    if (/ID=([^;^\s]+)/) {
        my $id = $1;
        s/(ID=[^;^\s]+)/$1;Name=$id/;
    }
    print;
}
