#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    s/\/home\/whitty\/gmap\/(\d+)\.fa/gmap/;
    my $acc = $1;
    s/^NA/$acc/;
    print $_;
}
