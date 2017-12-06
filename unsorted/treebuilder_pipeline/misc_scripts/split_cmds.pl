#!/usr/bin/perl

use strict;
use warnings;

my $counter = 0;
my $outfh; 
my $lc = 15;
while(<>) {
        if ($lc == 15) { 
            $lc = 0;
            open($outfh, ">cmds/cmds.".++$counter.".sh");
        }
        print $outfh $_; 
        $lc++; 
}
