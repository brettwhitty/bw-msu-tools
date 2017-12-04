#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
    if (/^cigar:/) {
        / \S+:subseq\((\d+),\d+\) (\d+) (\d+)/;
        my ($offset, $start, $end) = ($1, $2, $3);
        $start = $start + $offset;
        $end = $end + $offset;
        s/ (\S+):subseq\((\d+),\d+\) (\d+) (\d+)/ $1 $start $end/;
        print;
    } elsif (scalar(@_ = split("\t")) == 9) {
        /^(\S+):subseq\((\d+),\d+\)\t/;
        my ($acc, $offset) = ($1, $2);
        
        my @t = split("\t");
        
        $t[0] = $acc;
        $t[3] += $offset;
        $t[4] += $offset;
   
        my @a = split(" ; ", $t[8]);
        my @new = ();
        foreach my $a(@a) {
            if ($a =~ /^Align (\d+)/) {
                my $pos = $1;
                $pos += $offset;
                $a =~ s/Align \d+/Align $pos/;
            }
            push(@new, $a);
        }
        $t[8] = join(" ; ", @new);
        
        print join("\t", @t);
    } elsif (scalar(@_ = split("\t")) == 8) {
        /^(\S+):subseq\((\d+),\d+\)\t/;
        my ($acc, $offset) = ($1, $2);
        s/^\S+/$acc/;
        print;
    } else {
        print;
    }
    # (/^\S+:subseq\((\d+),\d+\)\t/;
        
}
