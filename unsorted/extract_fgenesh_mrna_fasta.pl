#!/opt/rocks/bin/perl

use strict;
use warnings;

## quick script to pull out FASTA from Fgenesh runs with -pmrna option selected
## 
## Brett Whitty
## whitty@msu.edu

use Getopt::Long;
use Carp;

my ($na, $aa);

GetOptions(
                'a|p!'  =>  \$aa,
                'n!'    =>  \$na,
          );

if (! ($na || $aa)) {
    confess "Must specify either -a or -n flag for aa or na fasta";
}

my $seq_id = '';
my $flag = 0;
while (<>) {
    chomp;
    if (/ Seq name: (\S+)/) {
        $seq_id = $1;
    }
    if (/>(FGENESH:(\[mRNA\])?\s+(\d+).*)/) {
        my ($defline, $type, $count) = ($1, $2, $3);
        if (! defined($type)) {
            $type = '[pep]';
        }
        if ($aa && $type eq '[pep]') {
            $flag = 1;
        } elsif ($na && $type eq '[mRNA]') {
            $flag = 1;
        } else {
            $flag = 0;
        }
        #print "$type\t$count\n";
        $_ = ">$seq_id.fgenesh.mRNA.$count $defline";
    }  
    if ($flag) {
        print $_."\n";
    }
}
