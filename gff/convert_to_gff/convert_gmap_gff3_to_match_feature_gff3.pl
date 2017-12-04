#!/usr/bin/perl

use warnings;
use strict;
use Carp;
use Getopt::Long;

my $input;
my $output;
my $min_pct_id = 0;

my $result = GetOptions(
                            'input|i=s'   => \$input,
                            'output|o=s'  => \$output,
                            'pct_id|p=i'  => \$min_pct_id,
                       );

unless ($input) {
    confess "Must provide input file with --input flag";
}

## filehandles
my ($infh, $outfh);

## will store match feature positions
my $matches = {};

## used to calculate weighted averages for the match features
my $avg_sum = {};
my $len_sum = {};

## do a pass through the input file to find match feature coordinates
$infh = open_read($input);    
while (<$infh>) {
    if (/^#/) { next;}
    
    chomp;

    my @col = split(/\t/);

    my $acc     = $col[0];
    my $start   = $col[3];
    my $end     = $col[4];
    my $pct_id  = $col[5];
    my $strand  = $col[6];
    my $att     = $col[8];
    $att =~ /ID=([^;]+)/ || confess "Couldn't parse ID attribute out of GFF3 file!";
    my $id      = $1;

    ## weighted average calculation
    my $len = $end - $start + 1;
    $len_sum->{$acc}->{$id} += $len;
    $avg_sum->{$acc}->{$id} += $pct_id * $len;
    
    ## set the strand of the match unless defined already
    unless (defined($matches->{$acc}->{$id}->[2])) {
        $matches->{$acc}->{$id}->[2] = $strand;
    }
    
    ## set min match start
    if (! defined($matches->{$acc}->{$id}->[0]) || $matches->{$acc}->{$id}->[0] > $start) {
        $matches->{$acc}->{$id}->[0] = $start;
    }
    ## set max match end
    if (! defined($matches->{$acc}->{$id}->[1]) || $matches->{$acc}->{$id}->[1] < $end) {
        $matches->{$acc}->{$id}->[1] = $end;
    }
}



## flag hash to determine whether match features have been written to output
my $match_ids = {};
my $match_counter = 0;

$infh  = open_read($input);
$outfh = open_write($output);
while (<$infh>) {
    if (/^#/) { print $outfh $_; next;}
    
    chomp;

    my @col = split(/\t/);

    my $acc     = $col[0];
    my $att     = $col[8];
    $att =~ /ID=([^;]+)/ || confess "Couldn't parse ID attribute out of GFF3 file!";
    my $id      = $1;
 
    unless (defined($match_ids->{$acc}->{$id})) {
        
        ## calculate percent identity for the match feature
        my $pct_id = sprintf("%.f", $avg_sum->{$acc}->{$id} / $len_sum->{$acc}->{$id});
       
        if ($pct_id < $min_pct_id) {
            next;
        }
        
        ## parse out the display name
        $att =~ /Name=([^;]+)/ || confess "Couldn't parse Name attribute out of GFF3 file!";
        my $name = $1;

        ## generate an ID for the match feature
        my $match_id = 'match.'.++$match_counter;

        print $outfh join("\t", (
                $col[0],
                $col[1],
                'match',
                $matches->{$acc}->{$id}->[0],
                $matches->{$acc}->{$id}->[1],
                $pct_id,
                $matches->{$acc}->{$id}->[2],
                '.',
                "ID=$match_id;Name=$name",
                                ))."\n";
        $match_ids->{$acc}->{$id} = $match_id;                            
    }

    ## skip outputting HSPs unless we've created a match feature
    ## (which won't happen if we've filtered out a match based on avg % identity)
    unless (defined($match_ids->{$acc}->{$id})) {
        next;
    }
    
    ## replace cDNA_match with HSP
    $_ =~ s/cDNA_match/HSP/ || confess "Failed to replace cDNA_match with HSP!";

    ## add in the Parent attribute
    my $match_id = $match_ids->{$acc}->{$id};
    $_ =~ s/^(.*ID=[^;]+;)(.*)$/$1Parent=$match_id;$2/ || confess "Failed to add Parent attribute!";

    print $_."\n";
}


## do a second pass through the file to add the match features
sub open_read {
    my ($in) = @_;

    if ($in) {
        unless (-f $in) {
            confess "Specified input file '$in' doesn't exist";
        }
        open (my $fh, "$in") || confess "Failed opening '$in' for reading";
        return $fh;
    }
}    

sub open_write {
    my ($out) = @_;

    if ($out) {
        open (my $fh, ">$out") || confess "Failed opening '$out' for writing";
        return $fh;
    } else {
        return \*STDOUT;
    }
}    
