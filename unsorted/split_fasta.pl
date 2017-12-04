#!/usr/bin/perl

# test

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename;
use POSIX qw(ceil);

my ($file_count, $chunk_size, $prefix, $input, $output_dir, $verbose, $header, $bases_per);

my $result = GetOptions(
             "input|i=s"        =>  \$input,
             "output_dir|o=s"   =>  \$output_dir,
             "chunk_size|c=i"   =>  \$chunk_size,
             "file_count|n=i"   =>  \$file_count,
             "bases_per|b=i"    =>  \$bases_per,
             "prefix|p=s"       =>  \$prefix,
             "verbose|v!"       =>  \$verbose,
             "header|h!"        =>  \$header, ## use header for file name
               );

unless ($chunk_size || $file_count || $bases_per) {
    confess "Must provide value to either --chunk_size|c, --file_count|n or --bases_per|b flags";
}
if ($output_dir && ! -d $output_dir) {
    confess "Specified output dir '$output_dir' doesn't exist";
}
unless ($output_dir) {
    $output_dir = dirname($input);
}

my $suffix = '';
if (! $prefix) {
    if ($input) {
        my $basename = basename($input);
        $basename =~ s/(\.[^\.]+(\.masked)?)$//; ## remove file suffix
        ## save the suffix for later
        $suffix = $1;
        unless ($basename) {    ## in case nothing's left
            $prefix = "split_fasta";    
        }
        $prefix = $basename;
    } else {
        $prefix = "split_fasta";
    }
}

if ($file_count) {
    my $infh = get_input_fh($input);
    my $header_count = header_count($infh);
    $chunk_size = ceil($header_count / $file_count);
}


if ($chunk_size || $file_count) {
    do_n_file_chunk();
} elsif ($bases_per) {
    do_bases_per();
}


sub do_n_file_chunk {
    my $infh = get_input_fh($input);

    my $file_id = 0;
    my $header_count = 0;
    my $outfile;
    my $outfh;
    while (<$infh>) {
        if (/^>(\S+)/) {
            my $acc = $1;
            $header_count++;
            if (($header_count - 1) % $chunk_size == 0) {
                if ($header) {
                    $outfile = $output_dir."/".$acc.$suffix;
                } else {
                    $outfile = $output_dir."/".$prefix."_".++$file_id.$suffix;
                }
                if ($verbose) {
                    print abs_path($outfile)."\n";
                }
                $outfh = get_output_fh($outfile);
            }
        }
        print $outfh $_;
    }
}

sub do_bases_per {
    my $infh = get_input_fh($input);

    my $file_id = 0;
    my $base_count = -1;
    my $header_count = 0;
    my $outfile;
    my $outfh;
    while (<$infh>) {
        if (/^>(\S+)/) {
            my $acc = $1;
            $header_count++;
            if ($base_count == -1 || $base_count > $bases_per) {
                $base_count = 0;
                if ($header) {
                    $outfile = $output_dir."/".$acc.$suffix;
                } else {
                    $outfile = $output_dir."/".$prefix."_".++$file_id.$suffix;
                }
                if ($verbose) {
                    print abs_path($outfile)."\n";
                }
                $outfh = get_output_fh($outfile);
            }
        } else {
            my $line = $_;
            $line =~ s/\s+//g;
            $base_count += length($line);
        }
        print $outfh $_;
    }

}

sub get_output_fh {
    my ($filename) = @_;

    open (my $outfh, ">".$filename) || confess "Can't open file '$filename' for writing: $!";

    return $outfh;
}

sub header_count {
    my ($infh) = @_;

    my $count = 0;
    while (<$infh>) {
        if (/^>/) {
            $count++;
        } else {
            next;
        }
    }

    return $count;
}


sub get_input_fh {
    
    my ($input) = @_;

    my $infh;

    if ($input) {
        open ($infh, $input) || confess "Failed opening '$input' for reading: $!";
    } else {
        open ($infh, \*STDIN) || confess "Failed opening STDIN for reading: $!";
    }

    return $infh;
}
