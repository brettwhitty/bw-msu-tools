#!/usr/bin/perl

use strict;
use warnings;

## this is a very specific purpose script that will take a .list file
## generated by my split_genbank_flat_files.pl script and use that
## to merge individual files in a directory into a single file based
## on the replicon they are from (4th column in .list file)

use Cwd qw( abs_path );
use DateTime;
use File::Basename;
use File::Copy;
use File::Temp qw( tempfile );
use Roman;

my $months = {
                'JAN' => 1,
                'FEB' => 2,
                'MAR' => 3,
                'APR' => 4,
                'MAY' => 5,
                'JUN' => 6,
                'JUL' => 7,
                'AUG' => 8,
                'SEP' => 9,
                'OCT' => 10,
                'NOV' => 11,
                'DEC' => 12
             };             

my $target_dir = shift @ARGV or die;

$target_dir = abs_path($target_dir);

my $base_dir = basename($target_dir);
my $parent_dir = basename(dirname($target_dir));

my $list_file = $target_dir."/.list";

my $timestamps = {};
my $accs = {};

## scan through to identify duplicate sequence files based on checksum (keep the newest)
open (my $list_fh, $list_file) || die "Failed to open list file '$list_file' for reading: $!";
while (<$list_fh>) {
    chomp;
    
    my @cols = split(/\t/);
    
    my $mol = $cols[3];
    my $crc = $cols[7];
    my $acc = $cols[2];

    my ($day, $month, $year) = split(/-/, $cols[6]);
    
    my $dt = new DateTime( day => $day, month => $months->{$month}, year => $year );

    if (defined($timestamps->{$crc})) {
        if (DateTime->compare($timestamps->{$crc}, $dt)) {
            print STDERR "current '$crc' is newer\n";
            next;
        }
    }
    $timestamps->{$crc} = $dt;
    $accs->{$crc} = $acc;
}
my $keep = {};

## this hash will determine which accessions we keep
foreach my $acc(values(%{$accs})) {
    $keep->{$acc} = 1;
}

## this will store the replicon mappings of the accessions we're keeping
my $mols = {};

## we'll write a new list file while removing the duplicate sequences
my ($out_fh, $out_file) = tempfile(CLEANUP => 1, OPEN => 1);

## scan through the list again to populate the hash of replicon mappings for files we're keeping
open ($list_fh, $list_file) || die "Failed to open list file '$list_file' for reading: $!";
while (<$list_fh>) {
    chomp;

    my @cols = split(/\t/);

    my $mol = $cols[3];
    my $crc = $cols[7];
    my $acc = $cols[2];
    
    ## convert from roman numerals if we have to
    if (isroman($mol)) {
        $mol = arabic($mol);
    }
    
    if (defined($keep->{$acc})) {
        print $out_fh $_."\n"; 
    }

    push(@{$mols->{$mol}}, $acc)
}
close $out_fh;
close $list_fh;

## copy over new list file
copy($out_file, $list_file);


## will create files that contain all molecules as well
my $all_out_prefix = "$target_dir/$base_dir.all";
my $all_gbk_out = $all_out_prefix.".gbk";
my $all_fsa_out = $all_out_prefix.".fsa";
open (my $all_gbk_fh, ">$all_gbk_out") || die "Failed to open '$all_gbk_out' for writing: $!";
open (my $all_fsa_fh, ">$all_fsa_out") || die "Failed to open '$all_fsa_out' for writing: $!";


## concatenate genbank and fasta files to new output files
foreach my $mol(sort {$a cmp $b} keys(%{$mols}) ) {

    my $out_prefix = "$target_dir/$base_dir.$mol";
    my $gbk_out = $out_prefix.".gbk";
    my $fsa_out = $out_prefix.".fsa";
    
    open (my $gbk_fh, ">$gbk_out") || die "Failed to open '$gbk_out' for writing: $!";
    open (my $fsa_fh, ">$fsa_out") || die "Failed to open '$fsa_out' for writing: $!";
    
    foreach my $acc(sort(@{$mols->{$mol}})) {
       
        my $in_fh;
        my $in_gbk = "$target_dir/$acc.gbk";
        open ($in_fh, $in_gbk) || die "Failed to open '$in_gbk' for reading: $!";
        while (<$in_fh>) {
            print $gbk_fh $_;
            print $all_gbk_fh $_;
        }
        close $in_fh;
        unlink($in_gbk);
        
        my $in_fsa = "$target_dir/$acc.fsa";
        open ($in_fh, $in_fsa) || die "Failed to open '$in_fsa' for reading: $!";
        while (<$in_fh>) {
            print $fsa_fh $_;
            print $all_fsa_fh $_;
        }
        close $in_fh;
        unlink($in_fsa);

    }
}
