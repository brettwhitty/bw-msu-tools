#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use File::Basename;

my $dir = shift @ARGV || getcwd();

my @glob = glob($dir."/*");

foreach my $item(@glob) {
    unless (-d $item) {
        next;
    }
    $dir = basename($item);
        
    my ($tax_id, $name) = split(/\./, $dir, 2);

    $name =~ s/_/ /g;
    
    my @gbk_files = glob($item."/*.gbk");
    
    foreach my $file(@gbk_files) {
        my $filename = basename($file);
        $filename =~ /^(.*)\.gbk/;
        my $acc = $1;

        my $chrom = grep_gbk_file_for_chromosome($file);
        $chrom =~ s/plastid\://;
        $chrom =~ s/^[0]+//;
#        if ($chrom eq 'UNKNOWN') {
#            print STDERR $file."\n";
#        }

        print "$tax_id\t$name\t$chrom\t$acc\n";
    }
}

sub grep_gbk_file_for_chromosome {
    my ($file) = @_;

    open (IN, $file) || die "Can't open '$file' for reading: $!";

    my $chrom = 'UNKNOWN';
    
    while (<IN>) {
        if (/\/chromosome=\"([^\"]+)\"/) {
            $chrom = $1;
            last;
        } elsif (/\/organelle=\"([^\"]+)\"/) {
            $chrom = $1;
            last;
        }
    }
    
    return $chrom;
}
