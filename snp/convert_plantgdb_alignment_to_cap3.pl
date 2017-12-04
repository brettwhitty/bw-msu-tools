#!/usr/bin/perl

use strict;
use warnings;

use DB_File;
use File::Basename qw{ basename dirname };
use Cwd qw{ abs_path };

my $infile = shift @ARGV || die;
my $outfile = shift @ARGV;

$infile = abs_path($infile);

unless ($outfile) {
    $outfile = dirname($infile).'/'.basename($infile, ".txt").'.cap3';
}

open(IN, $infile) || die "$!";
open(OUT, ">$outfile") || die "$!";

print STDERR "Converting PlantGDB file '".basename($infile)."' to CAP3 format '".basename($outfile)."'\n";

my $members = {};

my $msa_acc;
my $align;
while (<IN>) {
    chomp;
    if (/^>(\S+)/) {
        $msa_acc = $1;
    } elsif (/^\S+/) {
        if (/^consensus/) {
            next;
        }

        /^(.{22})(.*)/ || die "Failed to match sequence line properly";
        my ($acc, $seq) = ($1, $2);

        $acc =~ s/\s+//g;
        $members->{$msa_acc}->{$acc} = 1;
    }
}

print_fake_cap3_header();
dump_put_members();
print OUT "\nDETAILED DISPLAY OF CONTIGS\n";
dump_alignments();
print OUT "\n";

sub dump_put_members {
    
    my @put_ids = keys(%{$members});   

    foreach my $put_id(@put_ids) {
        print_contig_header($put_id);
        foreach my $member(keys(%{$members->{$put_id}})) {
            print OUT "$member\n";
        }
    }
}

sub dump_alignments {
    open (IN, $infile) || die "$!";
    while (<IN>) {
        if (/^>(\S+)/) {
            $msa_acc = $1;
            print_contig_header($msa_acc);
        } else {
            print OUT $_;
        }
    }
        
}

sub print_fake_cap3_header {
    print OUT <<END;
Number of segment pairs = 12345; number of pairwise comparisons = 6789
'+' means given segment; '-' means reverse complement

Overlaps            Containments  No. of Constraints Supporting Overlap    

END
}

sub print_contig_header {
    my ($id) = @_;

    print OUT <<END;
******************* $id ********************
END
}

sub print_dot_line {
    print OUT <<END;
                          .    :    .    :    .    :    .    :    .    :    .    :
END
}

sub print_bottom_line {
    print OUT <<END;
                      ____________________________________________________________
END
}

sub print_seq_line {
    my ($id, $seq) = @_;

    unless ($seq =~ /^\s+$/) {
        print OUT $id.(' ' x (22 - length($id)))."$seq\n";
    }
}
