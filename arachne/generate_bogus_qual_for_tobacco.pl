#!/usr/bin/perl

use Carp;
use strict;
use warnings;
use Bio::DB::Fasta;
use Cwd qw{ abs_path };
use File::Temp;

my $fname = shift @ARGV;

$fname = abs_path($fname);

my $out_fname = $fname.".qual";

#open (OUT, ">$out_fname") || confess "Failed to open '$out_fname' for writing";

my $db = new Bio::DB::Fasta($fname);
my @ids = $db->get_all_ids();

foreach my $id(@ids) {
   
    my $seq_len = $db->length($id);
    my $defline = $db->header($id);
   
    my ($gb_id, $trace_name, $lib_id);
    
    if ($defline =~ /^(\S+)\s+(\S+)\s+(\S+)/) {
        ($gb_id, $trace_name, $lib_id) = ($1, $2, $3);
    } else {
        confess "Failed to match defline '$defline'";
    }

    $trace_name =~ /^((.*)(\d{3})x([a-z]\d+))([a-z])\d+\.ab1/ || confess "Failed to match on trace name '$trace_name'";
    my $template_id = $1;
#    my $plate_id = $2;
    my $plate_id = $2.$3;
    my $well_id = $4;
    $well_id =~ tr/a-z/A-Z/;
    my $direction = $5;
    $direction =~ tr/a-z/A-Z/;
    
    print ">$defline\n";
    
    my @quals = ('15') x $seq_len;

    while (my @q_line = splice(@quals, 0, 19)) {
        print join(" ", @q_line)."\n";
    }
}
