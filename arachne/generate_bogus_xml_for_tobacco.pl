#!/usr/bin/perl

use Carp;
use strict;
use warnings;
use Bio::DB::Fasta;
use Cwd qw{ abs_path };

my $fname = shift @ARGV;

$fname = abs_path($fname);

my $out_fname = $fname.".xml";

open (OUT, ">$out_fname") || confess "Failed to open '$out_fname' for writing";

my $db = new Bio::DB::Fasta($fname);
my @ids = $db->get_all_ids();

print OUT <<XML;
<?xml version="1.0"?>
<trace_volume>
XML

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
    
print OUT <<XML;
    <trace>
       <trace_name>$gb_id</trace_name>
       <type>paired_production</type>
       <library_id>$lib_id</library_id>
       <plate_id>$plate_id</plate_id>
       <well_id>$well_id</well_id>
       <template_id>$template_id</template_id>
       <trace_end>$direction</trace_end>
       <clip_vector_left>1</clip_vector_left>
       <clip_vector_right>$seq_len</clip_vector_right>
       <insert_size>2500</insert_size>
       <insert_stdev>500</insert_stdev>
    </trace>
XML

}
print OUT "</trace_volume>\n";
