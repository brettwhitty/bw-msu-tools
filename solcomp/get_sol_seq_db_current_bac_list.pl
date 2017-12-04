#!/usr/bin/perl

$| = 1;

use lib "/home/whitty/SVN/lib";
use Sol::SeqDB;

my $db = new Sol::SeqDB();

my @gis = $db->get_gi_list('contig');

foreach my $gi(@gis) {
    my $atts = $db->get_gb_atts($gi);
    use Data::Dumper;
    print join("\t", (
            $atts->{'taxon_id'},
            $atts->{'replicon'},
            $atts->{'division'},
            $atts->{'date'},
            $gi,
            $atts->{'accession'}.'.'.$atts->{'version'},
            $atts->{'clone'},
            $atts->{'len'},
    ))."\n";
}
