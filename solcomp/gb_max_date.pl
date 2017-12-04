#!/opt/rocks/bin/perl

use strict;
use warnings;

use lib "/home/whitty/SVN/lib"; 
use Sol::SeqDB; 

my $type = shift @ARGV || die "provide type 'contig', 'read', etc.";

my $db = new Sol::SeqDB; 

my $d = $db->get_gb_max_date($type); 
print $d->format("%Y/%m/%d");
