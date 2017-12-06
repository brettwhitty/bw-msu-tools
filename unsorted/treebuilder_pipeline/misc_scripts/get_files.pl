#!/usr/bin/perl

use File::Basename; 
use File::Copy; 

while (<>) { 
    chomp; 
    $base = basename($_, ".xml"); 
    $dir = dirname($_); 
    unless (-d $base) { 
        mkdir($base);
    }
    my @exts = (
        '.codes.txt',
        '.fasta',
        '.msa',
        '.gfsa',
        '.gfsa.trim',
        '.gfsa.trim.txt',
        '.newick',
        '.final.newick',
        '.phyi',
        '.phyi.trim',
        '.xml',
        '.svg',
        '.png',
        '.inc.html',
        '.map',
    );
    foreach my $ext(@exts) { 
        copy($dir."/".$base.$ext, $base."/");
    }
}
