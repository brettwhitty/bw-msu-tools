#!/usr/bin/perl

$| = 1;
## make sure that set_gb_obsolete_not_current has been run after last load of new sequences


use lib "/home/whitty/SVN/lib"; 
use Sol::SeqDB; 
use File::Path qw{ mkpath };

my $db = new Sol::SeqDB(); 

my $tax = $db->get_taxonomy_hash();
my @list = $db->get_gi_list("read"); 

foreach my $gi(@list) {
    my $att = $db->get_gb_atts($gi);
    
    my $out_dir = $att->{'taxon_id'}
                  .'.'
                  . $tax->{$att->{'taxon_id'}}->{'scientific_name'};
    $out_dir =~ s/\s/_/g;
    my $out_file = $att->{'taxon_id'}
                  .'.'
                  . $tax->{$att->{'taxon_id'}}->{'scientific_name'}
                   . ".fna";

#    print "$out_dir/$out_file\t$att->{replicon}\n";

    unless (-e $out_dir) {
        mkpath($out_dir);
    }

    open ($out_fh, ">>$out_dir/$out_file") || die "Failed to open '$out_dir/$out_file' for writing: $!";
    my $fasta = $db->get_gb_fasta($gi) or die "Failed fetching fasta for '$gi'";
    
    $fasta =~ s/>gi\|\d+\|gb\|([^.]+)[^\|]+\|/>$1/g;

    print $out_fh $fasta;
    close $out_fh;
                  
}
