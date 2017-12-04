#!/opt/rocks/bin/perl


## quick script to map SNP positions from PGSC contigs to supercontigs
##
## genome GFF is missing mapping for contigs where the contig is the supercontig
## see: /projects/potato/genome/PGSC0003DMB.scaffold.info


use strict;
use warnings;

use Getopt::Long;

my ($gff, $snp_table);

GetOptions(
    'gff|g=s'       =>  \$gff,
    'snp_table|s=s' =>  \$snp_table,
);

my $infh;
my $snps = [];

open $infh, '<', $snp_table || die "$!";
while (<$infh>) {
    chomp;
    
    my @t = split(/\t/, $_);
    my ($snp_id, $contig_id, $position) = @t;

    push(@{$snps}, [$snp_id, $contig_id, $position]);
}
open $infh, '<', $gff || die "$!";

my $contig_pos = {};
while (<$infh>) {
    chomp;

    my @t = split(/\t/, $_);

    if (scalar(@t) != 9) {
        next;
    }

    if ($t[2] ne 'contig') {
        next;
    }

    my ($seqid, $start, $end, $strand, $atts) = @t[0, 3, 4, 6, 8];
    $atts =~ /ID=([^;]+)/;
    my $contig_id = $1;
    
    $contig_pos->{$contig_id} = {
        'supercontig'   =>  $seqid,
        'start'         =>  $start,
        'end'           =>  $end,
        'strand'        =>  $strand,
    };
}
foreach my $snp_ref(@{$snps}) {

    my ($snp_id, $contig_id, $position) = @{$snp_ref};
    
    my $contig_ref = $contig_pos->{$contig_id};
    if (! defined($contig_ref)) {
        die "$contig_id";
    }
    my $map_pos;
    if ($contig_ref->{'strand'} eq '+') {
        $map_pos = $position + ($contig_ref->{'start'} - 1);
    } else {
        $map_pos = ($contig_ref->{'end'} - 1) - $position;
    }
    print join ("\t", (
            $snp_id,
            $contig_ref->{'supercontig'},
            $map_pos,
    ))."\n";
}
