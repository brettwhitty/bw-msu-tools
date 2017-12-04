#!/usr/bin/perl

## converts the output from the script ssr.pl (S.Cartinhour 5/2000)
## to gff3 format

use warnings;
use strict;

use Getopt::Long;

my $input;
my $output;
my $method = 'SSR_putative';
my $feature = 'repeat_region';
my $add_prefix = 0;

my $result = GetOptions(
                            'input|i=s'     =>  \$input,
                            'output|o=s'    =>  \$output,
                            'method|m=s'    =>  \$method,
                            'feature|f=s'   =>  \$feature,
                            'add_prefix|p!' =>  \$add_prefix,
                       );

my %id_counter;

my $infh = get_infh($input);
my $outfh = get_outfh($output);

while (<$infh>) {
    chomp;
    my @line = split(/\t/);
    my $motif = "($line[3])$line[4]";
    unless ($line[0] && $line[5] && $line[6] && $motif){die "Error\n";}

    my $ssr_id;
    if ($add_prefix) {
        $id_counter{$line[0]}++;
        my $idpre = $line[0];
        if ($idpre =~ /^PGSC/) {
            $idpre =~ s/([A-Z]+)[0]+/$1/g;
        }
        $ssr_id = $idpre.".ssr.".$id_counter{$line[0]};
    } else {
        $id_counter{''}++;
        $ssr_id = "ssr.".$id_counter{''};
    }
    print $outfh "$line[0]\t$method\t$feature\t$line[5]\t$line[6]\t.\t+\t.\tID=$ssr_id;Name=$motif;motif=$line[3];repeat_count=$line[4]\n";
}

sub get_infh {
    my ($in) = @_;

    unless ($in) {
        return \*STDIN;
    }
    open (my $infh, $in) || die "Couldn't open input file '$in' for reading: $!";
    
    return $infh;
}   

sub get_outfh {
    my ($out) = @_;

    unless ($out) {
        return \*STDOUT;
    }
    open (my $outfh, ">$out") || die "Couldn't open output file '$out' for writing: $!";
    
    return $outfh;
}    
