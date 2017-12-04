#!/usr/bin/perl

package MyIO;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_infh get_outfh);  # symbols to export
  
sub get_infh {
    my ($in) = @_;

    unless ($in) {
        return \*STDIN;
    }
    my $infh;
    if ($in =~ /.gz$/) {
        open $infh, "<:gzip", $in || die "Couldn't open input file '$in' for reading: $!";
    } else {
        open ($infh, $in) || die "Couldn't open input file '$in' for reading: $!";
    }

    return $infh;
}

sub get_outfh {
    my ($out) = @_;

    unless ($out) {
        return \*STDOUT;
    }
    open (my $outfh, ">".$out) || die "Couldn't open output file '$out' for writing: $!";

    return $outfh;
}

1;
