#!/usr/bin/perl

package GFFTextEncoder;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(encode);  # symbols to export

use strict;
use warnings;

sub encode {
    my ($value) = @_;

    $value =~ s/([\t\n\r%&\=;,])/sprintf("%%%02X",ord($1))/ge;

    return $value;
}

1;
