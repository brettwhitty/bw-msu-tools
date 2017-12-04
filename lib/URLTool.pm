#!/usr/bin/perl

package URLTool;

##
## URLTool.pm
##
## Right now just has one method called key_replace for doing search and replace
## on placeholders in URLs
##
## eg:
##
## my $url = key_replace('http://some_url/script.cgi?tax_id={taxon_id}&type={some_type}', {
##      taxon_id    =>  '4081',
##      some_type   =>  'a_type',
## });
##
## This is to use for storing generic URLs in a DB (or elsewhere) and then subsequently
## instantiating real URLs with the values populated programmatically.
##
## Brett Whitty
## whitty@msu.edu
##

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(key_replace);  # symbols to export

use strict;
use warnings;
use Carp;

sub key_replace {
    my ($url, $args, $open, $close) = @_;

    unless(defined($open)) {
        $open = '{';
    }
    unless(defined($close)) {
        $close = '}';
    }

    while ($url =~ /\Q$open\E(.*?)\Q$close\E/gs) {
        my ($m_start, $m_end) = ($-[0], $+[0]);
        my $key = $1;
        if (defined($args->{$key})) {
            $url = substr($url, 0, $m_start) . $args->{$key} . substr($url, $m_end);
        } else {
            carp "Unmatched placeholder '$key' in:\n$url";
        }
    }

    return $url;
}

1;
