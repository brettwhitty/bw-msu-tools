#!/usr/bin/perl

## Fetches XML for a UniRef record and returns the species name of the record
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;

use CGI qw{ :standard };
use CGI::Carp('fatalsToBrowser');
use LWP::Simple;
use XML::Twig;
use Carp;

my $cgi = new CGI;

my $id = $cgi->param('id');

$id =~ /UniRef\d+_\w+/ || croak "Invalid access.";

my $xml = get("http://www.uniprot.org/uniref/$id.xml");
if (! defined($xml)) {
    croak "Request failed.";
}

my $twig = new XML::Twig(
                            twig_handlers => {
                                                entry => \&entry_handler,
                                             }
                        );

$twig->parse($xml);

sub entry_handler {
    my ($twig, $entry) = @_;

    foreach my $prop($entry->children('property')) {
        if ($prop->att('type') eq 'common taxon') {
            print $prop->att('value');
        }
    }
    $twig->purge();
}
