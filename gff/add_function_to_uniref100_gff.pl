#!/usr/bin/perl

use strict;
use warnings;

use DB_File;

tie my %db, 'DB_File', '/scratch/whitty/uniref100.function.txt.db', O_RDONLY or die;

while (<>) {
    if (/Name=([^;]+);/) {
        my $id = $1;
        if (defined($db{$id})) {
            my $note = $db{$id};
            $note =~ s/([^A-Za-z0-9])/sprintf("%%%02X", ord($1))/seg;
            s/(Name=[^;]+;)/$1Note=$note;/;
        }
    }
    print $_;
}
