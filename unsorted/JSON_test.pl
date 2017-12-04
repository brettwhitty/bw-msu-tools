#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $test = {
        'key1'  =>  'value1',
        'key2'  =>  'value2',
        'key3'  =>  ['a', 'b'],
    };

print encode_json($test);
