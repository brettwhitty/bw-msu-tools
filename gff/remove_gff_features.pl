#!/opt/rocks/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp qw{ tempfile };

my ($gff, $id_list, $remove_childless, $remove_parentless);

GetOptions(
    'gff|g=s'               =>  \$gff,
    'input|i=s'             =>  \$id_list,
    'remove_childless|c=s'  =>  \$remove_childless,
    'remove_parentless|p=s' =>  \$remove_parentless,
);

$remove_childless  ||= 'gene';
$remove_parentless ||= 'CDS,intron,exon';

my $feat_childless = {};
foreach my $feat_type(split(/,/, $remove_childless)) {
    $feat_childless->{$feat_type} = 1;
}

my $infh;

open $infh, '<', $id_list || die "Failed to open input list '$id_list' for reading: $!";

my $remove_ids = {};
while (<$infh>) {
    chomp;

    $remove_ids->{$_} = 1;
}

my $skip = {};

open $infh, '<', $gff || die "Failed to open GFF file '$gff' for reading: $!";
my ($tempfh, $tempfile) = tempfile(UNLINK => 1);

my $has_children = {};
while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    if (scalar(@t) != 9) {
        print $tempfh $_."\n";
        next;
    }

    my @atts = split(/;/, $t[8]);

    my $id = '';
    my @parents = ();
    my $noncoding = 0;
    my @new_atts;
    foreach my $att(@atts) {
        my ($key, $value_string) = split(/=/, $att);

        my @values = split(/,/, $value_string);

        if ($key eq 'ID') {
            $id = $values[0];
            if (defined($remove_ids->{$id})) {
                $skip->{$id} = 1;
                last;
            }
        } elsif ($key eq 'Parent') {
            @parents = @values;
            my @new_parents;
            foreach my $parent(@parents) {
                if (! defined($remove_ids->{$parent})) {
                    push(@new_parents, $parent);
                    $has_children->{$parent} = 1;
                }
            }
            @values = @new_parents;
            if (scalar(@new_parents) == 0) {
                $remove_ids->{$id} = 1;
                $skip->{$id} = 1;
            }
        }
        $value_string = join(",", @values);
        my $att_string = join("=", ($key, $value_string));
        push(@new_atts, $att_string);
    }
    $t[8] = join(";", @new_atts);

    if (! $skip->{$id}) {
        print $tempfh join("\t", @t)."\n";
    }
}
close $tempfh;


open $infh, '<', $tempfile || die "Failed to open temp GFF file '$tempfile' for reading: $!";

while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    if (scalar(@t) != 9) {
        print $_."\n";
        next;
    }

    $t[8] =~ /ID=([^;]+)/;
    my $id = $1;

    if ($feat_childless->{$t[2]} && defined($id) && ! defined($has_children->{$id})) {
        next;
    } else {
        print $_."\n";
    }
}
