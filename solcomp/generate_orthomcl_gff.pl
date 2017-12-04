#!/usr/bin/perl

$| = 1;

use strict;
use warnings;


## this script is meant to take gene models produced by MAKER
## which have been assigned to an orthomcl cluster in some way
## and produce GFF3 output suitable for loading to gbrowse
##
## it requires a table of gene model -> orthomcl mapping
## and the GFF3 output by maker that encodes the gene models
## and optionally a table of orthomcl id -> functional text
##
## it can optionally work on features other than gene models
##
## Brett Whitty
## whitty@msu.edu


use lib '/home/whitty/SVN/lib';
use MyIO;
use GFFTextEncoder;

use Getopt::Long;
use File::Find::Rule;
use Carp;

## store cluster memberships
my $cluster = {};
## for adding Note
my $cluster_function = undef;
## for adding sol_members
my $sol_members = undef;
## for generating IDs
my $id_counter = {};

## args
my $orthomcl_table;         ## model id -> orthomcl mapping table
my $gff_dir;                ## dir containing GFF files
my $suffix = '.gff';        ## suffix of GFF files
my $type = 'mRNA';          ## default GFF3 feature type
my $source = 'orthomcl';    ## source for output GFF3
my $feature = 'match';      ## feature for output GFF3
my $function_table = '';    ## orthomcl id -> function text table
my $function_unknown = 'uncharacterized protein';
my $sol_members_table = '';  ## table file listing sol member species in cluster
my $id_transform = '';

my $result = GetOptions(
    'input|i=s'     =>  \$orthomcl_table,
    'gff_dir|g=s'   =>  \$gff_dir,
    'suffix|x=s'    =>  \$suffix,
    'type|t=s'      =>  \$type,
    'source|s=s'    =>  \$source,
    'feature|f=s'   =>  \$feature,
    'function=s'    =>  \$function_table,
    'unknown=s'     =>  \$function_unknown,
    'sol_members=s' =>  \$sol_members_table,
    'id_transform=s'=>  \$id_transform,
);

unless ($orthomcl_table && -e $orthomcl_table) {
    confess "You must specify an orthomcl mapping table with --input";
}
unless ($gff_dir && -d $gff_dir) {
    confess "You must specify a directory containing GFF files with --gff_dir";
}

my $infh;

## read in the orthomcl assignments from simple blast table
$infh = get_infh($orthomcl_table);
while (<$infh>) {
    chomp;
    my @t = split("\t", $_);

    $cluster->{$t[0]} = $t[1];
}

## read in the cluster function table if it has been provided
if ($function_table && -e $function_table) {
    $infh = get_infh($function_table);
    while (<$infh>) {
        chomp;
        my @t = split("\t", $_);

        $cluster_function->{$t[0]} = $t[1];
    }
}

## read in the sol members table if it has been provided
if ($sol_members_table && -e $sol_members_table) {
    $infh = get_infh($sol_members_table);
    while (<$infh>) {
        chomp;
        my @t = split("\t", $_);

        my @members = split(/,/, $t[1]);
        $sol_members->{$t[0]} = scalar(@members);
    }
}

my @gff_files = File::Find::Rule->file()
                                ->name('*'.$suffix)
                                ->in($gff_dir);

## iterate through all GFF files and find features we recognize
foreach my $gff_file(@gff_files) {
    $infh = get_infh($gff_file);

    while (<$infh>) {

        ## find the feature type we're looking for
        if (/\t$type\t/) {

            ## match on the feature's ID
            my ($feature_id) = /\tID=([^;]+)/;

            ##SHOUD ALLOW ARBITRARY ID TRANSFORMATION HERE USING THE FLAG

            ## for grape
            if ($feature_id =~ /GSVIVT/) {
                $feature_id =~ s/GSVIVT/GSVIVP/;
            }

            ## only deal with things we have a cluster assignment for
            if (defined($cluster->{$feature_id})) {
                chomp;

                my $orthomcl_id = $cluster->{$feature_id};

                my @t = split("\t", $_);
                $t[1] = $source;
                $t[2] = $feature;
                my $id = join('-', (
                        $t[0],        ## sequence id
                        $orthomcl_id, ## orthomcl id
                        $id_counter->{$t[0].$orthomcl_id}++, ## counter
                ));

                ## prepare the GFF3 feature attributes
                my %attribs = (
                    'ID'    =>  $id,
                    'Name'  =>  $orthomcl_id,
                );

                ## add a note attribute if the user has specified a function table
                if ($cluster_function) {
                    my $note = defined($cluster_function->{$orthomcl_id}) ?
                        $cluster_function->{$orthomcl_id} : $function_unknown;
                    $attribs{'Note'} = $note;
                }
                ## add a note attribute if the user has specified a function table
                if ($sol_members) {
                    my $count = defined($sol_members->{$orthomcl_id}) ?
                        $sol_members->{$orthomcl_id} : 0;
                    $attribs{'sol_members'} = $count;
                }


                ## GFF3 formatting of attributes
                my @attribs = ();
                foreach my $key(sort keys(%attribs)) {
                    push(@attribs, $key.'='.encode($attribs{$key}) );
                }
                $t[8] = join(';', @attribs);

                print join("\t", @t)."\n";
            }
        } else {
            next;
        }
    }
}
