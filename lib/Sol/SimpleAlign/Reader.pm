#!/usr/bin/perl

package Sol::SimpleAlign::Reader;

## Retrieves Bio::SimpleAlign objects that have been serialized
## using XML::Dumper and stored in a DB_File database
##
##
## Example usage:
##
## use Sol::SimpleAlignDB::Reader; 
##
## my $db = new Sol::SimpleAlignDB::Reader(
##                   db => "Solanum_tuberosum/Solanum_tuberosum.alignment.with_subseq.db"
##                                        ); 
##
## my $aln_obj = $db->get_alignment("PUT-157a-Solanum_tuberosum-68554"); 
##

use DB_File;
use XML::Dumper;

use Bio::AlignIO;
use Bio::SimpleAlign;

use Carp;

use Set::Scalar;

my $dumper = new XML::Dumper;

my %db;
my $ids;

## constructor requires path to database file
sub new {
   my ($class, %args) = @_;

    my $db_file = '';
    if (defined($args{'db'})) {
        $db_file = $args{'db'};
    } else {
        croak "Constructor requires parameter 'db'";
    }

    _init_db($db_file);

    return bless {}, ref($class) || $class;
}

## initialize the alignment db
sub _init_db {
    my ($db_file) = @_;

    tie %db, 'DB_File', $db_file, O_RDONLY or croak "Failed to tie '$db_file'";


    $ids = new Set::Scalar(keys %db);

}

## returns an alignment object from the database
sub get_alignment {
    my ($self, $id) = @_;
    
    if ($self->has_alignment($id)) {
        return $dumper->xml2pl($db{$id});
    } else {
        carp "get_alignment returned undef";
        return undef;
    }
}

## returns xml in case there is ever a need for that
sub get_alignment_xml {
    my ($self, $id) = @_;

    if ($self->has_alignment($id)) {
        return $db{$id};
    } else {
        carp "get_alignment_xml returned undef";
        return undef;
    }
}

## check if alignment is in DB
sub has_alignment {
    my ($self, $id) = @_;
    
    if ($ids->has($id)) {
        return 1;
    } else {
        return 0;
    }
}

1;
