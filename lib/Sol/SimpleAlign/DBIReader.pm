#!/usr/bin/perl

package Sol::SimpleAlign::DBIReader;

## Retrieves Bio::SimpleAlign objects that have been serialized
## using XML::Dumper and stored in a DBI database
##
##
## Example usage:
##
## use Sol::SimpleAlignDB::DBIReader; 
##
## my $db = new Sol::SimpleAlignDB::DBIReader(
##                   db         =>  'dbi:Pg:dbname=solcomp_search;host=pg',
##                   user       =>  'access',
##                   password   =>  'access'
##                                        ); 
##
## my $aln_obj = $db->get_alignment("PUT-157a-Solanum_tuberosum-68554"); 
##

use DBI;
use XML::Dumper;

use Bio::AlignIO;
use Bio::SimpleAlign;

use Carp;

use Set::Scalar;

my $dumper = new XML::Dumper;

my $dbh;

## constructor requires path to database file
sub new {
   my ($class, %args) = @_;

    my $db = '';
    my $user = '';
    my $password = '';

    if (defined($args{'db'})) {
        $db = $args{'db'};
        $user = $args{'user'};
        $password = $args{'password'}
    } else {
        croak "Constructor requires parameter 'db'";
    }
    if (defined($args{'user'})) {
        $user = $args{'user'};
    } else {
        croak "Constructor requires parameter 'user'";
    }
    if (defined($args{'password'})) {
        $password = $args{'password'}
    } else {
        croak "Constructor requires parameter 'password'";
    }

    _init_db($db, $user, $password);

    return bless {}, ref($class) || $class;
}

## initialize the alignment db
sub _init_db {
    my ($db_file) = @_;

    $dbh = DBI->connect("dbi:Pg:dbname=solcomp_search;host=pg", 'access', 'access');
}

## returns an alignment object from the database
sub get_alignment {
    my ($self, $id) = @_;
    
    my $xml = $self->get_alignment_xml($id);

    if (defined($xml)) {
        return $dumper->xml2pl($xml);
    } else {
        carp "get_alignment returned undef";
        return undef;
    }
}

## returns xml in case there is ever a need for that
sub get_alignment_xml {
    my ($self, $id) = @_;

    my $msa_id = $self->has_alignment($id);
    if ($msa_id) {
        my $sth = $dbh->prepare("select m.value as xml from msa m where m.id = ?");
        $sth->execute($msa_id);
        my ($xml) = $sth->fetchrow_array();
        return $xml;
    } else {
        carp "get_alignment_xml returned undef";
        return undef;
    }
}

## check if alignment is in DB
sub has_alignment {
    my ($self, $id) = @_;
   
    my $sth = $dbh->prepare("select m.id as id from name n, msa m where n.id = m.name_id and n.value = ?");

    $sth->execute($id);
    my ($msa_id) = $sth->fetchrow_array();

    if ($msa_id) {
        return $msa_id;
    } else {
        return 0;
    }
}

1;
