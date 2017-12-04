#!/opt/rocks/bin/perl

package Sol::DBIDatastore;

use strict;
use warnings;
use Carp;

use DBI;
use DBIx::BLOB::Handle qw( :INTO_STATEMENT );

use Digest::MD5 qw{ md5_hex };
#use File::Path qw{ mkpath };
#use File::Copy qw{ };
#use File::Basename qw{ dirname };

use File::Temp qw{ tempfile };
$File::Temp::KEEPALL = 0;

my $db_table = 'sol_datastore';

sub new {
    my ($class, %args) = @_;

    unless (defined($args{'db_name'})) {
        croak "Constructor requires 'db_name' argument";
    }
    unless (defined($args{'db_host'})) {
        croak "Constructor requires 'db_host' argument";
    }
    unless (defined($args{'db_user'})) {
        croak "Constructor requires 'db_user' argument";
    }
    unless (defined($args{'db_pass'})) {
        croak "Constructor requires 'db_pass' argument";
    }
    if (defined($args{'db_table'})) {
        $db_table = $args{'db_table'};
    }

    my $dbh = DBI->connect("dbi:Pg:dbname=".$args{'db_name'}.";host=".$args{'db_host'}, $args{'db_user'}, $args{'db_pass'}, {RaiseError => 1, PrintError => 0 })
        || confess "Failed to connect to database!";
    $dbh->{LongTruncOk} = 1;

    my $self = {
                '_dbh'     => $dbh,
               };

    return bless $self, $class;
}

## copy or move a file to the datastore
sub copy {
    my ($self, $src_file_path, $target_name, $move, $force) = @_;

#    unless (-e $src_file_path) {
#        carp "Attempt to add non-existant file '$src_file_path' to datastore failed";
#        return undef;
#    }
#
#    my $dest_dir = $self->id_to_dir($target_name);
#
#    if ($move) {
#        File::Copy::move($src_file_path, $dest_dir);
#    } else {
#        File::Copy::copy($src_file_path, $dest_dir);
#    }

#    confess "copy not implemented by DBIDatastore";
    ## check whether file with the same name exists in the database first
    my $name_query = "select count(id) from $db_table where name = ?";
    my $sth= $self->{'_dbh'}->prepare($name_query);
    $sth->execute($target_name);
    my $row_count = ($sth->fetchrow_array())[0] || 0;
    
    ## if the file exists require that we force the overwrite 
    if ($row_count) {
        if ($force) {
            carp "Filename '$target_name' already exists in the database, and will be overwritten.";
        } else {
            carp "Filename '$target_name' already exists in the database, and will not be overwritten without force.";
            return 0;
        } 
    }
    my $oid = $self->{'_dbh'}->func($src_file_path, 'lo_import');
    #my $query = "insert into $db_table (name, data) values (?, lo_import(?))";
    #$sth = $self->{'_dbh'}->prepare($query) or carp $DBI::errstr;
    #$sth->execute($target_name, $src_file_path) or carp $DBI::errstr;
    if (! $oid) {
        carp "Object import failed";
        return 0;
    }
    my $query = "insert into $db_table (name, data) values (?, ?)";
    $sth = $self->{'_dbh'}->prepare($query) or carp $DBI::errstr;
    $sth->execute($target_name, $oid) or carp $DBI::errstr;

    if ($move) {
        unlink($src_file_path);
    }

    return 1;
}

## use copy to do a move
sub move {
    my ($self, $src_file_path, $target_name, $force) = @_;

#    confess "move not implemented by DBIDatastore";

    $self->copy($src_file_path, $target_name, 1, $force);
}

sub remove {
    my ($self, $file_name) = @_;

    #my $location = $self->locate($file_name);

    confess "remove not implemented by DBIDatastore";

    ##unlink($location);
    ## test that this succeeded    
}

sub locate {
    my ($self, $file_name) = @_;

#    my $location = $self->id_to_dir($file_name).$file_name;

#    return $location;
    confess "locate not supported by DBIDatastore";
}

sub create_fh {
    my ($self, $file_name) = @_;

    confess "create_fh is unimplemented";
#
#    my $location = $self->locate($file_name);
#
#    if (-e "$location") {
#        carp "Warning: File '$location' already exists and will be overwritten";
#    }
#
#    open(my $out_fh, ">$location") || croak "Failed to open file '$location' for writing: $!";
#
#    return $out_fh;
}

sub read_fh {
    my ($self, $file_name) = @_;

    my (undef, $temp_name) = tempfile(DIR => '/dev/shm');
    #my (undef, $temp_name) = tempfile(DIR => '/dev/shm', UNLINK => 1, CLEANUP => 1, OPEN => 0);

    my $query = "select data from $db_table where name = ?";
    my $sth= $self->{'_dbh'}->prepare($query) or carp $DBI::errstr;
    $sth->execute($file_name) or carp $DBI::errstr;
    my $oid = ($sth->fetchrow_array())[0] || 0;

    if ($oid) {

        $self->{'_dbh'}->func($oid, $temp_name, 'lo_export') or carp $DBI::errstr;

        open my $fh, '<', $temp_name or confess "Failed to open temp file '$temp_name' for reading: $!";

        ## trying a trick here to see if the file will get unlinked eventually
        ## because File::Temp-based cleanup isn't working
        unlink($temp_name);

        return $fh;
    } else {
        return undef;
    }

#    my $location = $self->locate($file_name);

#    unless (-e "$location") {
#        carp "Warning: File '$location' does not exist";
#        return undef;
#    }

#   open(my $in_fh, "$location") || croak "Failed to open file '$location' for reading: $!";

#    my $query = "select lo.data from $db_table t, pg_largeobject lo where t.data = lo.loid and t.name = ?";
#    my $sth= $self->{'_dbh'}->prepare($query);
#    $sth->execute($file_name);
    #$sth->fetch;
    #my $in_fh = new DBIx::BLOB::Handle->new($sth, 0, 4096);
#    my $in_fh = $sth->blob_as_handle(0, 4096);

    
}

## Borrowed from Datastore::MD5
sub id_to_dir {
    my ($self, $id) = @_;

    confess "Method id_to_dir is unsupported by DBIDatastore";

  #  my $digest = uc md5_hex($id); # the hex string, uppercased


  # make sure we have enough hex digits to build the dir string.
  #  if ($self->{'_depth'} > 16) {
  #      return(undef);
  #  }

  #  my $dir = $self->{'_root'};
  #  for (my $i = 0; $i < $self->{'_depth'}; $i++) {
  #      $dir .= substr($digest, $i * 2, 2) . "/";
  #  }
  # $dir .= $id . "/";

  #  unless (-d $dir) {
  #      mkpath($dir);
  #  }

  #  return $dir;
}

1;
