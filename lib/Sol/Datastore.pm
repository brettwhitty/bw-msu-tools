#!/usr/bin/perl

package Sol::Datastore;

use strict;
use warnings;
use Carp;

use Digest::MD5 qw{ md5_hex };
use File::Path qw{ mkpath };
use File::Copy qw{ };
#use File::Basename qw{ dirname };

sub new {
    my ($class, %args) = @_;

    unless (defined($args{'root'})) {
        croak "Constructor requires 'root' argument";
    }

    my $depth  = (defined($args{'depth'})) ? $args{'depth'}    : 3;

    my $root = $args{'root'};
    unless ($root =~ /\/$/) {
        $root .= "/";
    }

    my $self = {
                '_root'     => $root,
                '_depth'    => $depth,
               };

    return bless $self, $class;
}

## copy or move a file to the datastore
sub copy {
    my ($self, $src_file_path, $target_name, $move) = @_;

    unless (-e $src_file_path) {
        carp "Attempt to add non-existant file '$src_file_path' to datastore failed";
        return undef;
    }
   
    my $dest_dir = $self->id_to_dir($target_name);

    if ($move) {
        File::Copy::move($src_file_path, $dest_dir);
    } else {
        File::Copy::copy($src_file_path, $dest_dir);
    }

    return $dest_dir;
}

## use copy to do a move
sub move {
    my ($self, $src_file_path, $target_name) = @_;

    $self->copy($src_file_path, $target_name, 1);
}

sub remove {
    my ($self, $file_name) = @_;

    my $location = $self->locate($file_name);

    unlink($location);
    ## test that this succeeded    
}

sub locate {
    my ($self, $file_name) = @_;

    my $location = $self->id_to_dir($file_name).$file_name;

    return $location;
}

sub create_fh {
    my ($self, $file_name) = @_;

    my $location = $self->locate($file_name);

    if (-e "$location") {
        carp "Warning: File '$location' already exists and will be overwritten";
    }

    open(my $out_fh, ">$location") || croak "Failed to open file '$location' for writing: $!";

    return $out_fh;
}

sub read_fh {
    my ($self, $file_name) = @_;

    my $location = $self->locate($file_name);

    unless (-e "$location") {
        carp "Warning: File '$location' does not exist";
        return undef;
    }

    open(my $in_fh, "$location") || croak "Failed to open file '$location' for reading: $!";

    return $in_fh;
}

## Borrowed from Datastore::MD5
sub id_to_dir {
    my ($self, $id) = @_;

    my $digest = uc md5_hex($id); # the hex string, uppercased

  # make sure we have enough hex digits to build the dir string.
    if ($self->{'_depth'} > 16) {
        return(undef);
    }

    my $dir = $self->{'_root'};
    for (my $i = 0; $i < $self->{'_depth'}; $i++) {
        $dir .= substr($digest, $i * 2, 2) . "/";
    }
  # $dir .= $id . "/";

    unless (-d $dir) {
        mkpath($dir);
    }

    return $dir;
}

1;
