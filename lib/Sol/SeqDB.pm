#!/opt/rocks/bin/perl

package Sol::SeqDB;

## This class, along with Sol::Datastore and a MySQL database
## will be used for tracking and storing sequence files
## for the Solanaceae Genomics Resource site
##
## It has various utility methods for interacting with the
## database to store or fetch metadata related to the
## various sequences

use strict;
use warnings;
use Carp;

use DBI;
use Roman;
use Bio::SeqIO;
use Digest::MD5 qw{ md5_hex };
use Date::Simple;
#use Set::Scalar;

use Sol::DBIDatastore;

## for translating the date format
my $months = {
                'JAN' => '01',
                'FEB' => '02',
                'MAR' => '03',
                'APR' => '04',
                'MAY' => '05',
                'JUN' => '06',
                'JUL' => '07',
                'AUG' => '08',
                'SEP' => '09',
                'OCT' => '10',
                'NOV' => 11,
                'DEC' => 12
             };

    my $dbh;
{
    my $db_name     = 'sol_seq';
    my $db_server   = 'mysql';
    my $db_user     = $ENV{'DB_USER'};
    my $db_pass     = $ENV{'DB_PASS'};
   
    my $files_root = '';

    ## constructor
    sub new {
        my ($class) = @_;

        $dbh = DBI->connect("dbi:mysql:dbname=$db_name;host=$db_server", $db_user, $db_pass);
       
        my $ds = new Sol::DBIDatastore(
            #          'root'  =>  '/projects/compsol/.datastore',
            #                        'depth' =>  3,
            'db_host'   =>  'pg-dev.plantbiology.msu.edu',
            'db_name'   =>  'sol_datastore',
            'db_table'  =>  'file',
            'db_user'   =>  $db_user,
            'db_pass'   =>  $db_pass,
                                   );

        my $self = {'file' => $ds};
    
        bless $self, $class;

        return $self;
    }


## fetch GIs of a certain type
sub get_gi_hash {
    my ($self, $type, $obsolete) = @_;
    
    if (! _type_check($type)) {
        carp "Warning: get_gi_hash called with invalid type '$type'";
        return undef;
    }

    my %gi_list;

    my $query;
    if (defined($obsolete)) {
        $query = 'select gi from genbank g, type t where g.type_id = t.type_id and t.name = ?';
    } else {
        $query = 'select gi from genbank g, type t where g.type_id = t.type_id and t.name = ? and g.is_obsolete != 1';
    }
    my $sth = $dbh->prepare($query);

    $sth->execute($type);

    while (my $row = $sth->fetchrow_arrayref()) {
        $gi_list{$row->[0]} = 1;
    }

    return \%gi_list;
}

## returns a list of GIs
sub get_gi_list {
    my ($self, $type, $obsolete) = @_;

    return keys(%{$self->get_gi_hash($type, $obsolete)});
}

## fetch GIs of a certain type on a certain day
##
## WARNING: this sub ignores the obsolete flag
sub get_gi_hash_on {
    my ($self, $type, $date) = @_;
    
    my $o_date = new Date::Simple($date) 
        or confess "Fatal: Invalid date '$date', should be provided in format YYYY-MM-DD";

    if (! _type_check($type)) {
        carp "Warning: get_gi_hash called with invalid type '$type'";
        return undef;
    }

    my %gi_list;

    my $sth = $dbh->prepare('select gi from genbank g, type t where g.type_id = t.type_id and t.name = ? and date <= ?');

    $sth->execute($type, "$o_date");

    while (my $row = $sth->fetchrow_arrayref()) {
        $gi_list{$row->[0]} = 1;
    }

    return \%gi_list;
}

## returns a list of GIs available on a given date
##
## WARNING: this sub ignores the obsolete flag
sub get_gi_list_on {
    my ($self, $type, $date) = @_;

    return keys(%{$self->get_gi_hash_on($type, $date)});
}

## fetch GIs of a certain type after a certain date
##
## WARNING: this sub ignores the obsolete flag
sub get_gi_hash_after {
    my ($self, $type, $date) = @_;
    
    my $o_date = new Date::Simple($date) 
        or confess "Fatal: Invalid date '$date', should be provided in format YYYY-MM-DD";

    if (! _type_check($type)) {
        carp "Warning: get_gi_hash called with invalid type '$type'";
        return undef;
    }

    my %gi_list;

    my $sth = $dbh->prepare('select gi from genbank g, type t where g.type_id = t.type_id and t.name = ? and date > ?');

    $sth->execute($type, "$o_date");

    while (my $row = $sth->fetchrow_arrayref()) {
        $gi_list{$row->[0]} = 1;
    }

    return \%gi_list;
}

## returns a list of GIs available after a given date
##
## WARNING: this sub ignores the obsolete flag
sub get_gi_list_after {
    my ($self, $type, $date) = @_;

    return keys(%{$self->get_gi_hash_after($type, $date)});
}

## fetch list of GIs that have NULL values which would be populated
## if we had genbank files in the datastore for these records
sub get_gb_null_list {
    my ($self, $type) = @_;
    
    if (! _type_check($type)) {
        carp "Warning: get_null_list called with invalid type '$type'";
        return undef;
    }

    my @gi_list;

    my $sth = $dbh->prepare('select g.gi from genbank g, type t where g.type_id = t.type_id and t.name = ? and g.date is NULL and g.taxon_id is NULL');

    $sth->execute($type);

    while (my $row = $sth->fetchrow_arrayref()) {
        push(@gi_list, $row->[0]);
    }

    return \@gi_list;
}

## populates the checksum field of the genbank table
sub populate_gb_checksum {
    my ($self, $type) = @_;

    my @gi_list;

    my $sth = $dbh->prepare('select gi from genbank where checksum is NULL and date is NOT NULL and taxon_id is NOT NULL');

    $sth->execute();

    while (my $row = $sth->fetchrow_arrayref()) {
        push(@gi_list, $row->[0]);
    }

    foreach my $gi(@gi_list) {
        $self->update_gi($gi, {'checksum' => $self->get_gb_md5($gi)});
    }
}

## populates fields of the genbank table where value is currently null and values can be obtained from the flat file
sub populate_gb_atts {
    my ($self, $type) = @_;

    ## provide a list of attributes here to check for null values
    my @att_list = ('date', 'taxon_id', 'checksum', 'seqlen', 'accession', 'version', 'replicon');
  
    my $query = 'select gi from genbank where '.join(' is null or ', @att_list).' is null';

    my $sth = $dbh->prepare($query);

    $sth->execute();

    my @gi_list;
    while (my $row = $sth->fetchrow_arrayref()) {
        push(@gi_list, $row->[0]);
    }

    foreach my $gi(@gi_list) {
        my $gbk_atts = $self->get_gb_atts($gi);
        
        $self->update_gi($gi, $gbk_atts);
    }
}

## adds a GI to the database
sub add_gi {
    my ($self, $gi, $atts) = @_;
   
    my $type_id; 
    if (!($type_id = _type_check($atts->{'type'}))) {
        confess "Fatal: add_gi called with non-existant type '$atts->{type}'";
    }
    ## if date is provided do a date format check
    ## if taxon_id is provided do a valid taxon_id check
    ## etc.

    if (defined($self->get_gi_atts($gi))) {
        carp "Warning: Database record for GI '$gi' already exists";
        return 1;
    }

    my @cols = ('gi', 'type_id');
    my @vals = ($gi, $type_id);

    if (defined($atts->{'date'})) {
        push(@cols, 'date');
        push(@vals, $atts->{'date'});
    }
    if (defined($atts->{'taxon_id'})) {
        push(@cols, 'taxon_id');
        push(@vals, $atts->{'taxon_id'});
    }
    if (defined($atts->{'checksum'})) {
        push(@cols, 'checksum');
        push(@vals, $atts->{'checksum'});
    }
    if (defined($atts->{'accession'})) {
        push(@cols, 'accession');
        push(@vals, $atts->{'accession'});
    }
    if (defined($atts->{'version'})) {
        push(@cols, 'version');
        push(@vals, $atts->{'version'});
    }
    if (defined($atts->{'replicon'})) {
        push(@cols, 'replicon');
        push(@vals, $atts->{'replicon'});
    }
    if (defined($atts->{'seqlen'})) {
        push(@cols, 'seqlen');
        push(@vals, $atts->{'seqlen'});
    }

    my $query = "insert into genbank(".join(", ", @cols).") values ('".join("', '", @vals)."')";
        
    my $sth = $dbh->prepare($query);

    $sth->execute();

    ## do check here that query was successful
}

## add a list of GIs to database
## returns the number newly added GIs
sub add_gi_list {
    my ($self, $gi_arr_ref, $type) = @_;
    
    if (! _type_check($type)) {
        croak "Warning: add_gi_list called with invalid type '$type'";
    }
    
    my $db_gi = $self->get_gi_hash($type);

    my $insert_count = 0;
    foreach my $gi(@{$gi_arr_ref}) {
        unless ($db_gi->{$gi}) {
            $self->add_gi($gi, {'type' => $type});
            $insert_count++;
        }
    }

    return $insert_count;
}

## remove a gi number from the database
sub remove_gi {
    my ($self, $gi) = @_;

    if (! defined($gi)) {
        confess "Fatal: remove_gi called with no gi number specified";
    }

    my $sth = $dbh->prepare("delete from genbank where gi = ?");

    $sth->execute($gi);

    ## check here that query was successful
}


## update attributes of a gi record
sub update_gi {
    my ($self, $gi, $atts) = @_;
   
    ## if date is provided do a date format check
    ## if taxon_id is provided do a valid taxon_id check
    
    my @kvp = ();

    if (defined($atts->{'date'})) {
        push(@kvp, "date='$atts->{date}'");
    }
    if (defined($atts->{'taxon_id'})) {
        push(@kvp, "taxon_id=$atts->{taxon_id}");
    }
    if (defined($atts->{'checksum'})) {
        push(@kvp, "checksum='$atts->{checksum}'");
    }
    if (defined($atts->{'accession'})) {
        push(@kvp, "accession='$atts->{accession}'");
    }
    if (defined($atts->{'version'})) {
        push(@kvp, "version='$atts->{version}'");
    }
    if (defined($atts->{'replicon'})) {
        push(@kvp, "replicon='$atts->{replicon}'");
    }
    if (defined($atts->{'seqlen'})) {
        push(@kvp, "seqlen=$atts->{seqlen}");
    }

    if (scalar(@kvp) > 0) {

        my $query = "update genbank set ".join(", ", @kvp)." where gi=?";

        my $sth = $dbh->prepare($query);

        $sth->execute($gi);

        ## do check here that query was successful
    }
}

## return the type of the sequence by GI
sub get_gi_type {
    my ($self, $gi) = @_;
    
    my $sth = $dbh->prepare('select t.name from genbank g, type t where g.type_id = t.type_id and g.gi = ?');

    $sth->execute($gi);

    if (my $row = $sth->fetchrow_arrayref()) {
        return $row->[0];
    } else {
        confess "Warning: GI number '$gi' not found in the database";
        return undef;
    }
}

## returns the max date from the database
sub get_gb_max_date {
    my ($self, $type) = @_;
    
    if (! _type_check($type)) {
        croak "Warning: get_max_date called with invalid type '$type'";
    }
    
    my $sth = $dbh->prepare('select max(g.date) from genbank g, type t where g.type_id = t.type_id and t.name = ?');

    $sth->execute($type);

    if (my $row = $sth->fetchrow_arrayref()) {
        return new Date::Simple($row->[0]);
    } else {
        confess "Warning: Unable to fetch max date from database for type '$type'";
        return undef;
    }
}

## returns a Date::Simple object for today
## for use comparing against results of get_max_date
sub get_today_date {
    my ($self) = @_;

    return new Date::Simple;
}

## returns a string containing the genbank record for a specific GI number
sub get_genbank {
    my ($self, $gi) = @_;

    my $filename = "$gi.gbk";

    my $gbk_fh = $self->{'file'}->read_fh($filename);

    return do{ local $/; <$gbk_fh>; }
}

## finds genbank flat file for a given GI and returns FASTA format sequence with a GenBank format header
sub get_gb_fasta {
    my ($self, $gi) = @_;

    my $filename = "$gi.gbk";

    my $gbk_fh = $self->{'file'}->read_fh($filename);

    my $gbk_db = new Bio::SeqIO(-fh => $gbk_fh, -format => 'genbank');

    my $fasta = '';
    while (my $seq = $gbk_db->next_seq()) {
        my $head = '>gi|'.$seq->primary_id().'|gb|'.$seq->id().'.'.$seq->version().'| '.$seq->desc()."\n";
        my $seq  = $seq->seq();
        $seq =~ s/(\S{1,60})/$1\n/g;
        $fasta .= $head.$seq;
    }

    return $fasta;
}

## get chromosome assignment for a GI number
sub get_gb_atts {
    my ($self, $gi) = @_;

    my $atts = {};

    my $filename = "$gi.gbk";

    my $gbk_fh = $self->{'file'}->read_fh($filename);

    unless (defined($gbk_fh)) {
        return undef;
    }

    my $gbk_db = new Bio::SeqIO(-fh => $gbk_fh, -format => 'genbank');

    my $clone = undef;
    my $chrom = undef;
    
    my $seq = $gbk_db->next_seq() or confess "Fatal: Failed to read sequence object from '$filename'";

    my @dates = $seq->get_dates();
    for (my $i = 0; $i < scalar(@dates); $i++) {
        my ($d, $m, $y) = split("-", $dates[$i]);
        $m = $months->{$m};
        my $o_date = new Date::Simple("$y-$m-$d") or confess "Fatal: Oddly formatted date encountered '$dates[$i]' in GI '$gi'";
        $dates[$i] = $o_date;
    }
    ## sort most recent first    
    @dates = sort{$b cmp $a} @dates;
    my $date = $dates[0];

    my $taxon_id = $seq->species->id();
    my ($accession, $version) = ($seq->id(), $seq->version());
    my $seq_len = length($seq->seq());
    my $checksum = uc(md5_hex($seq->seq()));
    my $division = $seq->division();

    my @feats = $seq->get_SeqFeatures() or confess "Fatal: Failed to retrieve feature for sequence '$gi'";
        
    ## get clone attribute
    if ($feats[0]->has_tag('clone')) {
        $clone = ($feats[0]->get_tag_values('clone'))[0];
    } elsif ($seq->desc =~ /clone (\w+)/) {
        $clone = $1;
    }

    ## feature has chromosome tag
    if ($feats[0]->has_tag('chromosome')) {
        $chrom = lc(($feats[0]->get_tag_values('chromosome'))[0]);
    } elsif ($feats[0]->has_tag('organelle')) {
        $chrom = lc(($feats[0]->get_tag_values('organelle'))[0]);
    } elsif ($seq->desc =~ /chromosome (\w+)/) {
        $chrom = lc($1);
    } else {
            ## do something here for troubleshooting?
#            carp "Warning: Feature for GI '$gi' has no chromosome or organelle tags";
    }

    if (defined($chrom)) {
        if (isroman($chrom)) {
            $chrom = arabic($chrom);
        } elsif ($chrom eq '0') {
            $chrom = 'unknown';
        } elsif ($chrom =~ /^0[123456789]+/) {
            $chrom =~ s/^[0]+//;
        } elsif ($chrom =~ /plast/i) {
            $chrom = 'chloroplast';
        } elsif ($chrom =~ /mito/i) {
            $chrom = 'mitochondrion';
        }
        if ($chrom =~ /^\d+$/) {
            $chrom = "chr$chrom";
        }
    } else {
        $chrom = 'unknown';
    }

    ## this may turn out to be incorrect later, but unlikely
    my $type = ($seq_len > 10000) ? 'contig' : 'read';    

    $atts = {
                'type'      =>  $type,      ## for use by add_gi sub
                'date'      =>  $date, 
                'accession' =>  $accession,
                'version'   =>  $version,
                'taxon_id'  =>  $taxon_id,
                'checksum'  =>  $checksum,
                'seqlen'    =>  $seq_len,
                'replicon'  =>  $chrom,
                'clone'     =>  $clone,     ## not currently used in DB
                'division'  =>  $division,  ## not currently used in DB
            };

    return $atts;
}

## fetches attributes stored in the database for a given GI
sub get_gi_atts {
    my ($self, $gi, $quiet) = @_;

    my $sth = $dbh->prepare('select * from genbank where gi = ?');

    $sth->execute($gi);

    if (my $row = $sth->fetchrow_hashref()) {
        return $row;
    } else {
        unless ($quiet) {
            carp "Warning: GI number '$gi' not found in the database";
        }
        return undef;
    }
}

## fetch gi number for accession
sub get_gb_gi {
    my ($self, $accession, $quiet) = @_;

    my $sth = $dbh->prepare('select gi from genbank where accession = ? and is_obsolete = 0 order by version, date limit 1');

    $sth->execute($accession);

    if (my $row = $sth->fetchrow_arrayref()) {
        return $row->[0];
    } else {
        unless ($quiet) {
            carp "Warning: Accession '$accession' not found in the database";
        }
        return undef;
    }
}

## fetch accession number for gi
sub get_gb_accession {
    my ($self, $gi, $quiet) = @_;

    my $sth = $dbh->prepare('select accession from genbank where gi = ?');

    $sth->execute($gi);

    if (my $row = $sth->fetchrow_arrayref()) {
        return $row->[0];
    } else {
        unless ($quiet) {
            carp "Warning: GI '$gi' not found in the database";
        }
        return undef;
    }
}


## returns an MD5 hash for the specified sequence
sub get_gb_md5 {
    my ($self, $gi) = @_;

    my $filename = "$gi.gbk";

    my $gbk_fh = $self->{'file'}->read_fh($filename);

    my $gbk_db = new Bio::SeqIO(-fh => $gbk_fh, -format => 'genbank');

    my $seq = $gbk_db->next_seq() or confess "Fatal: Failed to retrieve sequence for '$gi'";

    return uc(md5_hex($seq->seq()));
}

## public method for getting types
## (may not be needed)
sub get_types {
    my ($self) = @_;
    return values %{_get_types()};
}

## update the contents of the taxon table based on the taxon_id's present in genbank
sub add_taxon {
    my ($self) = @_;
    
    my $taxon_hash = _get_taxon_hash();

    my @gi_taxon_ids = _get_gi_taxon_list();
    
    use NCBITaxonomyTree;
    my $tree = new NCBITaxonomyTree();

    my $sth = $dbh->prepare("insert into taxon (taxon_id, scientific_name, common_name) values (?, ?, ?)");
    foreach my $taxon_id(@gi_taxon_ids) {
        if (! defined($taxon_hash->{$taxon_id})) {
            my $scientific_name = $tree->get_scientific_name($taxon_id);
            $sth->execute($taxon_id, $scientific_name, '');
            ## check for success of query
        }
    }
}

## sets records in genbank to obsolete where taxon is non-solanaceae
sub set_gb_obsolete_not_sol {
    my ($self) = @_;
    
    ## get array of taxon_id from taxon table 
    my @taxon_ids = _get_taxon_list();

    use NCBITaxonomyTree;
    my $tree = new NCBITaxonomyTree();

    my @non_sol_taxon_id = ();
    foreach my $taxon_id(@taxon_ids) {
        unless ($tree->id_has_ancestor_name($taxon_id, 'Solanaceae')) {
            push(@non_sol_taxon_id, $taxon_id);
        }
    }
    if (scalar(@non_sol_taxon_id) > 0) {
        my $sth = $dbh->prepare("update genbank set is_obsolete = 1 where taxon_id in (".join(", ", @non_sol_taxon_id).")");
        $sth->execute();
    }
}

## sets records in genbank to obsolete where accession is not the max version 
sub set_gb_obsolete_not_current {
    my ($self) = @_;
    
    my $version_query = 'select accession, max(version), count(version) from genbank group by accession having count(version) > 1';

    my $sth = $dbh->prepare($version_query);
    $sth->execute();

    my $acc_ver = {};
    while (my $row = $sth->fetchrow_arrayref) {
        $acc_ver->{$row->[0]} = $row->[1];
    }

    $sth = $dbh->prepare("update genbank set is_obsolete = 1 where accession = ? and version != ?");
    foreach my $acc(keys(%{$acc_ver})) {
        $sth->execute($acc, $acc_ver->{$acc});
    }
}

sub get_gb_release_date {
    my ($self, $release_number) = @_;

    my %gi_list;

    my $sth = $dbh->prepare('select date from releases where release_id = ?');

    $sth->execute($release_number);

    if (my $row = $sth->fetchrow_arrayref()) {
        my $o_date = new Date::Simple($row->[0]) or confess "Fatal: Unable to create date object from string '$row->[0]'";
        return $o_date;
    } else {
        confess "Fatal: get_release_date was given an invalid release number '$release_number'";
    }
}


sub get_taxonomy_hash {
    my $taxonomy = {};
    
    my $sth = $dbh->prepare('select * from taxon');

    $sth->execute();

    while (my $row = $sth->fetchrow_hashref()) {
        my $taxon_id = $row->{'taxon_id'};
        delete $row->{'taxon_id'};
        $taxonomy->{$taxon_id} = $row;
    }

    return $taxonomy;
}

}

## returns a hash of taxon_ids from the taxon table
sub _get_taxon_hash {
    my $taxon_ids = {};
    
    my $sth = $dbh->prepare('select taxon_id from taxon');

    $sth->execute();

    while (my $row = $sth->fetchrow_arrayref()) {
        $taxon_ids->{$row->[0]} = 1;
    }

    return $taxon_ids;
}

sub _get_taxon_list {
    return keys(%{_get_taxon_hash()});
}

sub _get_gi_taxon_hash {
    my $taxon_ids = {};
    
    my $sth = $dbh->prepare('select distinct(taxon_id) from genbank');

    $sth->execute();
    
    while (my $row = $sth->fetchrow_arrayref()) {
        $taxon_ids->{$row->[0]} = 1;
    }

    return $taxon_ids;
}

sub _get_gi_taxon_list {
    return keys(%{_get_gi_taxon_hash()});
}

## returns a hash of types from the type table
sub _get_types {
    my ($self) = @_;

    my $types = {};
    
    my $sth = $dbh->prepare('select * from type');

    $sth->execute();

    while (my $row = $sth->fetchrow_arrayref()) {
        $types->{$row->[0]} = $row->[1];
    }

    return $types;
}

## checks whether a type is valid
sub _type_check {
    my ($type) = @_;

    unless (defined($type)) {
        confess "Fatal: _type_check called with null type";
    }

    my %types = reverse %{_get_types()};

    my $val = (defined($types{$type})) ? $types{$type} : 0;

    return $val;
}

1;
