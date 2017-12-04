#!/opt/rocks/bin/perl

## script to do fetching of PlantGDB PUTs for the solanaceae genomes in an automated fashion
##
## Brett Whitty 2008
## whitty@msu.edu

use lib "/home/whitty/SVN/lib";
use strict;
use warnings;
use Carp;
use Cwd;
use LWP::Simple;
use HTTP::Cache::Transparent;
use HTML::TreeBuilder;
use File::Fetch;
use File::Path;
use Archive::Extract;
use Bio::Taxon;
use NCBITaxonomyTree;
use DBM::Deep;

use Getopt::Long;

$File::Fetch::BLACKLIST = [qw|lwp|];

my $puts_repository;
my $append_taxon_id;
my $taxon_filter;

my $result = GetOptions(
                            'repository|r=s'        => \$puts_repository,
                            'append_taxon_id|t!'    => \$append_taxon_id,
                            'taxon_filter=s'        => \$taxon_filter,
                       );                    

## NCBI taxonomy tree object will be used to lookup taxon ids
my $tree = new NCBITaxonomyTree( path => '/home/whitty/db/NCBI/' );

## PUTs repository directory
$puts_repository = ($puts_repository) ? $puts_repository : "/home/whitty/projects/lineage_specific/puts/puts_db";

$taxon_filter = ($taxon_filter) ? $taxon_filter : 'Solanaceae';

## init the repository
initialize_repository($puts_repository);

my $puts_db = DBM::Deep->new( "$puts_repository/.puts.db" );

## enable caching of the LWP::Simple::get request
HTTP::Cache::Transparent::init( {
    BasePath => getcwd().'/.http_cache',
                                } );

## fetch the PlantGDB PUTs assembly status page                           
my $document = get('http://www.plantgdb.org/prj/ESTCluster/AssemblyRawData.php');

## parse the page to generate a hash of info for the PUTs
#my $puts_hash_ref = parse_assembly_data_html($document);
my $puts = parse_assembly_data_html($document);

### store in DBM::Deep db
#$puts_db->import($puts);

my $puts_fetched = {};
foreach my $puts_ref(values(%{$puts})) {
    if ($tree->id_has_ancestor_name($puts_ref->{'taxon_id'}, $taxon_filter)) {    
        fetch_puts($puts_repository, $puts_ref);
        $puts_fetched->{$puts_ref->{'taxon_id'}} = $puts_ref;
    }
}
### store in DBM::Deep db
#use Data::Dumper;
#print Dumper $puts_fetched;
$puts_db->import($puts_fetched);

sub parse_assembly_data_html {
    my ($document) = @_;
    
    my $parser = new HTML::TreeBuilder();
    $parser->parse($document);
    $parser->elementify();

    my @tables = $parser->find('table');
    #use Data::Dumper;
    #print Dumper @tables[1];
    my @rows = $tables[3]->find('tr');

#    use Data::Dumper;
#    print Dumper @rows;

    ## discard table headings
    shift @rows;
    shift @rows;

    my $puts = {};
    foreach my $row(@rows) {

        my @as = $row->find('a');
        unless (@as) {last;}

        my $href = $as[0]->attr('href');
        $href =~ /dir=(.*)/ || die "Failed parsing directory from '$href'";
        my $path = $1;

        $path =~ /\/([^\/]+)\/current_version/ || die "Failed matching species base name";
        my $species_base = $1;

        my $put_fasta_name = "$species_base.mRNA.PUT.fasta.bz2";
        my $put_member_name = "$species_base.PUT_member.txt.bz2";
        my $put_member_fasta_name = "$species_base.mRNA.PUTmemberSequence.fasta.bz2";
        my $put_duplicate_fasta_name = "$species_base.mRNA.Subsequence.fasta.bz2";
        my $put_dup_map_file_name = "$species_base.Sub_Sequence.txt.bz2";
        my $put_alignment_file_name = "$species_base.alignment.txt.bz2";
        my $put_fasta_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/$put_fasta_name";
        my $put_member_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/$put_member_name";
        my $put_member_fasta_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/raw_data/$put_member_fasta_name";
        my $put_alignment_file_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/$put_alignment_file_name";
        my $put_duplicate_fasta_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/raw_data/$put_duplicate_fasta_name";
        my $put_dup_map_file_url = "http://www.plantgdb.org/download/Download/Sequence/ESTcontig/$species_base/current_version/$put_dup_map_file_name";

        my @cells = $row->find('td');

        my @vals = ();
        foreach my $cell(@cells) {
            push(@vals, $cell->as_text());
        }
    
        my $taxon_id = $tree->get_taxon_id($vals[0]);
    
        $puts->{$species_base} = {
                                'name'          =>  $vals[0],
                                'current_ests'  =>  $vals[1],
                                'build_ests'    =>  $vals[2],
                                'base_name'     =>  $species_base,
                                'taxon_id'      =>  $taxon_id,
                                'member_count'  =>  $vals[8],
                                'puts_count'    =>  $vals[9],
                                'build_date'    =>  $vals[10],
                                'version'       =>  $vals[11],
                                'fasta_url'     =>  $put_fasta_url,
                                'member_url'    =>  $put_member_url,
                                'member_fasta_url'      =>  $put_member_fasta_url,
                                'duplicate_fasta_url'   =>  $put_duplicate_fasta_url,
                                'duplicate_file_url'    =>  $put_dup_map_file_url,
                                'alignment_file_url'    =>  $put_alignment_file_url,
                                 };
            
    #    print STDERR "Fetching '$put_fasta_url'...\n";    
    #    my $ff = new File::Fetch( uri => $put_fasta_url);
    #    my $loc = $ff->fetch() or die $ff->error;
    }

    return $puts;
}

sub initialize_repository {
    my ($dir) = @_;

    unless (-d $dir) {
        mkpath([$dir]);
    }
#    unless (-d $dir."/.by_taxon_id") {
#        mkpath([$dir."/.by_taxon_id"]);
#    }
}

sub fetch_puts {
    my ($repository, $puts_ref) = @_;
    
    my $taxon_id = $puts_ref->{'taxon_id'};

    my $dir;
    if ($append_taxon_id) {
        $dir = "$repository/$taxon_id.$puts_ref->{base_name}";
    } else {
        $dir = "$repository/$puts_ref->{base_name}";
    }
#    my $tdir = "$repository/.by_taxon_id/$taxon_id.$puts_ref->{base_name}";
    
    unless (-d $dir) {
        mkpath([$dir]);
    }
    
#    unless (-d $tdir) {
#        mkpath([$tdir]);
#    }
    
    print STDERR "Fetching: ".$puts_ref->{'fasta_url'}."\n";
    my $fasta_file = fetch($dir, $puts_ref->{'fasta_url'});
    
    ## symlink $dir to /.by_taxon_id/$taxon_id/
    ## symlink $dir/$fasta_file to /.by_taxon_id/$taxon_id/$taxon_id.fsa
#    my $symlink_exists = eval { symlink("$dir/$fasta_file","$tdir/$fasta_file"); 1 };
    
#    unless ($symlink_exists) {
#        confess "Failed to make symlink from '$dir/$fasta_file' to '$tdir/$fasta_file'";
#    }
    
    print STDERR "Fetching: ".$puts_ref->{'member_url'}."\n";
    my $member_file = fetch($dir, $puts_ref->{'member_url'});
    
    print STDERR "Fetching: ".$puts_ref->{'member_fasta_url'}."\n";
    my $member_fasta_file = fetch($dir, $puts_ref->{'member_fasta_url'});
    
    print STDERR "Fetching: ".$puts_ref->{'duplicate_fasta_url'}."\n";
    my $duplicate_fasta_file = fetch($dir, $puts_ref->{'duplicate_fasta_url'});
    print STDERR "Fetching: ".$puts_ref->{'duplicate_file_url'}."\n";
    my $duplicate_mapping_file = fetch($dir, $puts_ref->{'duplicate_file_url'});

    print STDERR "Fetching: ".$puts_ref->{'alignment_file_url'}."\n";
    my $alignment_file = fetch($dir, $puts_ref->{'alignment_file_url'});
    
    open (OUT, ">$dir/.version") || confess "Failed to open '$dir/.version' for writing";
    print OUT "$puts_ref->{name}\t$puts_ref->{taxon_id}\t$puts_ref->{version}\n";
    close OUT;
    
}
 
## will fetch and unarchive files
sub fetch {
    my ($dir, $url) = @_;
    
    ## fetch file
    my $ff = new File::Fetch( uri => $url );
    my $loc = $ff->fetch( to => $dir ) or carp $ff->error;
    
    unless ($loc) {
        return undef;
    }

    ## extract file
    my $ae = new Archive::Extract( archive => $loc );
    my $ok = $ae->extract( to => $dir ) or carp $ae->error;
    unlink($loc);

    return $ae->files->[0];
}
