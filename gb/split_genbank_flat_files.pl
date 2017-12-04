#!/opt/rocks/bin/perl

## simple script to split a file containing multiple genbank flat file records into individual files
## 
## will probably break easily

use strict;
use warnings;
use Carp;

use Getopt::Long;
use File::Basename qw( basename dirname);
use File::Path;
use Time::Piece;

use Bio::Taxon;
use Bio::SeqIO;

use FileCache qw( cacheout cacheout_close );

use lib '/home/whitty/SVN/lib';
use NCBITaxonomyTree;

my ($input, $output, $multi);

GetOptions(
    'input|i=s'     =>  \$input,
    'output|o=s'    =>  \$output,
    'multi|m!'      =>  \$multi,
);

my $lt = localtime();
my $current_date = $lt->ymd('');

## a hash to store FileCache filehandles
my $file_handles = {};

if (! defined($output)) {
    die "must provide a working directory with the --output flag";
}
my $work_dir = $output;

my $tree = new NCBITaxonomyTree();
#my $dbh = new Bio::DB::Taxonomy(                                
#                                -source        => 'flatfile',
#                                -directory     =>  '/home/whitty/projects/db/NCBI/Taxonomy_dump/index',
#                                -nodesfile     =>  '/home/whitty/projects/db/NCBI/Taxonomy_dump/nodes.dmp',
#                                -namesfile     =>  '/home/whitty/projects/db/NCBI/Taxonomy_dump/names.dmp',
#                               );

my $outfile;
my $outfh;
my $acc;
my $org;
my $tax_id;
my $def;
my $record_date;
my $division;
my $mol_type;
my $topology;
my $seq_len;
my $chrom = 'unknown';

my %stats = ();


my $gbk_out_string = '';

my $infh;
if (defined($input)) {
    open $infh, '<', $input || die "Failed to open '$input' for reading: $!";
} else {
    $infh = \*STDIN;
}


print STDERR "Processing input file '$input'...\n";

my $count = 0;
while (<$infh>) {
    if (/^$/) { next; }
   
    $gbk_out_string .= $_;
    
    ## match the LOCUS line and parse out accession
    if (/^LOCUS\s+(\S+)\s+(\d+)\sbp\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
        ($acc, $seq_len, $mol_type, $topology, $division, $record_date) = ($1, $2, $3, $4, $5, $6);   
    } 
    
    ## match on chromosome
    if (/\/chromosome=\"([^\"]+)\"/) {
        $chrom = $1;
    }
    ## match on organelle
    if (/\/organelle=\"([^\"]+)\"/) {
            $chrom = $1;
    }

    ## match the organism line and parse out organism name
    if (/^\s+ORGANISM\s+(.*)$/) {
        $org = $1;
        #my $taxon = $dbh->get_taxon(-name => $org) || confess "No such taxon: $org";
        #$tax_id = $taxon->id;
        $tax_id = $tree->get_taxon_id($org) || confess "No such taxon: $org";
    }

    ## match the definition line
    if (/^DEFINITION\s+(.*)$/) {
       $def = $1;
    }

    ## match the end of the record
    if (/^\/\/$/) {
        unless ($acc) {
            confess "Got to end of record without parsing accession";
        }
        unless ($org) {
            confess "Got to end of record without parsing organism";
        }
    
        ## clean up chrom if we got one
        if ($chrom) {
            $chrom =~ s/plastid\://;
            $chrom =~ s/^[0]+//;
        }
        
        ## generate an organism string to use as part of the directory name
        my $org_string = $org;
        $org_string =~ s/ /_/g;

        my $outdir = "$work_dir/$tax_id.$org_string";
       
        my $out_gbk;
        if ($multi) {
            $out_gbk = "$outdir/$tax_id.$org_string.gbk";
        } else {
            $out_gbk = "$outdir/$acc.gbk";
        }
    
        ## if the output dir doesn't exist, create it
        unless (-d $outdir) {
            mkpath([$outdir]) || confess "Failed to make output dir '$outdir': $!";
            print STDERR "Created directory '$outdir'\n";
            ## write taxon id and organism name to a .taxonomy file
            open (OUT, ">$outdir/.version") || confess "Failed to open .version file for writing: $!";
            print OUT "$org\t$tax_id\t$current_date\n";
        }
      
        my $gbk_fh = get_filehandle($out_gbk);
        
        no strict;
        print $gbk_fh $gbk_out_string; 
        $gbk_out_string = '';
        use strict;

        ## store the stats for outputting later
        $stats{$acc} = [$tax_id, $org, $acc, $chrom, $topology, "$acc.gbk", $record_date, '', $def];

        ## reset all variables
        $acc = '';
        $org = '';
        $tax_id = '';
        $def = '';
        $record_date = '';
        $division = '';
        $mol_type = '';
        $topology = '';
        $seq_len = '';
        $chrom = 'unknown';
       
    }
}

print STDERR "Converting genbank flat files to fasta...\n";

## create fasta files and summary stats tables for each genbank file
foreach my $gbk(keys(%{$file_handles})) {

    unless ($gbk =~ /\.gbk$/) { next; }
   
    cacheout_close $file_handles->{$gbk};
    
    my $outdir = dirname($gbk);
    my $base   = basename($gbk, '.gbk');
    my $out_fsa = "$outdir/$base.fsa";
    
    ## convert the file to fasta
    open(my $gbk_in, "$gbk") || die "Failed to open '$gbk' for reading: $!";
       
    open(my $fsa_fh, ">$out_fsa") || die "Failed to open '$out_fsa' for writing: $!";
    
    my $in  = Bio::SeqIO->new('-fh' => $gbk_in,  '-format' => 'genbank');
    my $out = Bio::SeqIO->new('-fh' => $fsa_fh, '-format' => 'fasta');

    while ( my $seq = $in->next_seq() ) {
        $out->write_seq($seq);
    }
    
    my $list_file_name = "$outdir/.list";
    my $list_fh = get_filehandle($list_file_name);
    
    ## get fasta checksums   
    my @checksums = `fastachecksum $out_fsa`;

    foreach my $checksum(@checksums) {
        chomp $checksum;
       
        $checksum =~ s/^(\S+) (\S+) (\S+)/$1_$2/;

        my $acc = $3;
        $stats{$acc}->[7] = $checksum;

        no strict;
        print $list_fh join("\t", @{$stats{$acc}})."\n";
        use strict;
    }
}

## use FileCache to provide filehandles for a specific string, either accession # or
## taxonomy-based (for multi) or any other string
sub get_filehandle {

    my ($path) = @_;

    unless (defined($file_handles->{$path})) {
        $file_handles->{$path} = cacheout '>', $path;
    }
    
    return $file_handles->{$path};    

}
