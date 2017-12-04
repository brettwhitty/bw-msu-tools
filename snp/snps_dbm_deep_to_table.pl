#!/opt/rocks/bin/perl 

## creates a summary table of SNP information from a DBM::Deep database

use strict;
use warnings;

use lib "/home/whitty/SVN/lib";
use NCBITaxonomyTree;
my $tree = new NCBITaxonomyTree(path => '/home/whitty/db/NCBI');

use Cwd qw{ abs_path };
use DBM::Deep;
use POSIX qw{ floor }; 
use File::Basename qw{ basename dirname };

my $snps_db_name         = shift @ARGV || die "Provide path to a SNPs DBM::Deep database";

$snps_db_name = abs_path($snps_db_name);
my $outdir = dirname($snps_db_name);
my $outbase = basename($snps_db_name, '.db');
my $outfile = "$outdir/$outbase.txt";

my $snps_db = new DBM::Deep($snps_db_name);

my $bins = [];
my $dbins = [];
my $total_snps;
my $count;
my $total_size;
my $put_with_snp = {};
my $snp_positions = {};
my $put_len = {};
my $put_depth = {};
my $alt_depth = {};

open(my $outfh, ">$outfile") || die "Failed to open '$outfile' for writing: $!";

print $outfh "# taxon_id\tscientific_name\tput_id\tmember_count\tseq_len\tsnp_loc\tsnp_cov\tref_base\tsnp_base\tbase_freq_a\tbase_freq_c\tbase_freq_g\tbase_freq_t\n";

foreach my $key(keys(%{$snps_db})) {

    my $species = $snps_db->{$key}->{'species'};
    my $taxon_id = $tree->get_taxon_id($species);  
    my $put_id = $snps_db->{$key}->{'put_id'};
    my $members = $snps_db->{$key}->{'member_count'};
    my $len = $snps_db->{$key}->{'size'};
 
    foreach my $snp(@{$snps_db->{$key}->{'snps'}}) {
        my $ref_base = $snp->{'ref'};
        my $snp_base = $snp->{'snp'};
        my $snp_loc = $snp->{'loc'};

        my $depth = 0;
        my @cov = ();
        foreach my $base( qw{a c g t} ) {
            $depth += $snp->{$base};
            push(@cov, uc($base).':'.$snp->{$base});
        }

        print $outfh join("\t", (
                $taxon_id,
                $species,
                $put_id,
                $members,
                $len,
                $snp_loc,
                $depth,
                $ref_base,
                $snp_base,
#                reverse sort{$a cmp $b} @cov,
                @cov,                
            ))."\n";
    }
}
