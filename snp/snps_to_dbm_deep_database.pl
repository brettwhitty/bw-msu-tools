#!/opt/rocks/bin/perl

## processes the SNP script output and generates a DBM::Deep database

use strict;
use warnings;

use DBM::Deep;
use Cwd qw{ abs_path };
use File::Basename qw{ basename dirname };

use lib "/home/whitty/SVN/lib";
use NCBITaxonomyTree;

my $tree = new NCBITaxonomyTree(path => '/home/whitty/db/NCBI');

my $snp_file = shift @ARGV || die "Provide output file from SNPs script as input";

$snp_file = abs_path($snp_file);
my $outdir = dirname($snp_file);
my $outname = basename($snp_file).".db";
my $outdb = "$outdir/$outname";
my $outlog = $outdb.".log";

my $snps = {};

#my $flag=0;
my ($id, $size, $members, $taxon_id, $species);

open (my $logfh, ">$outlog") || die "'$outlog': $!";

print $logfh "Processing '$snp_file'\n";
open(my $infh, $snp_file) || die "'$snp_file': $!";
while (<$infh>) {
    chomp;
    
    if (/^>(\S+) size (\d+) with (\d+) members/) {
        ($id, $size, $members) = ($1, $2, $3);

        $id =~ /^PUT-\d+[a-z]+-([^-]+)/;
        $species = $1;
        my $put_id = $id;
        $id = lc($id);

#        unless($flag) {
#            print STDERR $species."\n";
#            $flag = 1;
#        }

        $species =~ s/subsp_/subsp\./;
        $species =~ s/_/ /g;
        $taxon_id = $tree->get_taxon_id($species);
        $snps->{$taxon_id}->{$id} = {
                                        'put_id'        =>  $put_id,
                                        'species'       =>  $species,
                                        'size'          =>  $size,
                                        'member_count'  =>  $members,
                                        'snps'     =>  [],
                                    };

    }

    if (/^\s+has no SNP/) {
        delete $snps->{$taxon_id}->{$id};
    }

    if (/^total (\d+) SNP/) {
        my $snp_count = $1;
        my $counter = 0;
        for (my $i = 0; $i < $snp_count; $i++) {
            
            $_ = <$infh>;
            chomp;
            
            if (/(\d+)\s+:\[(.)\/(.)\]\s+A:\s+(\d+)\s+C:\s+(\d+)\s+G:\s+(\d+)\s+T:\s+(\d+)/) {
                $counter++;
                my ($loc, $ref, $snp, $a, $c, $g, $t) = ($1, $2, $3, $4, $5, $6, $7);

                push(@{$snps->{$taxon_id}->{$id}->{'snps'}}, {
                                                    'loc'   =>  $loc,
                                                    'ref'   =>  $ref,
                                                    'snp'   =>  $snp,
                                                    'a'     =>  $a,
                                                    'c'     =>  $c,
                                                    'g'     =>  $g,
                                                    't'     =>  $t,
                                                             });
            } else {
                if ($snp_count != $counter) {
                    print $logfh "WARNING: $id should have $snp_count SNPs but only has data for $counter\n";
                    last;
                }
            }
                                           
        }
    }

}

## write the database file
if (-e $outdb) {
    print $logfh "Database '$outdb' exists, removing.\n";
    unlink($outdb);
}
my $db = new DBM::Deep($outdb);
$db->import($snps->{$taxon_id});
