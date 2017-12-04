#!/usr/bin/perl

use strict;
use warnings;

## genome to CAP3 assembly alignment tool
##
## Uses exonerate to align genome segments to the consensus
## sequence of a CAP3 alignment, and incorporate the genomic
## regions into the multiple sequence alignment
##
## This is to be used to produce output compatible with the
## legacy TIGR SNP-calling script
##
## Brett Whitty
## whitty@msu.edu
##
## Required input:
##
##  hit_table
##
##      A table file of EST to genomic region mappings resulting 
##      from a precompute (eg: using gmap) that takes less time 
##      to run than exonerate. Alignments done by this script
##      will be constrained to within these regions (and some
##      number of extended bases, default 5000) using exonerate.
##      This is a tab-delimited text file with at least the 
##      following columns:
##          EST_id  assembly_id  assembly_start  assembly_end
##      
##  est_fasta   
##  
##      fasta file containing the est assemblies
##
##  genome_fasta
##
##      fasta file containing the genomic sequence
##
##  cap3_alignment
##
##      the CAP3 alignment file
##


use FindBin;
use lib "$FindBin::Bin";
use CigarAligner;

use Getopt::Long;
use Cwd qw{ abs_path getcwd };
use File::Temp qw/ :POSIX /;
use String::CRC32;
use DBM::Deep;

my $consensus_temp = tmpnam();
my $genome_fasta;
my $est_fasta;
my $cap3_align;
my $hit_table;
my $depth = 4; ## depth of genomic seq to add into MSA
my $chars = 'abcdefghijklmnopqrstuvwxyz'; ## used to create unique id
my $extend = 5000;
$| = 1;

my $result = GetOptions(
    'hit_table|h=s'         =>  \$hit_table,
    'cap3_alignment|a=s'    =>  \$cap3_align,
    'genome_fasta|g=s'      =>  \$genome_fasta,
    'est_fasta|e=s'         =>  \$est_fasta,
    'depth|d=i'             =>  \$depth,
    'extend=i'              =>  \$extend,
);

$est_fasta = abs_path($est_fasta);
$genome_fasta = abs_path($genome_fasta);

unless ($hit_table && -e $hit_table) {
    die "Provide a table of genome alignment segments with --hit_table flag";
}
unless ($cap3_align && -e $cap3_align) {
    die "Provide a CAP3 alignment file with --cap3_alignment flag";
}
unless ($est_fasta && -e $est_fasta) {
    die "Provide an est fasta file with --est_fasta flag";
}
unless ($genome_fasta && -e $genome_fasta) {
    die "Provide a genome fasta file with --genome_fasta flag";
}

my $aln_seq = {};

## create a new aligner object
my $align = new CigarAligner(
    'genome_fasta'  =>  $genome_fasta,
    'est_fasta'     =>  $est_fasta,
    'extend'        => $extend, ## extend by this amount (default 5000)
);

my $con_db = make_gapped_consensus_db($cap3_align);

open my $infh, "<", $hit_table || die "Failed to open '$hit_table' for reading: $!";

while (<$infh>) {
    chomp;

    my @t = split("\t", $_);

    my $est_id = $t[0];
    my $genome_id = $t[1];
    my $genome_start = $t[2];
    my $genome_end = $t[3];

    if ($genome_end < $genome_start) {
        die "genome end position < genome start position: this is a problem";
    }

    my ($genome_segment, $cigar) = $align->do_alignment($est_id, $genome_id, $genome_start, $genome_end);

    my @cigar_atts = split(/\s+/, $cigar);

    my $score = $cigar_atts[9];

    my $consensus_seq = $con_db->seq($est_id);
   
    my @consensus_bases = split(//, $consensus_seq);
    my @genome_bases    = split(//, $genome_segment); 

    my $gapped_genome_segment = '';
    foreach my $base(@consensus_bases) {
        if ($base eq '-') {
            $gapped_genome_segment .= '-';
        } elsif (scalar(@genome_bases) == 0) {
            ## add trailing end gaps
            $gapped_genome_segment .= '-';
        } else {
            $gapped_genome_segment .= shift @genome_bases;
        }
    }
    if (@genome_bases) {
        die "There were still some bases left in the stack --- shouldn't happen";
    }
    
    my $short_id = 'g'.uc(sprintf("%x", crc32($cigar)));


    my $ggs_fasta = $gapped_genome_segment;
    $ggs_fasta =~ s/(.{1,60})/$1\n/g;

    $consensus_seq =~ s/(.{1,60})/$1\n/g;

    print STDERR ">$short_id $cigar\n$ggs_fasta";
    
    ## replace start gaps with space
    if ($gapped_genome_segment =~ /^([-]+)/) {
        my $temp_space = $1;
        $temp_space =~ tr/-/ /;
        $gapped_genome_segment =~ s/^[-]+//;
        $gapped_genome_segment = $temp_space . $gapped_genome_segment;
    }
    ## remove end gaps
    $gapped_genome_segment =~ s/[-]+$//;
    ## store sequence for mapping to CAP3 file 
    push(@{$aln_seq->{$est_id}}, [$short_id, $gapped_genome_segment, $score]);
}

foreach my $key(keys(%{$aln_seq})) {
    my $arr_ref = $aln_seq->{$key};
    if (scalar @{$arr_ref} > 1) {
        #use Data::Dumper;
        #print Dumper $arr_ref;
        @{$arr_ref} = sort {$b->[2] <=> $a->[2]} @{$arr_ref};
        #print Dumper $arr_ref;
        #die();
        @{$arr_ref} = ($arr_ref->[0]);
    }
}

#use Data::Dumper;
#print Dumper $aln_seq;



## read CAP3 alignment and add in genomic sequences
process_cap3($cap3_align);


## cleanup
unlink($consensus_temp);


## add genome sequence to CAP3 alignment
sub process_cap3 {
    my ($cap3_alignment) = @_;

#    my $depth = 4;

    open my $infh, "<", $cap3_alignment || die "Failed to open CAP3 alignment '$cap3_alignment'";

    my $align_flag = 0;
    my $contig_id = '';
    my $consensus = '';
    my @genome_seqs = ();
    while (<$infh>) {
        print;
        if (/^DETAILED/) {
            $align_flag = 1;
            next;
        }
        if (/\*+ (\S+) \*+/ && ! $align_flag) {
            my $temp_id = $1;
            @genome_seqs = ();
            if (defined($aln_seq->{$temp_id})) {
                @genome_seqs = @{$aln_seq->{$temp_id}};
                foreach my $seq_ref(@genome_seqs) {
                    for (my $i = 0; $i < $depth; $i++) {
                        my $seq_id = $seq_ref->[0] . substr($chars, $i, 1);
                        print $seq_id."\n";
                    }
                }
            }
        }
        if (! $align_flag) {
            next;
        }
        if (/\*+ (\S+) \*+/) {
            $contig_id = $1;
            @genome_seqs = ();
            if (defined($aln_seq->{$contig_id})) {
                @genome_seqs = @{$aln_seq->{$contig_id}};
            }
        } elsif (/^Number of/) {
            $align_flag = 0;
            next;
        } elsif (/^\s+\.\s{4}:/) {
            # print STDERR "Hit alignment for $contig_id\n";
            if (scalar(@genome_seqs) < 1) {
                next;
            }

            foreach my $genome_seq(@genome_seqs) {
                for (my $i = 0; $i < $depth; $i++) { 
                    my $seq_id = $genome_seq->[0] . substr($chars, $i, 1);
                    my $seq_line = $seq_id 
                        . ' ' x (22 - length($seq_id))
                        . substr($genome_seq->[1], 0, 60)
                        ."\n";
                
                    print $seq_line;
                }

                if (length($genome_seq->[1]) > 60) {  
                    $genome_seq->[1] = substr($genome_seq->[1], 60);
                }
            }
        }
    }
}


## extra consensus sequence with gaps from the CAP3 alignment
sub make_gapped_consensus_db {
    my ($cap3_alignment) = @_;

    open my $infh, "<", $cap3_alignment || die "Failed to open CAP3 alignment '$cap3_alignment'";
    open my $outfh, ">", $consensus_temp || die "Failed to open consensus temp '$consensus_temp' for writing";

    my $align_flag = 0;
    my $contig_id = '';
    my $consensus = '';
    while (<$infh>) {
        chomp;
        if (/^DETAILED/) {
            $align_flag = 1;
            next;
        }
        if (! $align_flag) {
            next;
        }
        if (/\*+ (\S+) \*+/) {
            my $id = $1;
            if ($consensus) {
                $consensus =~ s/(.{1,60})/$1\n/g;
                print $outfh ">$contig_id\n$consensus";
                $consensus = '';
            }
            $contig_id = $id;
        } elsif (/^Number of/) {
            $align_flag = 0;
            next;
        } elsif (/^consensus\s+(\S+)/) {
            $consensus .= $1;
        }
    }
    if ($consensus) {
        $consensus =~ s/(.{1,60})/$1\n/g;
        print $outfh ">$contig_id\n$consensus";
    }

    close $outfh;

    my $db = new Bio::DB::Fasta($consensus_temp);

    return $db;
}
