#!/usr/bin/perl

use strict;
use warnings;

package SynNon;

## 
## A class for comparing codons to determine whether they are synonymous or non-synonymous.
##
## 
## Brett Whitty
## whitty@msu.edu
##

use strict;
use warnings;
use Carp;

my $aa;

## support ambiguous base codes
my $base_codes = {
    'A' =>  ['A'],
    'C' =>  ['C'],
    'G' =>  ['G'],
    'T' =>  ['T'],
    'X' =>  ['X'],
    'N' =>  ['N'],
    'K' =>  [
        'G',
        'T',
    ],
    'M' =>  [
        'A',
        'C',
    ],
    'R' =>  [
        'A',
        'G',
    ],
    'S' =>  [
        'C',
        'G',
    ],
    'W' =>  [
        'A',
        'T',
    ],
    'Y' =>  [
        'C',
        'T',
    ],
};

sub new {

    my ($this, %args) = @_;

    my $class = ref($this) || $this;
    
    my $self = {};

    my $translation_table = $args{'translation_table'} || 0;
    
    ## initialize the codon lookup hash
    my $codon_table_ref = _initialize_codon_table($translation_table);
    
    ## ref to the aa translation hash
    $aa = $codon_table_ref->{'aa'};

    bless $self, $class;

    return $self;
}


## returns a floating point value between 0 and 1 related to the truthiness
## of whether two codons are synonymous
##
## NB: since ambiguous codons are supported, two codons can be partially synonymous
## so using that has implications for using the return value of this method 
## in a boolean context
sub is_synonymous {
    my ($self, $ref_codon, $snp_codon) = @_;

    my $result = $self->compare_codons($ref_codon, $snp_codon);

    return $result->[0];
}

## translates and compares two codons, returns results in an array reference
sub compare_codons {

    my ($self, $ref_codon, $snp_codon) = @_;

    my @ref_codons = _translate_ambiguous_codon(uc($ref_codon));
    my $ref_ambiguous = (scalar(@ref_codons) > 1) ? 1 : 0;
    my @snp_codons = _translate_ambiguous_codon(uc($snp_codon));
    my $snp_ambiguous = (scalar(@snp_codons) > 1) ? 1 : 0;

    my $comps = 0;
    my $syn = 0;

    my @ref_aas = ();
    my @snp_aas = ();
    foreach my $aref_codon(@ref_codons) {
        foreach my $asnp_codon(@snp_codons) {
            $comps++;
            
            my $ref_aa = $aa->{$aref_codon} || '?';
            my $snp_aa = $aa->{$asnp_codon} || '?';
            
            push(@ref_aas, $ref_aa);
            push(@snp_aas, $snp_aa);

            if ($ref_aa eq '?') {
                carp "Reference codon '$aref_codon' not in translation table";
            }
            if ($snp_aa eq '?') {
                carp "SNP codon '$asnp_codon' not in translation table";
            }

            if ($ref_aa eq $snp_aa) {
                $syn++;
            }
        }
    }

    return [ sprintf("%.1f", $syn/$comps), \@ref_aas, \@snp_aas, $ref_ambiguous, $snp_ambiguous ];
}

## translates a codon with possible ambiguity codes to one or more codons
sub _translate_ambiguous_codon {

    my ($codon) = @_;

    $codon = uc($codon);

    my $bases = [];
    for (my $i = 0; $i < 3; $i++) {
        my $base = substr($codon, $i, 1);
        my @codes = (defined($base_codes->{$base})) ?
            @{$base_codes->{$base}} : ( 'N' );
        $bases->[$i] = [ @codes ] ; 
    }
    my $codons = {};
    foreach my $b1(@{$bases->[0]}) {
        foreach my $b2(@{$bases->[1]}) {
            foreach my $b3(@{$bases->[2]}) {
                $codons->{$b1.$b2.$b3} = 1;
            }
        }
    };
    my @codons = keys(%{$codons});

    #use Data::Dumper;
    #print Dumper $bases;
    #die();

    return @codons;
}

## read a codon table and initialize a hash structure for storing it
sub _initialize_codon_table {
    my ($table) = @_;

    my $table_ref = [];
    my $aa_ref = {};
    my $start_ref = {};
    my $stop_ref = {};
    my $table_string = '';
    
    my $parsed_code;
 
    my %supported_tables;
    
    ## initialize a flag hash for the internally supported codon tables
    foreach my $table_id ((0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)) {
        $supported_tables{$table_id} = 1;
    } 
    
    ## table 0 
    $table_ref->[0] = [split('', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')];
    $table_ref->[1] = [split('', '-----------------------------------M----------------------------')];
    $table_ref->[2] = [split('', 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG')];
    $table_ref->[3] = [split('', 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG')];
    $table_ref->[4] = [split('', 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG')];
     
    ## temp hash for finding degenerate codons
    my $ncodons_ref = {};
    while (scalar(@{$table_ref->[0]})) {
        my $aa    = shift(@{$table_ref->[0]});
        my $start = shift(@{$table_ref->[1]});
        my $base1 = shift(@{$table_ref->[2]}); 
        my $base2 = shift(@{$table_ref->[3]}); 
        my $base3 = shift(@{$table_ref->[4]}); 
        my $codon = $base1 . $base2 . $base3;

        $aa_ref->{$codon} = $aa;
        
        $ncodons_ref->{$base1.$base2.'N'}->{$aa} = 1;
        $ncodons_ref->{$base1.'N'.$base3}->{$aa} = 1;
        $ncodons_ref->{'N'.$base2.$base3}->{$aa} = 1;
        $ncodons_ref->{$base1.'N'.'N'}->{$aa}    = 1;
        
        ## store if not '-' --- may want to put other conditions on here
        if ($start ne '-') {
            $start_ref->{$codon} = $start;
        }

        ## might as well store a hash of stop codons as 'stop'
        if ($aa eq '*') {
            $stop_ref->{$codon} = $aa;
        }
        
    }
   
    ## add codons for unknown bases
    foreach my $codon(keys(%{$ncodons_ref})) {
        
#        ## if it could be a stop, we'll assume it is
#        if ($options{'assume_stops'}) {
#            if ($ncodons_ref->{$codon}->{'*'}) {
#                $stop_ref->{$codon} = '*';
#            }
#        }
        
        ## @aas = the set of aas coded by the ambiguous codon
        my @aas = keys(%{$ncodons_ref->{$codon}});
        
        ## add codon lookups for unknowns at degenerate positions
        if (scalar(@aas) == 1) {
            $aa_ref->{$codon} = shift(@aas);
        }
    }

    ## NNN is to be treated as a hard stop
    $stop_ref->{'NNN'} = '*';
     
    return { 'aa' => $aa_ref, 'start' => $start_ref, 'stop' => $stop_ref, 'table' => $table_string};
}

1;


##################################################### MAIN #########################################


package main;

## Given a table of SNP positions and genome GFF and fasta
## will output a table of reference and SNP codons
##
## Brett Whitty
## whitty@msu.edu

use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use FindBin;
use POSIX;

my ($input, $output, $genome_gff, $genome_fasta, $source, $process_org, $process_ref, $ignore_errors);

GetOptions(
    'input|i=s'     =>  \$input,
    'output|o=s'    =>  \$output,
    'gff|g=s'       =>  \$genome_gff,
    'fasta|f=s'     =>  \$genome_fasta,
    'source|s=s'    =>  \$source,
    'org=s'         =>  \$process_org,
    'ref=s'         =>  \$process_ref,
    'errors|e!'     =>  \$ignore_errors,
);

$source ||= 'MSU';

my $synnon = new SynNon();

## open for writing
my $outfh;
if (defined($output)) {
    open $outfh, '>', $output || croak "Failed to open '$output' for writing: $!";
} else {
    $outfh = \*STDOUT;
}

if (! -e $genome_gff) {
    croak "Unable to access genome GFF!";
}
my $db = new Bio::DB::Fasta($genome_fasta);

if (! $db) {
    croak "Unable to access genome fasta database!";
}


## these globals are going to store all the GFF relationship data
my $gene_coords;
my $cds_coords;
my $feat_types;
my $parents;
my $children;


my $infh;
if (defined($input) && -e $input) {
    ## open the input table for reading
    open $infh, '<', $input || croak "Failed to open SNP table file '$input' for reading: $!";
} else {
    $infh = \*STDIN;
}

## parse the GFF to get the coordinates files
parse_gff($genome_gff);

my %id_counter = ();
while (<$infh>) {
    chomp;

    ## skip comment lines
    if (/^#/) { next; }

    my @t = split(/\t/, $_);

    ## reading in the SNP table here
    my ($org, $ref_genome, $snp_genome, $ref_seq_id, $snp_loc, $ref_base, $snp_base) = @t;

    ## only process the organism specified by the user
    if (defined($process_org) && $org ne $process_org) { next; }
    
    ## use gene coordinates to quickly find which genes SNP affects
    my @gene_refs = do_snp_gene_scan($ref_seq_id, $snp_loc);

    ## iterate through any genes the SNP hits
    foreach my $gene_ref(@gene_refs) {
        
        ## get all CDS children of gene features
        my @cdses = get_children($gene_ref->[2], 'CDS');

        ## find what region of each CDS is affected by the SNP
        foreach my $cds_id(@cdses) {
            my %atts = do_snp_cds_scan($ref_seq_id, $snp_loc, $ref_base, $snp_base, $cds_id, $gene_ref);

            ## store the origin of the SNP in the atts
            $atts{'org'} = $org;
            $atts{'ref_genome'} = $ref_genome;
            $atts{'snp_genome'} = $snp_genome;

            ## set CDS as Parent
            $atts{'Parent'} = $cds_id;

            $atts{'ref_base'} = $ref_base;
            $atts{'snp_base'} = $snp_base;

            ## if we have a codon
            if (defined($atts{'snp_codon'})) {
                my $result = $synnon->compare_codons($atts{'ref_codon'}, $atts{'snp_codon'});
                $atts{'synonymous'} = $result->[0];
                $atts{'ref_aa'} = join(',', @{$result->[1]});
                $atts{'snp_aa'} = join(',', @{$result->[2]});
                $atts{'ref_ambiguous'} = $result->[3];
                $atts{'snp_ambiguous'} = $result->[4];
            }

            ## output the result
            my @atts = ();
            push (@atts, 'ID='.$ref_seq_id.'.'.$snp_loc.'.'.++$id_counter{$ref_seq_id.'.'.$snp_loc});
            foreach my $key(sort(keys(%atts))) {
                push(@atts, $key.'='.$atts{$key});
            }
            my $atts = join(';', @atts);
            print join("\t", (
                $ref_seq_id,
                $source,
                'SNP',
                $snp_loc,
                $snp_loc,
                '.',
                '+', ## SNP itself will always be on +ve strand
                '.',
                $atts,
            ))."\n";
        }
    }
}

## scans for genes containing the SNP
sub do_snp_gene_scan {
    my ($seq_id, $loc) = @_;

    my @gene_refs = ();

    foreach my $gene_ref(@{$gene_coords->{$seq_id}}) {
        if ($loc < $gene_ref->[0]) {
            ## SNP is upstream of the current (and future) genes
            last;
        } elsif ($loc > $gene_ref->[1]) {
            ## SNP coord is beyond the end of this gene
            next;
        } elsif ($loc >= $gene_ref->[0]) {
            ## SNP is inside this gene
            push(@gene_refs, $gene_ref);
        }
    }

    return @gene_refs;
}

sub parse_gff {
    my ($gff) = @_;

    ## store the CDS coords
    $cds_coords = {};
    $gene_coords = {};

    ## these will consume a lot of memory but in case we want to dump out corresponding
    ## transcript / gene information I'm going to keep them for now
    $feat_types = {};
    $parents = {};
    $children = {};

    open my $infh, '<', $gff || croak "Failed to open GFF file '$gff' for reading: $!";

    my $cds = {};
    while (<$infh>) {
        chomp;

        my @t = split("\t", $_);

        ## only consider CDS feature GFF lines
        if (scalar(@t) != 9) {
            next;
        }
        
        my $feat_type = $t[2];
        my $feat_id = '';

        my $att_string = $t[8];
        my @atts = split(/;/, $att_string);
        foreach my $att(@atts) {

            my ($key, $value_string) = split(/=/, $att);
            my @values = split(/,/, $value_string);

            if ($key eq 'ID') {
                $feat_id = $values[0];
                ## this should never happen (NB: CDSs can span multiple lines in GFF3)
                if (defined($feat_types->{$feat_id}) && $feat_type !~ /^CDS$/i) {
                    croak "Refusing to stomp on feat_type hash for '$feat_id'";
                }
                $feat_types->{$feat_id} = $feat_type;
            } elsif ($key eq 'Parent') {
                foreach my $parent_id(@values) {
                    $parents->{$feat_id}->{$parent_id} = 1;
                    $children->{$parent_id}->{$feat_id} = 1;
                }
            }
        }
        
        my ($seq_id, $start, $end, $strand, $phase) = ($t[0], $t[3], $t[4], $t[6], $t[7]);

        if ($feat_type =~ /^CDS$/i) {
            ## initially store CDS segments in a hash by CDS ID
            push(@{$cds_coords->{$seq_id}->{$feat_id}}, [
                $start,
                $end,
                $strand,
                $phase,
            ]);
        } elsif ($feat_type =~ /^gene$/i) {
            ## store only coordinates and feature ID of genes for  
            push(@{$gene_coords->{$seq_id}}, [
                $start,
                $end,
                $feat_id,
            ]);
       }
    }
 
    ## sort the gene coordinate hash 
    foreach my $seq_id(keys(%{$gene_coords})) {
        @{$gene_coords->{$seq_id}} = sort { $a->[0] <=> $b->[0] } @{$gene_coords->{$seq_id}};
    }

    ## sort the cds coordinate hash 
    foreach my $seq_id(keys(%{$cds_coords})) {
        foreach my $cds_id(keys(%{$cds_coords->{$seq_id}})) {
            @{$cds_coords->{$seq_id}->{$cds_id}} = sort { $a->[0] <=> $b->[0] } @{$cds_coords->{$seq_id}->{$cds_id}};
        }
    }

}

## scans a CDS for where a SNP falls and categorizes it
sub do_snp_cds_scan {
    my ($ref_seq_id, $snp_loc, $ref_base, $snp_base, $cds_id, $gene_ref) = @_;

    ## to store output attributes
    my %atts = ();

    ## gene coordinates
    my $gene_start = $gene_ref->[0];
    my $gene_end   = $gene_ref->[1];

    ## reference to an array of CDSes
    my $cds_ref = $cds_coords->{$ref_seq_id}->{$cds_id};

    my $loc_code = '';
    my $cds_start = -1;
    my $cds_end = -1;
    my $cds_strand = '';
    my $cds_phase = -1;

    my $cds_count = scalar(@{$cds_ref});
    ## iterate through all CDS regions     
    for (my $i = 0; $i < $cds_count; $i++) {

        ## for convenience
        $cds_start = $cds_ref->[$i]->[0];
        $cds_end = $cds_ref->[$i]->[1];
        $cds_strand = $cds_ref->[$i]->[2];
        $cds_phase = $cds_ref->[$i]->[3];

        ## sequential CDS number
        my $cds_number = ($cds_strand eq '+') ? $i + 1 : $cds_count - $i;

        ## sequential intron number (meaningless when =0 or =cds_count)
        my $intron_number = ($cds_strand eq '+') ? $i : $cds_count - $i;
        
        ## if the SNP is before the first CDS region
        if ($snp_loc < $cds_start) {

            if ($i == 0) {
                ## we're in the upstream UTR region

                ## F = 5' UTR / T = 3' UTR
                $loc_code = ($cds_strand eq '+') ? 'F' : 'T';

                ## calculate the offset value from the start of the gene
                my $pos_offset = calc_offset($gene_start, $cds_start, $snp_loc, $cds_strand);
                $atts{'region_offset'} = $pos_offset;
            } else {
                ## SNP is in an intron

                ## I = intron
                $loc_code = 'I' . ($intron_number);

                ## calculate offset between previous CDS end and current CDS start
                my $pos_offset = calc_offset($cds_ref->[$i - 1]->[1], $cds_start, $snp_loc, $cds_strand);
                $atts{'region_offset'} = $pos_offset;
            }

            last;

        } elsif ($snp_loc > $cds_end) {
            ## past end of CDS so on to the next one
            next;
        } elsif ($snp_loc <= $cds_end) {
            ## OK, we are in a CDS so we need to get the codon
            $loc_code = 'C'.$cds_number;
            
            ## calculate the offset from the CDS boundaries
            my $pos_offset = calc_offset($cds_start, $cds_end, $snp_loc, $cds_strand);
            $atts{'region_offset'} = $pos_offset;


            ## calculate the position of the SNP in relative coordinates within the CDS
            ## and polypeptide sequences
            ($atts{'pos_in_cds'}, $atts{'pos_in_pep'}) = calc_cds_pep_pos($loc_code, $pos_offset, $cds_strand, $cds_ref); 


            my $codon_start = -1;
            if ($cds_strand eq '+') {
                ## calculate codon start position on forward strand
                $codon_start = $snp_loc - ($snp_loc - ($cds_start + $cds_phase)) % 3;
            } else {
                ## calculate codon start position on the reverse strand
                $codon_start = $snp_loc + (($cds_end - $cds_phase) - $snp_loc) % 3;
            }

            ## set 'codon_start' attribute
            $atts{'codon_start'} = $codon_start;

            ## OK, now get each base individually
            my $base_loc = -1;
            my $base_position = -1;
            my $ref_codon = '';
            my $snp_codon = '';

            for (my $base_index = 0; $base_index < 3; $base_index++) {

                my $missing_base = 0;

                if ($cds_strand eq '+') {
                    ## forward strand
                    $base_position = $codon_start + $base_index;
                } else {
                    ## reverse strand
                    $base_position = $codon_start - $base_index;
                }
                    
                ## adjust base position if it crosses CDS segment region boundaries
                if ($base_position < $cds_start) {
                    ## base position beyond the CDS start
                    
                    if (! defined($cds_ref->[$i - 1])) {
                        ## in the case of partial CDSes we might be missing a CDS region
                        $missing_base = 1;
                    } else {
                        $base_loc = $cds_ref->[$i - 1]->[1] - ($cds_start - $base_position - 1);
                    }

                } elsif ($base_position > $cds_end) {
                    ## base position beyond the CDS end
                    
                    #   print "DEBUG: cds_ref->[i + 1]->[1]=".$cds_ref->[$i - 1]->[1]."  cds_start=$cds_start  base_position=$base_position\n";
                    if (! defined($cds_ref->[$i + 1])) {
                        use Data::Dumper;
                        print "$i + 1\n";
                        $missing_base = 1;
                    } else {
                        $base_loc = $cds_ref->[$i + 1]->[0] + ($base_position - $cds_end - 1);
                    }
                } else {
                    ## base position is just fine
                    $base_loc = $base_position;
                }
                
                my $fetched_base;

                if (! $missing_base) {
                    ## fetch the base
                    $fetched_base = uc($db->seq($ref_seq_id, $base_position, $base_position));
                } else {
                    $fetched_base = 'N';
                }

                ## lower-case the SNP base
                if ($base_position == $snp_loc) {
                    ## complain if there's a reference base mismatch between input table and fasta file
                    if (uc($ref_base) ne $fetched_base) {
                        carp "Expected reference base '$ref_base' at $ref_seq_id:$snp_loc but found '$fetched_base'";
                        if (! $ignore_errors) {
                            confess "Refusing to ignore ref base in genome fasta != ref base in SNP table, override with --errors flag";
                        }
                    }
                    $fetched_base = lc($fetched_base);
                    $snp_codon .= lc($snp_base);
                } else {                
                    $snp_codon .= $fetched_base;
                }
                $ref_codon .= $fetched_base;
            }

            ## complement the codon if its on the negative strand
            ## (it was already reversed when we built it base by base)
            if ($cds_strand eq '-') {
                $ref_codon = complement_dna($ref_codon);
                $snp_codon = complement_dna($snp_codon);
            }

            $atts{'ref_codon'} = $ref_codon;
            $atts{'snp_codon'} = $snp_codon;
            $atts{'codon_strand'} = $cds_strand;

            last;
        }
    }

    ## if we iterated through all CDSes without hitting the interval containing the SNP
    if ($loc_code eq '') {
        ## so we're in downstream UTR
        $loc_code = ($cds_strand eq '+') ? 'T' : 'F';

        ## calculate offset base on end coordinate of gene region
        my $pos_offset = calc_offset($cds_end, $gene_end, $snp_loc, $cds_strand);
        $atts{'region_offset'} = $pos_offset;
    }
    $atts{'region_code'} = $loc_code;

#    print $atts."\n";
    return %atts;
}

## calculates the offset of the SNP from the nearest end of the feature (or region)
## -ve means bases from the left side of the feature
## +ve means bases from the right side of the feature
sub calc_offset {
    my ($start, $end, $snp, $strand) = @_;

    my $left_offset = $snp - $start;
    my $right_offset = $end - $snp;

    my $feature_offset;
    if ($left_offset < $right_offset) {
        $feature_offset = -($left_offset);
    } else {
        $feature_offset = $right_offset;
    }

    if ($strand eq '-') {
        $feature_offset = -($feature_offset);
    }

    return $feature_offset;
}


## uses the already determined location code and offset to calculate
## the postition of the SNP within the CDS and polypeptide sequence records
## (in their relative coordinates)
sub calc_cds_pep_pos {
    my ($loc_code, $offset, $strand, $cds_segs) = @_;

    ## base position within the CDS sequence record
    my $cds_pos;
    ## aa position within the polypeptide sequence record
    my $pep_pos;

    ## we expect that this will only be run for CDSes
    $loc_code =~ /^C(\d+)/ or confess "calc_cds_pep_pos called with bad params";
    my $cds_index = $1 - 1;

    my $cdses = [ @{$cds_segs} ];
    my $cds_count = @{$cdses};
    
    ## we'll just flip the order of the CDS segments if we're dealing with the -ve strand
    if ($strand eq '-') {
        @{$cdses} = reverse @{$cdses};
    }

    $cds_pos = 0;
    ## iterate through CDS segments until we reach the one the SNP occurs in
    for (my $i = 0; $i <= $cds_index; $i++) {

        if ($i == $cds_index) {
            ## if we're in the CDS with the SNP, we need to do special math

            if ($offset < 0) {
                ## calc # of bases if the offset is from the left side of the current CDS segment
                $cds_pos += -1 * $offset + 1;
            } else {
                ## calc # of bases if the offset is from the right side of the current CDS segment
                $cds_pos += $cdses->[$i]->[1] - $cdses->[$i]->[0]  - $offset + 1;
            }

        } else {
            ## otherwise we're going to add the length of this whole CDS segment to the total
            $cds_pos += $cdses->[$i]->[1] - $cdses->[$i]->[0] + 1;
        }
    }

    ## calc which aa within the polypeptide is affected by the SNP
    $pep_pos = ceil($cds_pos / 3);

    return ($cds_pos, $pep_pos);
}


## gets all children of a certain type for a feature
sub get_children {
    my ($parent_id, $child_type) = @_;

    my $found_children;
    my @to_search = ($parent_id);

    while (@to_search) {
        my $search_id = shift @to_search;

        ## if a feature is the type we're looking for, push it to found children
        if ($feat_types->{$search_id} eq $child_type) {
            $found_children->{$search_id} = 1;
        ## otherwise we need to search its children
        } else {
            push (@to_search, keys(%{$children->{$search_id}}));
        }
    }

    return keys(%{$found_children});
}


sub complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;

    return $r_seq;
}

__END__

## some codon tables are attached here
__DATA__
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differs from Genetic Code [1] only in that the initiation sites have been
# changed to only 'AUG'.
#
# This is intended to be the default genetic code where the use of
# the rare initiation sites (CUG, UUG) is not intended.

Genetic Code [0]

Standard with AUG start only
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.4
#     Added CTG,TTG as allowed alternate start codons in Standard code.
#        Prats et al. 1989, Hann et al. 1992
#
# Initiation Codon:
#
# AUG 
# 
# Alternative Initiation Codons
#
# In rare cases, translation in eukaryotes can be initiated from codons
# other than AUG.  A well documented case (including direct protein
# sequencing) is the GUG start of a ribosomal P protein of the fungus
# Candida albicans (Abramczyk et al.).  Other examples can be found in the
# following references: Peabody 1989; Prats et al.  1989; Hann et al. 
# 1992; Sugihara et al.  1990. 
#
# GUG, CUG, UUG 

Genetic Code [1]

Standard
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 2          Standard
# 
#  AGA    Ter  *          Arg  R
#  AGG    Ter  *          Arg  R
#  AUA    Met  M          Ile  I
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codon:
# 
# Bos: AUA 
# Homo: AUA, AUU
# Mus: AUA, AUU, AUC
# Coturnix, Gallus: also GUG (Desjardins and Morais, 1991)
# 
# Systematic Range:
# 
# Vertebrata
# 
# Comment: 
# 
# The transcripts of several vertebrate mitochondrial genes end in U or
# UA, which become termination codons (UAA) upon subsequent
# polyadenylation. 

Genetic Code [2]
 
Vertebrate Mitochondrial

AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
Starts = --------------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 3          Standard
# 
#  AUA    Met  M          Ile  I
#  CUU    Thr  T          Leu  L
#  CUC    Thr  T          Leu  L
#  CUA    Thr  T          Leu  L
#  CUG    Thr  T          Leu  L
#  UGA    Trp  W          Ter  *
# 
#  CGA    absent          Arg  R
#  CGC    absent          Arg  R
# 
# 
# Systematic Range: 
# 
# Saccharomyces cerevisiae, Candida glabrata, Hansenula saturnus,
# and Kluyveromyces thermotolerans
# (Clark-Walker and Weiller, 1994)
# 
# Comments:
# 
# The remaining CGN codons are rare in Saccharomyces cerevisiae and
# absent in Candida glabrata (= Torulopsis glabrata). 
#
# The AUA codon is common in the gene var1 coding for the single
# mitochonLIial ribosomal protein, but rare in genes encoding the enzymes. 
#
# The coding assignments of the AUA (Met or Ile) and CUU (possibly Leu,
# not Thr) are uncertain in Hansenula saturnus. 
#
# The coding assignment of Thr to CUN is uncertain in Kluyveromyces
# thermotolerans (Clark-Walker and Weiller, 1994). 

Genetic Code [3]
 
Yeast Mitochondrial
 
AAs  =   FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ----------------------------------MM----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 4         Standard
# 
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codons: 
# 
# Trypanosoma: UUA, UUG, CUG
# Leishmania: AUU, AUA 
# Tertrahymena: AUU, AUA, AUG 
# Paramecium: AUU, AUA, AUG, AUC, GUG, GUA(?) 
# (Pritchard et al., 1990)
# 
# Systematic Range: 
# 
# Mycoplasmatales: Mycoplasma, Spiroplasma (Bove et al., 1989); 
# 
# Fungi: Emericella nidulans, Neurospora crassa, Podospora anserina,
# Acremonium (Fox, 1987), Candida parapsilosis (Guelin et al., 1991),
# Trichophyton rubrum (de Bievre and Dujon, 1992), Dekkera/Brettanomyces,
# Eeniella (Hoeben et al., 1993), and probably Ascobolus immersus,
# Aspergillus amstelodami, Claviceps purpurea, and Cochliobolus
# heterostrophus. 
# 
# Protozoa: Trypanosoma brucei, Leishmania tarentolae, Paramecium
# tetraurelia, Tetrahymena pyriformis and probably Plasmodium gallinaceum
# (Aldritt et al., 1989)]. 
# 
# Metazoa: Coelenterata (Ctenophora and Cnidaria) 
# 
# Comments: 
# 
# This code is also used for the kinetoplast DNA (maxicircles,
# minicircles).  Kinetoplasts are modified mitochondria (or their parts). 
# 
# This code is not used in the Acholeplasmataceae and plant-pathogenic
# mycoplasma-like organisms (MLO) (Lim and Sears, 1992)

Genetic Code [4]
  
Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --MM---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.3 - 10/13/95
#     Added alternate intiation codon ATC to code 5
#        based on complete mitochondrial genome of honeybee
#        Crozier and Crozier (1993)
#
#  Version 3.2 - 6/24/95
#     GTG allowed as alternate initiator
#
# Comment: 
# 
# The codon AGG is absent in Drosophila. 
# 
# Differences from the Standard Code: 
# 
# 
#         Code 5          Standard
# 
#  AGA    Ser  S          Arg  R
#  AGG    Ser  S          Arg  R
#  AUA    Met  M          Ile  I
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codons: 
# 
# AUA, AUU
# AUC: Apis (Crozier and Crozier, 1993)
# GUG: Polyplacophora (Boore and Brown, 1994 GenBank Accession Number: U09810)
# UUG: Ascaris, Caenorhabditis
# 
# Systematic Range:
# 
# Nematoda: Ascaris, Caenorhabditis;
# Mollusca: Bivalvia (Hoffmann et al., 1992); Polyplacophora (Boore and
# Brown, 1994)
# Arthropoda/Crustacea: Artemia (Batuecas et al., 1988);
# Arthropoda/Insecta: Drosophila [Locusta migratoria (migratory locust),
# Apis mellifera (honeybee)]
# 
# Comments: 
# 
# GUG may possibly function as an initiator in Drosophila (Clary and
# Wolstenholme, 1985; Gadaleta et al., 1988).  AUU is not used as an
# initiator in Mytilus (Hoffmann et al., 1992). 
# 
# "An exceptional mechanism must operate for initiation of translation of
# the cytochrome oxidase subunit I mRNA in both D.  melanogaster (de
# Bruijn, 1983) and D.  yakuba (Clary and Wolstenholme 1983), since its
# only plausible initiation codon, AUA, is out of frame with the rest of
# the gene.  Initiation appears to require the "reading" of of an AUAA
# quadruplet, which would be equivalent to initiation at AUA followed
# immediately by a specific ribosomal frameshift.  Another possible
# mechanism ...  is that the mRNA is "edited" to bring the AUA initiation
# into frame." (Fox, 1987)

Genetic Code [5]
    
Invertebrate Mitochondrial
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
Starts = ---M----------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 6       Standard
# 
#  UAA      Gln  Q        Ter  *
#  UAG      Gln  Q        Ter  *
# 
# 
# Systematic Range: 
# 
# Ciliata: Oxytricha and Stylonychia (Hoffman et al.  1995), Paramecium,
# Tetrahymena, Oxytrichidae and probably Glaucoma chattoni. 
# 
# Dasycladaceae: Acetabularia (Schneider et al., 1989) and Batophora
# (Schneider and de Groot, 1991). 
# 
# Diplomonadida: 
# Scope: Hexamita inflata, Diplomonadida ATCC50330, and ATCC50380. 
# Ref.: Keeling, P.J.  and Doolittle, W.F.  1996.  A non-canonical genetic
# code in an early diverging eukaryotic lineage.  The EMBO Journal 15,
# 2285-2290. 
# 
# Comment: 
# 
# The ciliate macronuclear code has not been determined completely.  The
# codon UAA is known to code for Gln only in the Oxytrichidae. 

Genetic Code [6]
 
Ciliate, Dasycladacean and Hexamita Nuclear
    
AAs  =   FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 
# Genetic Code Table 
# 
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# 
#  Version 3.8
#     Added GTG start to Echinoderm mitochondrial code, code 9
#
# Differences from the Standard Code: 
# 
# 
#           Code 9        Standard
# 
#  AAA      Asn  N        Lys K
#  AGA      Ser  S        Arg R
#  AGG      Ser  S        Arg R
#  UGA      Trp  W        Ter *
# 
# 
# Systematic Range: 
# 
# Asterozoa (starfishes) (Himeno et al., 1987) Echinozoa (sea urchins)
# (Jacobs et al., 1988; Cantatore et al., 1989)

Genetic Code [9]
   
Echinoderm and Flatworm Mitochondrial
       
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.2 - 6/24/95
#     Alternative Ciliate Macronuclear renamed to Euplotid Macronuclear
#
# Differences from the Standard Code:
# 	Code 10     Standard
# UGA 	Cys  C        Ter  *
# 
# Systematic Range: 
# Ciliata: Euplotidae (Hoffman et al. 1995). 

Genetic Code [10]
   
Euplotid Nuclear
    
AAs  =   FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
  
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.1 - 1995
#     Addition of Eubacterial by J.Ostell at NCBI
#  Version 3.2 - 6/24/95
#     Eubacterial renamed to Bacterial as most alternate starts
#                 have been found in Archaea
#
# Differences from the Standard Code: 
# 
# None 
# 
# Alternative Initiation Codons: 
# 
# GUG, UUG, AUU, CUG 
# 
# Systematic Range and Comments: 
# 
# Table 11 is used for Bacteria, Archaea, prokaryotic viruses and
# chloroplast proteins.  As in the standard code, initiation is most
# efficient at AUG.  In addition, GUG and UUG starts are documented in
# Archaea and Bacteria (Kozak 1983, Fotheringham et al.  1986, Golderer et
# al.  1995, Nolling et al.  1995, Sazuka & Ohara 1996, Genser et al. 
# 1998, Wang et al.  2003).  In E.  coli, UUG is estimated to serve as
# initiator for about 3% of the bacterium's proteins (Blattner et al. 
# 1997).  CUG is known to function as an initiator for one plasmid-encoded
# protein (RepA) in Escherichia coli (Spiers and Bergquist, 1992).  In
# addition to the NUG initiations, in rare cases Bacteria can initiate
# translation from an AUU codon as e.g.  in the case of poly(A) polymerase
# PcnB and the InfC gene that codes for translation initiation factor IF3
# (Polard et al.  1991, Liveris et al.  1993, Sazuka & Ohara 1996, Binns &
# Masters 2002).  The internal assignments are the same as in the standard
# code though UGA codes at low efficiency for Trp in Bacillus subtilis
# and, presumably, in Escherichia coli (Hatfiled and Diamond, 1993). 

Genetic Code [11]

Bacterial and Plant Plastid
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#      Addition of Alternative Yeast by J.Ostell at NCBI
#
# Differences from the Standard Code: 
# 
#            Code 12      Standard
# 
#  CUG       Ser          Leu
#        
# 
# Alternative Initiation Codons: 
# 
# CAG may be used in Candida albicans (Santos et al., 1993). 
# 
# Systematic Range: 
# 
# Endomycetales (yeasts): Candida albicans, Candida cylindracea, Candida
# melibiosica, Candida parapsilosis, and Candida rugosa (Ohama et al.,
# 1993). 
# 
# Comment: 
# 
# However, other yeast, including Saccharomyces cerevisiae, Candida azyma,
# Candida diversa, Candida magnoliae, Candida rugopelliculosa, Yarrowia
# lipolytica, and Zygoascus hellenicus, definitely use the standard
# (nuclear) code (Ohama et al., 1993). 


Genetic Code [12]
  
Alternative Yeast Nuclear

AAs  =   FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -------------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 13     Standard
# 
#  AGA      Gly  G        Arg  R
#  AGG      Gly  G        Arg  R
#  AUA      Met  M        Ile  I
#  UGA      Trp  W        Ter  *
# 
# 
# Systematic Range and Comments: 
# 
# There is evidence from a phylogenetically diverse sample of tunicates
# (Urochordata) that AGA and AGG code for glycine.  In other organisms,
# AGA/AGG code for either arginine or serine and in vertebrate
# mitochondria they code a STOP.  Evidence for glycine translation of
# AGA/AGG has been found in Pyura stolonifera (Durrheim et al.  1993),
# Halocynthia roretzi (Kondow et al.  1999, Yokobori et al., 1993,
# Yokobori et al.  1999) and Ciona savignyi (Yokobori et al.  2003).  In
# addition, the Halocynthia roretzi mitochondrial genome encodes an
# additional tRNA gene with the anticodon U*CU that is thought to enable
# the use of AGA or AGG codons for glycine and the gene has been shown to
# be transcribed in vivo (Kondow et al.  1999, Yokobori et al.  1999). 


Genetic Code [13]

Ascidian Mitochondrial
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 14      Standard
# 
#  AAA      Asn  N       Lys  K
#  AGA      Ser  S       Arg  R
#  AGG      Ser  S       Arg  R
#  UAA      Tyr  Y       Ter  *
#  UGA      Trp  W       Ter  *
# 
# 
# Systematic Range: 
# 
# Platyhelminthes (flatworms) 
# 
# Comments:
#
# Code 14 differs from code 9 only by translating UAA to Tyr rather than
# STOP.  A recent study [PMID:11027335] has found no evidence that the
# codon UAA codes for Tyr in the flatworms but other opinions exist. 
# There are very few GenBank records that are translated with code 14 but
# a test translation shows that retranslating these records with code 9
# can cause premature terminations.  Therefore, GenBank will maintain code
# 14 until further information become available. 

Genetic Code [14]

Alternative Flatworm Mitochondrial
  
AAs  =   FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
# Version 3.2 - 6/24/95
#    Blepharisma Macronuclear code added
#
# Differences from the Standard Code: 
# 
#           Code 10       Standard
# 
# UAG       Gln  Q        Ter  *
# 
# 
# Systematic Range: 
# 
# Ciliata: Blepharisma (Liang and Heckman, 1993) 


Genetic Code [15]
  
Blepharisma Nuclear
  
AAs  =   FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
     
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.5
#     Added code 16, Chlorophycean Mitochondrial
#       (TAG can translated to Leucine instaed to STOP in chlorophyceans
#        and fungi)
#
# Systematic Range:
#
# Chlorophyceae: Hayashi-Ishiimaru, Y, T.  Ohama, Y.  Kawatsu, K. 
# Nakamura, S.  Osawa, 1996.  Current Genetics 30: 29-33

Genetic Code [16]

Chlorophycean Mitochondrial

AAs  =   FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.5
#     Added code 21, Trematode Mitochondrial
#       (as deduced from: Garey & Wolstenholme,1989; Ohama et al, 1990)
#
# Systematic Range:
#
# Trematoda: Ohama, T, S.  Osawa, K.  Watanabe, T.H.  Jukes, 1990.  J. 
# Molec Evol.  30 Garey, J.R.  and D.R.  Wolstenholme, 1989.  J.  Molec. 
# Evol.  28: 374-387 329-332. 
# 
# Other Alternative Initiation Codons
#
# GUG, UUG (and possibly CUG) in the Archaea (Noelling et al., unpublished) 
#
# AUA, GUG, UUG, and AUC or AAG may be used (at least in experimental
# systems) by the yeasts Saccharomyces cerevisiae (Olsen, 1987, and references
# therein). 
#
# ACG initiates translation of certain proteins in the adeno-associated virus
# type 2 (Becerra et al., 1985), the phage T7 mutant CR17 (Anderson and
# Buzash-Pollert, 1985), Sendai virus (Gupta and Patwardhan, 1988), and rice
# chloroplast (Hiratsuka et al., 1989). Also, it is the most effective non-AUG
# initiation codon in mammalin cells (Koepke and Leggatt, 1991). 
#
# CUG is the initiation codon for one of the two alternative products of the
# human c-myc gene (Hann et al., 1987). 
#

Genetic Code [21]

Trematode Mitochondrial

AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.6
#     Added code 22 TAG-Leu, TCA-stop
#        found in mitochondrial DNA of Scenedesmus obliquus
#        submitted by Gertraude Berger, Ph.D.
#        Organelle Genome Megasequencing Program, Univ Montreal
#
# Other Alternative Initiation Codons
#
# GUG, UUG (and possibly CUG) in the Archaea (Noelling et al., unpublished) 
# AUA, GUG, UUG, and AUC or AAG may be used (at least in experimental
# systems) by the yeasts Saccharomyces cerevisiae (Olsen, 1987, and references
# therein). 
#
# ACG initiates translation of certain proteins in the adeno-associated virus
# type 2 (Becerra et al., 1985), the phage T7 mutant CR17 (Anderson and
# Buzash-Pollert, 1985), Sendai virus (Gupta and Patwardhan, 1988), and rice
# chloroplast (Hiratsuka et al., 1989). Also, it is the most effective non-AUG
# initiation codon in mammalin cells (Koepke and Leggatt, 1991). 
#
# CUG is the initiation codon for one of the two alternative products of the
# human c-myc gene (Hann et al., 1987). 
#
# Systematic Range:
#
# Scenedesmus obliquus: Nedelcu A, Lee RW, Lemieux C, Gray MW and Burger G.
# "The complete mitochondrial DNA sequence of Scenedesmus obliquus reflects an
# intermediate stage in the evolution of the green algal mitochondrial genome."
# Genome Research (in press). 

Genetic Code [22]

Scenedesmus obliquus Mitochondrial

AAs  =   FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.7
#     Added code 23 Thraustochytrium mitochondrial code
#        formerly OGMP code 93
#        submitted by Gertraude Berger, Ph.D.
#
# This code has been created for the mitochondrial genome of the labyrinthulid
# Thraustochytrium aureum sequenced by the The Organelle Genome Megasequencing
# Program (OGMP).
#
# It is the similar to the bacterial code (trans_table 11) but it contains an
# additional stop codon (TTA) and also has a different set of start codons.

Genetic Code [23]

Thraustochytrium Mitochondrial

AAs  =   FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --------------------------------M--M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
