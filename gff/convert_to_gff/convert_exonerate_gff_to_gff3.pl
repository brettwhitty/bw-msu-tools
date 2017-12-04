#!/usr/bin/perl

##
## Quick script to convert exonerate GFF output into GFF3
##
## Requires that exonerate be run with the following flags:
##
## exonerate -Q [dna|protein] -T dna --model [est2genome|protein2genome] ... --verbose 0 --showalignment no --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no --ryo %ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n ...
##
## cigar and ryo output are used to populate GFF3 with additional attributes
##

use strict;
use warnings;
use Carp;

use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Getopt::Long;

use lib "/home/whitty/SVN/lib";
use MyIO;

my $input;
my $output;
my $version;
my $source;
my $no_skip;
my $match_type;
my $zero_strip;

my $result = GetOptions(
                           'input|i=s'      =>  \$input,
                           'output|o=s'     =>  \$output,
                           'source|s=s'     =>  \$source,
                           'no_skip!'       =>  \$no_skip,
                           'match_type|m=s' =>  \$match_type,
                           'zero_strip|z!'  =>  \$zero_strip,
                       );

my $infh = get_infh($input);
my $outfh = get_outfh($output);

my $in_gff = new Bio::Tools::GFF();

my $out_gff = new Bio::Tools::GFF(
    #                    -fh             =>  $outfh,
                        -gff_version    =>  3,
                    );

my $skip_feature = {
                        'utr5'      =>  1,
                        'utr3'      =>  1,
                        'utr3b'     =>  1,
                        'splice5'   =>  1,
                        'splice3'   =>  1,
                        'intron'    =>  1,
                        'cds'    =>  1,
                   };

## these are valid SO terms for match types
my $valid_match_type = {
    'match'                                 =>  1,
        'nucleotide_match'                  =>  1,
            'cross_genome_match'            =>  1,
            'expressed_sequence_match'      =>  1,
                'cDNA_match'                =>  1,
                'EST_match'                 =>  1,
                'RST_match'                 =>  1,
                'UST_match'                 =>  1,
            'primer_match'                  =>  1,
            'translated_nucleotide_match'   =>  1,
        'protein_match'                     =>  1,
};
## set the default match type
my $match_type_default = 'match';

## check that a valid match type has been provided
if (defined($match_type) && ! $valid_match_type->{$match_type}) {
    croak("Value provided to --match_type flag is not a valid SO match type, valid types are:\n"
          . join("\n", sort keys %{$valid_match_type})."\n");
}

$match_type = (defined($match_type)) ? $match_type : $match_type_default;

my $gff_flag;
my $query_seq_id;
my $gene_feature_id;
my ($target_acc, $query_acc, $qlen, $qlen_aln, $rank, $score, $pct_id, $pct_sim);
my ($target, $cigar);
my %feature_counter;
my @gff_features;
my $utr5_flag = 0;
while (<$infh>) {
    chomp;
    if (/^# seqname source/) {
        $gff_flag = 1;
        @gff_features = ();
        next;
    } elsif (/^# --- END OF/) {
        $gff_flag = 0;
        next;
    } elsif (/^#/) { next; }
    if ($gff_flag) {
        my $feature = new Bio::SeqFeature::Generic;
        $in_gff->from_gff_string($feature, $_);
#        use Data::Dumper;
#        print Dumper $feature;
#        print $out_gff->gff_string($feature)."\n";
        
        ## skip some feature types
        if (! $no_skip && $skip_feature->{$feature->primary_tag}) { 
            next;
        }
        if ($feature->has_tag('sequence')) {
            $query_seq_id = ($feature->get_tag_values('sequence'))[0];
        }
        
        if ($feature->has_tag('gene_id')) {
            $feature->remove_tag('gene_id');        
        }
        if ($feature->has_tag('alignment_id')) {
            $feature->remove_tag('alignment_id');        
        }
        if ($feature->has_tag('gene_orientation')) {
            $feature->remove_tag('gene_orientation');        
        }
        if ($feature->has_tag('Align')) {
            $feature->remove_tag('Align');        
        }
        if ($feature->has_tag('Query')) {
            $feature->remove_tag('Query');        
        }
        push(@gff_features, $feature);
    } elsif (scalar(my @t = split("\t")) == 8) {
        ($target_acc, $query_acc, $qlen, $qlen_aln, $rank, $score, $pct_id, $pct_sim) = split("\t");

        my $pct_cov = sprintf("%.1f", $qlen_aln / $qlen * 100);
        
        foreach my $feature(@gff_features) {
            
            if ($source) {
                $feature->source_tag($source);
            }
           

            if ($feature->primary_tag eq 'utr5') {
                $utr5_flag = 1;
                next;
            }
#            if ($utr5_flag) {
#                if ($feature->primary_tag eq 'exon') {
#                    $feature->primary_tag('five_prime_UTR');
#                } else {
#                    use Data::Dumper;
#                    print STDERR Dumper $feature;
#                    die "Unexpected output";   
#                }
#                $utr5_flag=0;
#            }
            
            if ($feature->primary_tag eq 'similarity') {
                ## change feature type
                $feature->primary_tag($match_type); 
            }
            if ($feature->primary_tag eq 'cds') {
                ## change feature type
                $feature->primary_tag('CDS'); 
            }

            if ($feature->primary_tag eq 'gene' || $feature->primary_tag eq $match_type) {
                $feature->add_tag_value('Name', $query_acc);
                my $feature_id = join("-", (
                        $feature->source_tag,
                        ($zero_strip) ?
                            ( map{s/(\D)[0]+(\d+)/$1$2/g; $_; } ($feature->seq_id) )[0]
                            : $feature->seq_id,
                        $feature->primary_tag,
                        ++$feature_counter{$feature->primary_tag}
                    ));
                $feature->add_tag_value('ID', $feature_id);
                if ($feature->primary_tag eq 'gene') {
                    $gene_feature_id = $feature_id;
                }
                $feature->add_tag_value('Rank', $rank);
                $feature->add_tag_value('identity', $pct_id);
                $feature->add_tag_value('similarity', $pct_sim);
                $feature->add_tag_value('coverage', $pct_cov);
            } else {
                $feature->add_tag_value('Parent', $gene_feature_id);
            }
        }

        ## print features
        foreach my $feature(@gff_features) {
            print $outfh $out_gff->gff_string($feature);
            if ($feature->primary_tag eq $match_type) {
#                $feature->add_tag_value('Target', $target);
#                $feature->add_tag_value('Gap', $cigar);
                print $outfh ";Target=$target;Gap=$cigar";
            }
            print $outfh "\n";
        }
 
    } elsif (/^cigar: (.*)/) {
        my @new_cigar;
        my $cigar_string = $1;
        my @cigar_arr = split(/\s+/, $cigar_string);
        $target = "$cigar_arr[0] $cigar_arr[1] $cigar_arr[2] $cigar_arr[3]";
        @cigar_arr = @cigar_arr[9 .. $#cigar_arr];
        while (my @pair = splice(@cigar_arr, 0, 2)) {
            push(@new_cigar, "$pair[0]$pair[1]");
        }
        $cigar = join(" ", @new_cigar);
   }
}
