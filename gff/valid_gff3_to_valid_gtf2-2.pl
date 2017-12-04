#!/opt/rocks/bin/perl 

=head1  NAME 

valid_gff3_to_valid_gtf2-2.pl - Tries to convert spec-compliant GFF3 to spec-compliant GTF 2.2

=head1 SYNOPSIS

USAGE: valid_gff3_to_valid_gtf2-2.pl 
        --input=/path/to/gene_describing.gff3
        --output=/path/to/gene_describing.gtf
       [--keep-features] [--exons] [--utrs] 
       [--infer-start-stop] [--ignore_frame]
       [--keep-attributes] [--keep-features] [--keep-comments]

=head1 OPTIONS

B<--input,-i> 
    A valid GFF3 file containing holistic GFF3-encoded gene features.

B<--output,-o> 
    A hopefully valid GTF 2.2 file.

B<--exons,-e>
    Output exons features if provided in input (optional valid GTF 2.2 features)

B<--utrs,-u>
    Output UTR features if provided in input (optional valid GTF 2.2 features)

B<--infer-start-stop,-m>
    If a gene doesn't have a start codon defined, presume it is the first 3 bases
    at the start of the first CDS. If a gene doesn't have a stop codon defined,
    presume it is the 3 bases after the end of the last CDS (or first 3 bases
    of upstream 3' UTR region if present).

B<--ignore_frame,-f>
    If frame information is missing for features where it is required, suck it up.

B<--keep-comments,-c>
    Keep comments from the top of the GFF3 file.

B<--keep-features,-k>
    Keep features that are non-GTF-2.2-spec.

B<--keep-attributes,-a>
    Keep attributes that are not required in GTF.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

Does some dirty business to convert GFF3 into GTF for supplying as input to 
any tools that inflexibly require GTF input, such as TopHat and Cufflinks.

*** WARNING ***
Not all optional features are implemented yet. Automatic addition of missing
features such as start and stop codons is not yet functional. If GFF3 contains
valid GTF 2.2 mapped features then they should be converted without problems.

=head1 INPUT

A spec-compliant GFF3 file that encodes genes (eg: a genome annotation file, or genefinder output)

=head1 OUTPUT

A hopefully spec-compliant GTF file representing the input gene annotation.

=head1 CONTACT

Brett Whitty
whitty@msu.edu

=cut

use strict;
use warnings;

use Carp;
use Pod::Usage;
use Getopt::Long;

my ($input, $output, $help, $exons, $introns, $utrs, $start_stop, $ignore_frame, $keep_atts, $keep_feats, $keep_comments);

GetOptions(
    'input|i=s'             =>  \$input,
    'output|o=s'            =>  \$output,
    'exons|e!'              =>  \$exons,
    'introns|n!'            =>  \$introns,
    'utrs|u!'               =>  \$utrs,
    'infer-start-stop|m!'   =>  \$start_stop,
    'ignore-frame|f!'       =>  \$ignore_frame,
    'keep-attributes|a!'    =>  \$keep_atts,
    'keep-features|k!'      =>  \$keep_feats,
    'keep-comments|c!'      =>  \$keep_comments,
    'help|h!'               =>  \$help,
);

if ($help) {
    pod2usage(verbose => 2);
}

## hash to store arrays of feats keyed by sequence id
my $gtf_feats = {};
my $other_feats = {};
my $gff_feat_types = {};
my $gff_feat_parents = {};

my $feat_map = {
    'start_codon'       =>  'start_codon',
    'five_prime_UTR'    =>  '5UTR',
    'CDS'               =>  'CDS',
    'intron'            =>  'intron_CNS', 
    'exon'              =>  'exon',
    'three_prime_UTR'   =>  '3UTR',
    'stop_codon'        =>  'stop_codon',
#    'xxx'               =>  'inter',           ## valid GTF features
#    'xxx'               =>  'inter_CNS',       ## valid GTF features
};

## filehandle for reading from input file
my $infh;

open $infh, '<', $input
    or croak "Failed to open input GFF3 file '$input' for reading: $!";

## parse the file once to get the feature relationships
while (<$infh>) {
    chomp;
    
    my @cols = split(/\t/, $_);

    if (scalar(@cols) != 9) {
        next;
    }
   
    my @atts = split(/;/, $cols[8]);
    my $id = '';
    my $parents = [];
    foreach my $att(@atts) {
        my ($key, $vals) = split(/=/, $att);
        my @vals = split(/,/, $vals);
        if ($key eq 'ID') {
            $id = $vals[0];
        } elsif ($key eq 'Parent') {
            $parents = \@vals;
        }
    }

    if ($id ne '') {
        ## ignore features without IDs
        ## store parents for feature
        foreach my $parent_id(@{$parents}) {
            push(@{$gff_feat_parents->{$id}}, $parent_id);
        }
        ## store GFF3 feature type
        $gff_feat_types->{$id} = $cols[2];
    }
}
#use Data::Dumper;
#print Dumper $gff_feat_parents;
#die();
open $infh, '<', $input
    or croak "Failed to open input GFF3 file '$input' for reading: $!";

my @comments;
my $flag = 0;
while (<$infh>) {
    chomp;
    
    my @cols = split(/\t/, $_);

    if (scalar(@cols) == 9) {
        ## row describing feature
        $flag = 1;
        process_feature(@cols);
    } elsif ($keep_comments && ! $flag && /^#/) {
        ## comment row at top of file
        push(@comments, $_);
    } else {
        ## some other content, eg: /^##$/ or fasta
        next;
    }

}

sub process_feature {
    my ($seq_id, $source, $feat_type, $start, $end, $score, $strand, $frame, $atts) = @_;

    if (defined($feat_map->{$feat_type})) {
        ## recognized valid GTF feature
        my $gtf_type = $feat_map->{$feat_type};

        if ($gtf_type eq 'exon' && ! $exons) {
            return 0;
        } elsif ($gtf_type =~ /[53]UTR/ && ! $utrs) {
            return 0;
        } elsif ($gtf_type eq 'intron_CNS' && ! $introns) {
            return 0;
        }

        my @gtf_att_other = ();
        my @atts = split(/;/, $atts);
        my $id = '';
        foreach my $att(@atts) {
            my ($key, $vals) = split(/=/, $att);
            my @vals = split(/,/, $vals);
            if ($key eq 'ID') {
                $id = $vals[0];
            } elsif ($keep_atts && $key !~ /^Parent$/) {
                ## support keeping atts outside the GTF spec
                ## not sure of the correct way to provide more than
                ## one val per att in GTF 2.2, so keeping it like GFF3
                push(@gtf_att_other, "$key \"".join(',', @vals)."\";");
            }
        }

        ## in GFF3 features could have more than one parent, so these return an array
        my @gtf_transcript_ids = get_parent_ids('mRNA', $id);

        foreach my $gtf_transcript_id(@gtf_transcript_ids) {
            my ($gtf_gene_id) = get_parent_ids('gene', $gtf_transcript_id);

            print join("\t", (
            #push(@{$gtf_feats->{$seq_id}}, [
                    $seq_id,    ## added
                $source,
                $gtf_type,
                $start,
                $end,
                $score,
                $strand,
                $frame,
                #$gtf_atts,
                "gene_id \"$gtf_gene_id\"; transcript_id \"$gtf_transcript_id\";",
#            ]); 
            ))."\n";
        }

    } elsif ($keep_feats) {
            print join("\t", (
#        push(@{$gtf_feats->{$seq_id}}, [
                    $seq_id,    ## added
            $source,
            $feat_type,
            $start,
            $end,
            $score,
            $strand,
            $frame,
            "gene_id \"\"; transcript_id \"\";",
            #$gtf_atts,
#        ]);  
        ))."\n";
    }
}

#use Data::Dumper;
#print Dumper $gtf_feats;

sub get_parent_ids {
    my ($parent_type, $id) = @_;

    my $transcript_ids = {};

    if (! defined($gff_feat_types->{$id})) {
        ## this shouldn't happen, but just in case
        carp "WARNING: No feat_type defined in GFF for ID '$id'";
        return ('');
#    } elsif (! defined($feat_map->{$gff_feat_types->{$id}})) {
#        ## don't return transcript IDs for non-GTF feature types
#        carp "WARNING: No GTF feature type for feature '$id'";
#        return ('');
    } else {
        my @parent_ids = (defined($gff_feat_parents->{$id})) ? @{$gff_feat_parents->{$id}} : (); 
        while (scalar(@parent_ids) > 0) {
            my @grandparents = ();
            foreach my $parent_id(@parent_ids) {
                my $parent_feat_type = $gff_feat_types->{$parent_id} || '';
                if ($parent_feat_type eq $parent_type) {
                    $transcript_ids->{$parent_id} = 1;
                } else {
                    if (defined($gff_feat_parents->{$parent_id})) {
                        push(@grandparents, @{$gff_feat_parents->{$parent_id}});
                    }
                }
            }
            @parent_ids = @grandparents;
        }
        return sort keys(%{$transcript_ids});
    }
}
