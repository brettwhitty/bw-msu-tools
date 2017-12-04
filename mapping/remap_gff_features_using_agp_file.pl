#!/opt/rocks/bin/perl 

=head1  NAME 

remap_gff_features_using_agp_file.pl - Short description

=head1 SYNOPSIS

USAGE: remap_gff_features_using_agp_file.pl 
        --input=/path/to/input
        --output=/path/to/output

=head1 OPTIONS

B<--input,-i> 
    Description of input.

B<--output,-o> 
    Description of output.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

Description of script.

=head1 INPUT

Description of input.

=head1 OUTPUT

Description of output.

=head1 CONTACT

Brett Whitty
whitty@msu.edu

=cut

use strict;
use warnings;

use Carp;
use Pod::Usage;
use Getopt::Long;

use Bio::DB::Fasta;

use File::Temp qw{ tempfile };
use DB_File;

my ($help, $output_fasta, $child_fasta, $agp_file, $child_gff, $truncate);

GetOptions(
#    'parent_fasta|p=s'  =>  \$parent_fasta,
    'output_fasta|f=s'  =>  \$output_fasta,
    'child_fasta|c=s'   =>  \$child_fasta,
    'agp_file|a=s'      =>  \$agp_file,
    'child_gff|g=s'     =>  \$child_gff,
    'truncate|t!'      =>  \$truncate,
#    'input|i=s'     => \$input,
#    'output|o=s'    => \$output,
    'help|h!'       => \$help,
);

if ($help) {
    pod2usage(verbose => 2);
}

my (undef, $tempseqdb) = tempfile();
tie my %outseq, 'DB_File', $tempseqdb or confess "Failed to tie hash to temp file '$tempseqdb'";

#my $parent_db = new Bio::DB::Fasta($parent_fasta);
my $child_db = new Bio::DB::Fasta($child_fasta);

#my @parent_seq_ids = $parent_db->ids();
my @child_seq_ids  =  $child_db->ids();

my $c_lens = {};
foreach my $id(@child_seq_ids) {
    $c_lens->{$id} = $child_db->length($id);
}

my $infh;

open $infh, '<', $agp_file
    or confess "Failed to open AGP file '$agp_file' for reading: $!";

my $remap = {};
while (<$infh>) {
    chomp;

    my @t = split(/\t/, $_);

    my ($p_seq_id, $p_start, $p_end, undef, $i, $c_seq_id, $c_start, $c_end, $c_orient) = @t;
    
    my $subseq = '';
    if ($i eq 'W') {
        if (defined($output_fasta)) {
            $subseq = $child_db->seq($c_seq_id, $c_start, $c_end);
            if (defined($c_orient) && $c_orient eq '-') {
                $subseq = reverse($subseq);
            }
            $outseq{$p_seq_id} .= $subseq;    
        }
    } elsif ($i eq 'N') {
        if (defined($output_fasta)) {
            my $n_len = $p_end - $p_start + 1;
            $subseq = 'N' x $n_len;
            $outseq{$p_seq_id} .= $subseq;
        }
        next;
    } else {
        confess "AGP file '$agp_file' has unexpected content, unrecognized code '$i'";
    }
    
    if (! defined($c_orient)) {
        carp "No orientation information for '$c_seq_id' in AGP file '$agp_file', setting to '0'";
        $c_orient = '0';
    }


    push(@{$remap->{$c_seq_id}}, [ $c_start, $c_end, $c_orient, $p_seq_id, $p_start, $p_end ]);
}

if (defined($output_fasta)) {

## write the parent fasta file
open my $outfh, '>', $output_fasta
    or confess "Failed to open output fasta file '$output_fasta' for writing: $!";

foreach my $id(keys %outseq) {
    $outseq{$id} =~ s/(.{1,60})/$1\n/g;
    print $outfh ">$id\n".$outseq{$id};
}
}
##

## sort intervals on each child sequence by start position
foreach my $id(keys(%{$remap})) {
    @{$remap->{$id}} = sort {$a->[0] <=> $b->[0]} @{$remap->{$id}};
}
#use Data::Dumper;
#print Dumper $remap;
#print Dumper $c_lens;

open $infh, '<', $child_gff
    or confess "Failed to open '$child_gff' for reading: $!";

while (<$infh>) {
    chomp;

    my @gff = split("\t", $_);

    ## discard anything but feature lines
    if (scalar(@gff) != 9) { next; }

    ## we're only interested in features on molecules we need to remap
    if (! defined($remap->{$gff[0]})) {
        next;
    }

    ## for testing store ID
    $gff[8] =~ /ID=([^;]+)/;
    my $id = $1;
    ####

    ## strip trailing ';' in case the GFF is bad
    $gff[8] =~ s/;$//; 

    my @map = @{$remap->{$gff[0]}};
    my $handled = 0;

    foreach my $seg(@map) {
        my @remapped = undef;
        if ($gff[4] < $seg->[0]) {
            ## feature is completely outside of the left bound of the subseq included in the pseudomolecule
            #    carp "Feature '$id' on '$gff[0]' excluded from pseudomolecule, outside of left bound";
            $handled = 1;
        } elsif ($gff[3] < $seg->[0] && $gff[4] > $seg->[1]) {
            ## for supercontigs this will be the case if only a part is mapped to the chromosome
            if ($truncate) {
                @remapped = remap_coords( ($seg->[0], $gff[4], $gff[6], @{$seg}) );
                if ($gff[8] !~ /partial_left/ && $gff[8] !~ /partial_right/) {
                    $gff[8] .= ';partial_left=1;partial_right=1';
                }
            } else {
                @remapped = remap_coords( ($gff[3], $gff[4], $gff[6], @{$seg}) );
            }
        } elsif ($gff[3] < $seg->[0] && $gff[4] >= $seg->[0] && $gff[4] <= $seg->[1]) {
            ## feature starts out of the left bound of the subseq included in the pseudomolecule
            #carp "Feature '$id' extends out of left bound";
            ## I'm going to remap anyway
            if ($truncate) {
                ## truncate the feature
                @remapped = remap_coords( ($seg->[0], $gff[4], $gff[6], @{$seg}) );
                if ($gff[8] !~ /partial_left/) {
                    $gff[8] .= ';partial_left=1';
                }
            } else {
                ## remap ignoring that the coordinates fall outside of the child bounds
                @remapped = remap_coords( ($gff[3], $gff[4], $gff[6], @{$seg}) );
            }
            $handled = 1;
        } elsif ($gff[3] >= $seg->[0] && $gff[4] <= $seg->[1]) {
            ## feature is in the segment, we need to remap
            #carp "Feature '$id' should be remapped";
            @remapped = remap_coords( ($gff[3], $gff[4], $gff[6], @{$seg}) );
            $handled = 1;
        } elsif ($gff[3] >= $seg->[0] && $gff[3] <= $seg->[1] && $gff[4] > $seg->[1]) {
            ## feature extends out of the right bounds of the subseq included in the pseudomolecule
            #carp "Feature '$id' extends out of right bound: ".$seg->[0]."-".$seg->[1]." ".$gff[3]."-".$gff[4];
            ## I'm going to remap anyway
            if ($truncate) {
                ## truncate the feature
                @remapped = remap_coords( ($gff[3], $seg->[1], $gff[6], @{$seg}) );
                if ($gff[8] !~ /partial_right/) {
                    $gff[8] .= ';partial_right=1';
                }
            } else {
                ## remap ignoring that the coordinates fall outside of the child bounds
                @remapped = remap_coords( ($gff[3], $gff[4], $gff[6], @{$seg}) );
            }
            $handled = 1;
        }
        if (scalar(@remapped) == 4) {
            $gff[0] = $remapped[0];
            $gff[3] = $remapped[1];
            $gff[4] = $remapped[2];
            $gff[6] = $remapped[3];
            print join("\t", @gff)."\n";
        }
    }
    if (! $handled) {
            ## feature is presumably outside of the right bounds of the subseq in the pseudomolecule
            #carp "Feature '$id' on '$gff[0]' excluded from pseudomolecule, outside of right bound";
    }
}

## remap the feature coordinate based on the child to parent mappings from the AGP file
sub remap_coords {
    my ($f_start, $f_end, $f_orient, $c_start, $c_end, $c_orient, $p_seq_id, $p_start, $p_end) = @_;

    my @remapped = ();

    if ($c_orient eq '-') {
        ## should be correct for reverse orientation
        my $p_f_start = $p_end - ($f_start - $c_start);
        my $p_f_end   = $p_end - ($f_end   - $c_start);

        my $p_f_orient;
       
        if ($f_orient eq '-') {
            $p_f_orient = '+';
        }
        elsif ($f_orient eq '+') {
            $p_f_orient = '-';
        } else {
            $p_f_orient = $f_orient;
        }

        ## in GFF3 start < end so 
        @remapped = ($p_seq_id, $p_f_end, $p_f_start, $p_f_orient);

    } else {
        ## should be correct for forward orientation
        my $p_f_start = $f_start - $c_start + $p_start;
        my $p_f_end   = $f_end   - $c_start + $p_start;

        my $p_f_orient = $f_orient;

        @remapped = ($p_seq_id, $p_f_start, $p_f_end, $p_f_orient);
    }    

    if ($remapped[1] < 0) {
        carp "WARNING: remapped start coordinate is <1: $remapped[0]";
    }
    if ($remapped[2] < 0) {
        carp "WARNING: remapped end coordinate is <1: $remapped[1]";
    }

    return @remapped;
}

