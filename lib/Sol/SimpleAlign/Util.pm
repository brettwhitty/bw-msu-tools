#!/usr/bin/perl

package Sol::SimpleAlign::Util;

## Can be used to generate HTML representing a Bio::SimpleAlign multiple sequence alignment
## 
##
## my @coords = $db->get_mismatch_coords($aln_obj); 
##
## print $db->get_alignment_html(
##          $aln, 
##          "each_character"    => 1, 
##          "char_class"        => {
##                  "mismatch" => \@coords
##                                 },  
##       );
##

use strict;
use warnings;

use Carp;

sub new {
    my ($class) = @_;

    return bless {}, ref($class) || $class;
}

## returns a list containing the coordinates of mismatched bases
## for use with coloring the alignment
## can be passed an ID or an alignment object
sub get_mismatch_coords {
    my ($self, $aln_obj) = @_;

#    my $aln_obj = (ref($ref) eq 'Bio::SimpleAlign') ? $ref : $self->get_alignment($ref);

    my $gap_line    = $aln_obj->gap_line();
    my $match_line  = $aln_obj->match_line();

    my $aln_len = $aln_obj->length();
   
    ## mismatches at the ends are due to gaps (presumably) 
    my $start_gap_pos = -1;
    if ($match_line =~ /(^\s+)/) {
        $start_gap_pos = $+[1];
    }
    ## mismatches at the ends are due to gaps (presumably) 
    my $end_gap_pos = $aln_len;
    if ($match_line =~ /(\s+$)/) {
        $end_gap_pos = $-[1];
    }

    my @matches = split(//, $match_line);

    my @mismatch_coords = ();
    for (my $i = 0; $i < $aln_len; $i++) {
        if ($matches[$i] ne '*' && $i > $start_gap_pos && $i < $end_gap_pos) {
            push(@mismatch_coords, $i + 1);
        }
    }

    return @mismatch_coords;
}

sub get_gap_coords {
    my ($self, $id) = @_;

    croak "Not implemented yet";
}

## returns a html for an alignment
sub get_alignment_html {
    my ($self, $aln_obj, %args) = @_;

    # my $aln_obj = (ref($ref) eq 'Bio::SimpleAlign') ? $ref : $self->get_alignment($ref);

    ## class of MSA container div
    my $container_class = (defined($args{'msa'}))          
        ? $args{'msa'}           : 'msa';
    ## id of MSA container div
    my $container_id    = (defined($args{'msa_container'}))
        ? $args{'msa_container'} : 'msa_container';
    ## class of MSA legend
    my $legend_class    = (defined($args{'legend_class'}))
        ? $args{'legend_class'}         : 'msa_legend';
    ## class of MSA legend
    my $legend_class_seq    = (defined($args{'legend_class_seq'}))
        ? $args{'legend_class_seq'}     : 'msa_legend_seq';
    ## class of MSA legend
    my $legend_class_con    = (defined($args{'legend_class_con'}))
        ? $args{'legend_class_con'}     : 'msa_legend_con';
    ## class of MSA member sequences
    my $seqs_class      = (defined($args{'seqs_class'}))
        ? $args{'seqs_class'}    : 'msa_seqs';
    ## class of MSA member sequence
    my $seq_class       = (defined($args{'seq_class'}))
        ? $args{'seq_class'}     : 'msa_seq';
    ## class of MSA consensus sequence
    my $con_class       = (defined($args{'con_class'}))
        ? $args{'con_class'}     : 'msa_con_seq';
    ## ID of consensus sequence
    my $con_id          = (defined($args{'con_id'}))
        ? $args{'con_id'}        : 'msa_con_seq';
    ## spacer class
    my $spacer_class    = (defined($args{'spacer_class'}))
        ? $args{'spacer_class'}     : 'msa_spacer';

    ## add class to each character of sequence
    my $each_character  = (defined($args{'each_character'})) 
        ? $args{'each_character'}   :   0;
    ## allow user to define character class lookup table for use with each_character option
    my $char_class_lookup = (defined($args{'char_class_lookup'}))
        ? $args{'char_class_lookup'}    :   undef;

    ## user can provide a hash of sequence labels to use in the MSA legend
    ## eg: [ 'consensus', 'something1', 'something2', ...]
    my $legend_labels = (defined($args{'legend_labels'}))
        ? $args{'legend_labels'}    :   undef;


    ## take user-specified class coordinates and rearrange to coordinate classes
    ## FROM:
    ## 'char_class' => { 'class' => [ pos1, pos2, ... ], ... }
    ##
    my $coord_class = {};
    if (defined($args{'char_class'})) {
        my $char_class_ref = $args{'char_class'};
        my @char_classes = keys %{$char_class_ref};
        foreach my $char_class(@char_classes) {
            foreach my $pos(@{$char_class_ref->{$char_class}}) {
                push(@{$coord_class->{$pos}}, $char_class);
            }
        }
    }

    my $seq_count = $aln_obj->no_sequences();

    my $html = '';

    ## msa div
    $html .= _html_line("<div class='$container_class' id='$container_id'>", 0, 1);
    ## legend div
    $html .= _html_line("<div class='$legend_class' id='$legend_class'>", 1, 1);

    ## consensus is the first sequence, so prepare an iterator array
    my @iter = ( 2 .. $seq_count, 1);
    foreach my $i(@iter) {
        my $text;
        if (defined($legend_labels)) {
            $text = $legend_labels->[$i - 1];
        } else {
            ## if consensus, seq_name is expected to be accession for alignment
            ## otherwise it's the display_name of the member sequence
            my $seq_name = ($i == 1) ? $aln_obj->accession() : $aln_obj->get_seq_by_pos($i)->display_name;
            ## we'll append the description for the sequence if there is one
            my $seq_desc = $aln_obj->get_seq_by_pos($i)->desc;

            if ($seq_desc) {
                $text = "$seq_name $seq_desc";
            } else {
                $text = $seq_name;
            }
        }

        ## sequence is not the consensus
        if ($i != 1) {
            my $id = $legend_class_seq . '_' . ($i - 1);
            $html .= _html_line("<div class='$legend_class_seq' id='$id'>", 2, 0);
        ## sequence is the consensus
        } else {
            my $id = $legend_class_con;
            ## add a spacer line
            $html .= _html_line("<div class='$spacer_class'></div>", 2, 1);
            $html .= _html_line("<div class='$legend_class_con' id='$id'>", 2, 0);
        }
        $html .= _html_line($text.'</div>', 0, 1);
    }
    ## close legend div
    $html .= _html_line("</div>", 1, 1);
    ## seqs div
    $html .= _html_line("<div class='$seqs_class' id='$seqs_class'>", 1, 1);
    ## member sequences
    for (my $i = 2; $i <= $seq_count; $i++) {
        my $count = $i - 1;
        my $id = $seq_class.'_'.$count;

        $html .= _html_line("<div class='$seq_class' id='$id'>", 2, 0);
        $html .= _seq_to_div(
                    $aln_obj->get_seq_by_pos($i)->seq,
                        'each_character'    =>  $each_character,
                        'char_class_lookup' =>  $char_class_lookup,
                        'coord_class'       =>  $coord_class,
                        'id_prefix'         =>  $id,
                            );
        $html .= _html_line("</div>", 0, 1);
    }
    ## consensus sequence
    #### add a spacer line
    $html .= _html_line("<div class='$spacer_class'></div>", 2, 1);
    $html .= _html_line("<div class='$con_class' id='$con_id'>", 2, 0);
    $html .= _seq_to_div(
                    $aln_obj->get_seq_by_pos(1)->seq,
                        'each_character'    =>  $each_character,
                        'char_class_lookup' =>  $char_class_lookup,
                        'coord_class'       =>  $coord_class,
                        'id_prefix'         =>  $con_id,
                );
    $html .= _html_line("</div>", 0, 1);
    ## close seqs div
    $html .= _html_line("</div>", 1, 1);
    ## close msa div
    $html .= _html_line("</div>", 0, 1);

}

## convert a sequence string to a div
sub _seq_to_div {
    my ($seq, %args) = @_;

    ## character class prefix
    my $char_class_prefix  = (defined($args{'char_class_prefix'})) 
        ? $args{'char_class_prefix'}    :   'msa_char';
    ## classes to be added at specific character positions
    my $coord_class = (defined($args{'coord_class'}))
        ? $args{'coord_class'}          :   undef;
    ## prefix for adding ids to the character divs
    my $id_prefix   = (defined($args{'id_prefix'}))
        ? $args{'id_prefix'}            :   undef;
    ## optionally add class to each character
    my $each_character = (defined($args{'each_character'})) 
        ? $args{'each_character'}       :   0;
    ## allow character class lookup as argument
    my $char_class_lookup = (defined($args{'char_class_lookup'}))
        ? $args{'char_class_lookup'}    :                  {
                'A'   =>    'a',
                'a'   =>    'a',
                'C'   =>    'c',
                'c'   =>    'c',
                'G'   =>    'g',
                'g'   =>    'g',
                'T'   =>    't',
                't'   =>    't',
                '-'   =>    'gap',
                'N'   =>    'n',
                'n'   =>    'n',
                'X'   =>    'x',
                'x'   =>    'x',
                                                           };
 
    my $html = '';

    my $seq_len = length($seq);

    for (my $i = 0; $i < $seq_len; $i++) {
        
        ## get character at the position
        my $char = substr($seq, $i, 1);

        ## for 1-based coords on sequence
        my $pos = $i + 1;
        
        ## for storing the divs
        my $left_html = '';
        my $right_html = '';
        
        ## optionally add class by character type
        if ($each_character) {
            if (defined($char_class_lookup->{$char})) {
                my $class = $char_class_prefix.'_'.$char_class_lookup->{$char};
                $left_html  .= "<span class='$class'>";
                $right_html .= "</span>";
            }
        }

        ## optionally add classes by position
        if (defined($coord_class->{$pos})) {
            foreach my $class(@{$coord_class->{$pos}}) {
                
                my $class_type = $char_class_prefix.'_'.$class;

                my $id_string = '';
                if (defined($id_prefix)) {
                    $id_string = " id='".$id_prefix.'_'.$class.'_'.$pos."'";
                }

                $left_html .= "<span class='$class_type'$id_string>";
                $right_html .= "</span>";
            }
        }
        $html .= $left_html . $char . $right_html;
    }

    return $html;
}

## used for formatting the HTML lines
sub _html_line {
    my ($text, $level, $newline) = @_;

    my $pad_string = '    ';

    $level = int($level);

    unless ($level) {
        $level = 0;
    }

    my $pad = $pad_string x $level;

    my $string = "$pad$text";
    if ($newline) {
        $string .= "\n";
    }

    return $string;
}


## this is to fix Bio::SimpleAlign's version for alignments
## with end gaps in the alignment
sub sort_by_start {
    my ($self, $aln_obj, $has_reference) = @_;

    ## 'has_reference' flag is used to indicate whether
    ## the MSA contains a reference sequence (0 = no, 1 = yes)
    ## in which case the reference will be assumed to be currently
    ## at position 0 and will be left out of the sort
    if (! defined($has_reference)) {
        $has_reference = 1;
    }

    my %hash;
    my $ref_flag = 0;
    foreach my $seq ( $aln_obj->each_seq() ) {
        ## if there's a reference sequence, it's going to retain
        ## it's position at 0
        if (! $ref_flag && $has_reference) {
            $ref_flag = 1;
            next;
        }
        
        my $nse = $seq->get_nse;

        ## my replacement for LocatableSeq::get_nse
        ## that takes end gaps into account
        my $seq_name = $seq->display_name;
        $seq->seq =~ /^[-]*/;
        my $left_pos = $+[0] + 1;
        $seq->seq =~ /[-]*$/;
        my $right_pos = $-[0] - 1;
        my $adj_nse = "$seq_name/$left_pos-$right_pos";

        $hash{$adj_nse} = $nse;
    }


    my $count = 0;
    ## start at 1 if there's a reference sequence
    if ($has_reference) {
        %{$aln_obj->{'_order'}} = (0 => $aln_obj->{'_order'}->{'0'}); # reset the hash;
        $count = 1;
    } else {
        %{$aln_obj->{'_order'}} = (); # reset the hash;
    }
    foreach my $nse_adj ( sort _startend keys %hash ) {
        $aln_obj->{'_order'}->{$count} = $hash{$nse_adj};
        $count++;
    }

    return 1;
}

## returns an array of arrays containing member accessions, start, end, and
## Bio::SeqFeature::Generic format strand positions in the alignment in
## whatever the current order is in the alignment object
sub get_member_coords {
    my ($self, $aln_obj, $has_reference) = @_;

    ## 'has_reference' flag is used to indicate whether
    ## the MSA contains a reference sequence (0 = no, 1 = yes)
    ## in which case the reference will be assumed to be currently
    ## at position 0 and will be left out of the sort
    if (! defined($has_reference)) {
        $has_reference = 1;
    }

    my @members;
    my $ref_flag = 0;
    foreach my $seq ( $aln_obj->each_seq() ) {
        ## if there's a reference sequence, it's going to retain
        ## it's position at 0
        if (! $ref_flag && $has_reference) {
            $ref_flag = 1;
            next;
        }

        my $nse = $seq->get_nse;

        ## my replacement for LocatableSeq::get_nse
        ## that takes end gaps into account
        my $seq_name = $seq->display_name;
        $seq->seq =~ /^[-]*/;
        my $left_pos = $+[0] + 1;
        $seq->seq =~ /[-]*$/;
        my $right_pos = $-[0] - 1;

        my $orient = $seq->desc; ## + or -
        
        push(@members, [
            $seq_name,
            $left_pos,
            $right_pos,
            ($orient eq '-') ? -1 : 1,
        ]);
    }

    return @members;
}

sub _startend {
    my ($aname, $arange) = split (/[\/]/, $a);
    my ($bname, $brange) = split (/[\/]/, $b);
    my ($astart, $aend)  = split(/\-/, $arange);
    my ($bstart, $bend)  = split(/\-/, $brange);

    return $astart <=> $bstart;
}

1;
