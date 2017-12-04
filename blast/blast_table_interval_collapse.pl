#!/usr/local/bin/perl

##
## Takes NCBI BLAST tabular output (-m 8)
## and generates collapsed intervals 
## with the specified number of tiers 
## in an output format similar to
## filter.pl from the AAT component
##
## Brett Whitty
## whitty@msu.edu
##

use warnings;
use strict;

my $qlen = 0;
my $slen = 0;

my $result_hash = {};

my @blast_results = ();

while (<>) {
    chomp;
    my @t = split("\t");
    ## query_id  subject_id  percent_id  alignment_len  mismatches  gap_openings  query_start  query_end  subject_start  subject_end  e_value  bit_score
    ##  0           1           2           3           4               5           6           7           8               9           10      11
    
    my ($q_comp, $s_comp) = (0, 0);
    if ($t[6] > $t[7]) {
        ($t[7], $t[6]) = ($t[6], $t[7]);
        $q_comp = 1;
    }
    if ($t[8] > $t[9]) {
        ($t[9], $t[8]) = ($t[8], $t[9]);
        $s_comp = 1;
    }
    my $orientation;
    if (($q_comp && $s_comp) || (! $q_comp && ! $s_comp)) {
        $orientation = 0;
    } else {
        $orientation = 1;
    }
    push(@{$result_hash->{$t[0]}}, pad(8, $t[6])." ".pad(8, $t[7])." ".pad(6,$t[11])." ".pad(7, $t[8])." ".pad(5,$t[9])." ".$orientation." ".pad(5,'0')." ".pad(5,'0')." ".$t[1]."\n");
}

sub pad {
    my ($width, $string) = @_;
    return " " x ($width - length($string)) . $string;
}

## extCollapse.pl

my $collapsed_hash = {};

foreach my $query(keys(%{$result_hash})) {
    my @chains;
foreach (@{$result_hash->{$query}}) {
    unless (/\w/) { next;}
    chomp;
    $_ =~ s/^\s+//; #rm leading whitespace
    ## using var names as in ext.c
    my ($dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc) = split (/\s+/);
    
    push (@chains, 
      { dstart=> $dstart,
        dend => $dend,
        score => $score,
        astart => $astart,
        aend => $aend,
        orient => $orient,
        zero1 => $zero1,
        zero2 => $zero2,
        acc => $acc } );
    
}

## output collapsed chains:
unless (@chains) {
    # no output to process
    exit(0);
}

@chains = sort { 
    $a->{acc} cmp $b->{acc}
    ||
        $a->{orient} <=> $b->{orient} 
    ||
        $a->{dstart} <=> $b->{dstart}
} @chains;

my @collapsedChains  = shift @chains;

while (@chains) {
    my $prevChain = $collapsedChains[$#collapsedChains];
    
    my ($prev_dstart, 
        $prev_dend, 
        $prev_score, 
        $prev_astart, 
        $prev_aend, 
        $prev_orient, 
        $prev_acc) = ($prevChain->{dstart},
                      $prevChain->{dend},
                      $prevChain->{score},
                      $prevChain->{astart},
                      $prevChain->{aend},
                      $prevChain->{orient},
                      $prevChain->{acc});
    
    my $nextChain = shift @chains;
    my ($next_dstart, 
        $next_dend,
        $next_score,
        $next_astart,
        $next_aend,
        $next_orient,
        $next_acc) = ($nextChain->{dstart},
                      $nextChain->{dend},
                      $nextChain->{score},
                      $nextChain->{astart},
                      $nextChain->{aend},
                      $nextChain->{orient},
                      $nextChain->{acc});
    
    if ( $next_acc eq $prev_acc 
         && 
         $next_orient == $prev_orient
         && 
         $next_dstart <= $prev_dend) {
        
        ## merge overlapping entry:
        my @dcoords = sort {$a<=>$b} ($prev_dstart, $prev_dend, $next_dstart, $next_dend);
        my $dstart = shift @dcoords;
        my $dend = pop @dcoords;
        
        my @acoords = sort {$a<=>$b} ($prev_astart, $prev_aend, $next_astart, $next_aend);
        my $astart = shift @acoords;
        my $aend = pop @acoords;
        
        my @scores = sort {$a<=>$b} ($prev_score, $next_score);
        my $score = pop @scores;
        
        ## adjust prevChain contents
        $prevChain->{dstart} = $dstart;
        $prevChain->{dend} = $dend;
        $prevChain->{astart} = $astart;
        $prevChain->{aend} = $aend;
        $prevChain->{score} = $score;
        
    } else {
        
        ## terminate previous chain
        push (@collapsedChains, $nextChain);
        
    }
}

## Sort collapsed chains by dstart and score:

@collapsedChains = sort { $a->{dstart} <=> $b->{dstart}
|| 
    $b->{score} <=> $a->{score} } @collapsedChains;

my @collapsed_result = ();

foreach my $chain (@collapsedChains) {
    my ($dstart, 
        $dend,
        $score,
        $astart,
        $aend,
        $orient,
        $acc, 
        $zero1,
        $zero2 ) = ($chain->{dstart},
                    $chain->{dend},
                    $chain->{score},
                    $chain->{astart},
                    $chain->{aend},
                    $chain->{orient},
                    $chain->{acc},
                    $chain->{zero1},
                    $chain->{zero2}
                    );

push(@{$collapsed_hash->{$query}}, sprintf("%8d %8d %6d %7d %5d %1d %5d %5d %s\n",
       $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc));

}

}

## filter.pl

main: {
    
    my $tier_num = 1;
   

foreach my $query(keys(%{$collapsed_hash})) {
 
    
    chomp @{$collapsed_hash->{$query}};
    
    my @features;
    my $count = 0;
    foreach my $chain (@{$collapsed_hash->{$query}}) {
        $count++;
        my $text = $chain;
        $text =~ s/^\s+//;
        my @x = split (/\s+/, $text);
        my ($genome_lend, $genome_rend, $chain_score) = ($x[0], $x[1], $x[2]);
        
        my $feature = new TierFeatures::Feature($genome_lend, $genome_rend, $count);
        
        $feature->{chain_text} = $chain; # extract later
        $feature->{chain_score} = $chain_score;

        push (@features, $feature);
    }

    my $feature_tierer = new TierFeatures($tier_num);
    
    ## sort features by score, desc, so best features are tiered first:
    @features = reverse sort {$a->{chain_score} <=> $b->{chain_score}} @features; 

    my @tiers = $feature_tierer->tier_features(@features);

    my @tiered_features;
    for (my $i = 1; $i <= $tier_num; $i++) {
        my $tier = shift @tiers;
        if (ref $tier) {
            push (@tiered_features, @$tier);
        } 
        else {
            # no tiers left.  all done.
            last;
        }
    }

    ## put back in same order from which they were derived:
    @tiered_features = sort {$a->{feat_id}<=>$b->{feat_id}} @tiered_features;

    ## print report:
    ##print $header;
    foreach my $feature (@tiered_features) {
        print "$query " . $feature->{chain_text} . "\n";
    }
    }
    exit(0);
}

package TierFeatures;
use strict;
use warnings;

sub new {
    my $packagename = shift;

    my $max_tier = shift;
    
    
    my $self = { tiers => [ [] ], # give it one empty tier to start
                 max_tier => $max_tier,  # if undef, then no maximum.
                 
                 max_percent_overlap_allowed => 20,  #if overlap by this percent more than either feature length, cannot coexist on same tier.

};

    bless ($self, $packagename);

    return ($self);
}


sub set_max_percent_overlap {
    my $self = shift;
    
    my $max_percent = shift;
    if ($max_percent < 0 || $max_percent > 100) {
        die "Error, percent ($max_percent) is not allowed.\n";
    }

    $self->{max_percent_overlap_allowed} = $max_percent;
}





sub tier_features {
    my $self = shift;
    my @features = @_;
            
    my $max_tier = $self->{max_tier};
    my $max_percent_overlap_allowed = $self->{max_percent_overlap_allowed};


    ## start at first tier:
    foreach my $feature (@features) {
        my ($feat_lend, $feat_rend) = ($feature->{lend}, $feature->{rend});
        my @tiers = @{$self->{tiers}};
        
        my $feature_length = $feat_rend - $feat_lend + 1;

        my $tiered_feature_flag = 0;
        
      tiers:
        foreach my $tier (@tiers) {
            my @tiered_feats = @$tier;
          feats:
            foreach my $feat (@tiered_feats) {
                my ($lend, $rend) = ($feat->{lend}, $feat->{rend});
                # check for overlap
                my $feat_length = $rend - $lend + 1;

                if ($lend <= $feat_rend && $rend >= $feat_lend) {
                    # got overlap
                    
                    my $nucs_common = $self->_nucs_in_common($feat_lend, $feat_rend, $lend, $rend);

                    my $percent_feature = $nucs_common / $feature_length * 100;
                    
                    my $percent_feat = $nucs_common / $feat_length * 100; ## yes, should use something other than feature and feat for discrimination facilitation.
                    if ($percent_feat > $max_percent_overlap_allowed || $percent_feature > $max_percent_overlap_allowed) {
                        next tiers;
                    }
                }
            }
            
            # if got here, no overlap in current tier.  Just add it:
            push (@$tier, $feature);
            $tiered_feature_flag = 1;
            last tiers;
        }
        
        unless ($tiered_feature_flag) {
            # no current tier can accommodate it.  Add another tier with this element
            unless (defined($max_tier) && scalar (@{$self->{tiers}} >= $max_tier) ) {
                push (@{$self->{tiers}}, [$feature]);
            }
        }
    }
    
    ## return tiers:

    return (@{$self->{tiers}});
    
}


# private
sub _nucs_in_common {
    my $self = shift;
    
    my ($e5, $e3, $g5, $g3) = @_;
    ($e5, $e3) = sort {$a<=>$b} ($e5, $e3);
    ($g5, $g3) = sort {$a<=>$b} ($g5, $g3);
    my $length = abs ($e3 - $e5) + 1;
    my $diff1 = ($e3 - $g3);
    $diff1 = ($diff1 > 0) ? $diff1 : 0;
    my $diff2 = ($g5 - $e5);
    $diff2 = ($diff2 > 0) ? $diff2 : 0;
    my $overlap_length = $length - $diff1 - $diff2;
    return ($overlap_length);
}


package TierFeatures::Feature;
use strict;
use warnings;


sub new {
    my $packagename = shift;
    my ($lend, $rend, $feat_id) = @_;
    
    if ($lend > $rend) {
        ($lend, $rend) = ($rend, $lend);
    }

    my $self = { lend => $lend,
                 rend => $rend,
                 feat_id => $feat_id,
             };

    bless ($self, $packagename);
    return ($self);
}
                
1; #EOM
