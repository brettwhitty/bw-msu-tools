package TierFeatures;

use strict;
use warnings;
use TierFeatures::Feature;

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

1;
