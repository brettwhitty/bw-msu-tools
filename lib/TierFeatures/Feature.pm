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
                
1; 
