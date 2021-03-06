#!/usr/bin/perl

use warnings;
use strict;

## will convert PlantGDB PUTs.alignment.txt files into fasta format MSA
## for feeding into muscle profile-profile alignments

use Data::Dumper;

my $msa_acc;
my $align;
my $block_counter;
while (<>) {
    chomp;
    if (/^>(\S+)/) {
        if ($align) {
            process_alignment($msa_acc, $align);
        }    
        $msa_acc = $1;
        $align = {};
        $block_counter = 0;
    } elsif (/^\S+/) {
        /^(.{22})(.*)/ || die "Failed to match sequence line properly";
        my ($acc, $seq) = ($1, $2);
        $acc =~ s/\s+//g;
        $seq =~ s/\s/-/g;
        $align->{$acc}->[$block_counter] = $seq;
        if ($acc eq 'consensus') {
            $block_counter++;
        }
    }
}       
if ($align) {
    process_alignment($msa_acc, $align);
}    

sub process_alignment {
    my ($msa_acc, $align) = @_;

    my @consensus = (@{$align->{'consensus'}});

    ## how many alignment blocks
    my $block_count = scalar(@consensus);
    
    delete($align->{'consensus'});

    print "## $msa_acc\n";
    
    foreach my $acc(keys(%{$align})) {
        
        ## print the FASTA header
        print ">$acc\n";
        
        for (my $i = 0; $i < $block_count; $i++) {
            ## what is the length of the alignment block
            my $block_length = length($consensus[$i]);
            
            ## add a line of -s if there was no line for this sequence in the alignment
            if (! defined($align->{$acc}->[$i])) {
                $align->{$acc}->[$i] = '-' x $block_length;
            }
            
            ## for convenience
            my $seq_len = length($align->{$acc}->[$i]);
            
            ## pad on the right with -s
            if ($seq_len < $block_length) {
                $align->{$acc}->[$i] .= '-' x ($block_length - $seq_len);
            } 
        
            ## print the line of sequence
            print $align->{$acc}->[$i]."\n";
        }
    }

#    print Dumper $align;
        
}
