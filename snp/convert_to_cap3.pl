#!/usr/bin/perl

use strict;
use warnings;

use DB_File;

my $prefix = shift @ARGV || die;
my $put_id = shift @ARGV || '';

tie my %aln, 'DB_File', "$prefix.full_alignment.db", O_RDONLY;

my @put_ids = keys %aln;

#print scalar(@put_ids)."\n";

print_fake_cap3_header();
dump_put_members();
print "\nDETAILED DISPLAY OF CONTIGS\n";
dump_alignments();


sub dump_put_members {
    foreach $put_id(@put_ids) {
        print_contig_header($put_id);
        $_ = $aln{$put_id};
        while (/>(\S+)\n/g) {
            my $acc = $1;
            unless ($acc eq 'consensus') {
                print $acc."\n";
            }
        }
    }
}

sub dump_alignments {
    foreach $put_id(@put_ids) {
        print_contig_header($put_id);
        
        my $seqs = {};
        my $consensus = '';
        $_ = $aln{$put_id};

        while (/>(\S+)\n([^>]+)/g) {
            my $acc = $1;
            my $seq = $2;
            $seq =~ s/\s+//g;
            if ($acc eq 'consensus') {
                $consensus = [split("", $seq)];
            } else {
                #$seq =~ s/^-+|-+$/ /g;
                
                ## too tired to come up with a less verbose solution to this
                if ($seq =~ /^(-+)/) {
                    my $leading_dash = $1;
                    $leading_dash =~ tr/-/ /;
                    $seq =~ s/^-+/$leading_dash/;
                }
                if ($seq =~ /(-+)$/) {
                    my $lagging_dash = $1;
                    $lagging_dash =~ tr/-/ /;
                    $seq =~ s/-+$/$lagging_dash/;
                }
                
                $seqs->{$acc} = [split("", $seq)];
            }
           #print "$acc:\n$seq";
#            unless ($acc eq 'consensus') {
#                print $acc."\n";
#            }
            
        }
        my $block_count = 0;
        while (scalar(@{$consensus}) > 0) {
            $block_count++;
            print_dot_line();
            foreach my $acc(keys(%{$seqs})) {
                if (scalar(@{$seqs->{$acc}}) > 0) {
                    my $seq_line = join("", splice(@{$seqs->{$acc}}, 0, 60));
#                    if ($block_count == 1) {
#                        $seq_line =~ s/^-+|-+$//g;
#                    $seq_line =~ s/^-+|-+$//g;
                    print_seq_line($acc, $seq_line);
                }
            }
            print_bottom_line();
            print_seq_line('consensus', join("", splice(@{$consensus}, 0, 60)) );
            print "\n";
        }
    }

    ## print and extra line after consensus
}

sub print_fake_cap3_header {
    print <<END;
Number of segment pairs = 12345; number of pairwise comparisons = 6789
'+' means given segment; '-' means reverse complement

Overlaps            Containments  No. of Constraints Supporting Overlap    

END
}

sub print_contig_header {
    my ($id) = @_;

    print <<END;
******************* $id ********************
END
}

sub print_dot_line {
    print <<END;
                          .    :    .    :    .    :    .    :    .    :    .    :                          
END
}

sub print_bottom_line {
    print <<END;
                      ____________________________________________________________
END
}

sub print_seq_line {
    my ($id, $seq) = @_;

    unless ($seq =~ /^\s+$/) {
        print $id.(' ' x (22 - length($id)))."$seq\n";
    }
}
