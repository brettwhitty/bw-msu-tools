#!/usr/bin/perl

use strict;
use warnings;

use Devel::Size qw{ total_size };
use DB_File;
use Bio::DB::Fasta;
use File::Temp qw{ tempfile };
use File::Copy qw{ copy };

$| = 1;

my $prefix = shift @ARGV || die "Provide PUTs prefix";

my $members_txt_file = "$prefix.PUT_member.txt";
my $subsequence_txt_file = "$prefix.Sub_Sequence.txt";
my $members_fasta_file = "$prefix.mRNA.PUTmemberSequence.fasta";
my $subsequence_fasta_file = "$prefix.mRNA.Subsequence.fasta";
my $puts_fasta_file = "$prefix.mRNA.PUT.fasta";

my $mem = {};
my $in_members_file = {};

print STDERR "Assigning PUT memberships...\n";
open (IN, $members_txt_file) || die "$!";
while (<IN>) {
    chomp;
    
    my @t = split(/\t/);

    $mem->{$t[2]}->{$t[1]} = 1;

    $in_members_file->{$t[2]} = 1;
}

#print "".(total_size($mem)/1000)."KB\n";

my @nonmember = ();
#my @nonunique;
#my $smem = {};

my $zero_count = 0;
do {
    my $mput = {};
    @nonmember = ();

    open (IN, $subsequence_txt_file) || die "$!";
    while (<IN>) {
        chomp;

        my @t = split(/\t/);
    
        unless(defined($mput->{$t[0]})) {
            $mput->{$t[0]} = {};
        }
    
        if (defined($mem->{$t[1]})) {
            #$mput->{$t[0]}->{$mem->{$t[1]}} = 1;
            foreach my $k(keys(%{$mem->{$t[1]}})) {
                $mput->{$t[0]}->{$k} = 1;
            }
        }
    }


    foreach my $key(keys(%{$mput})) {
        my $put_count = scalar(keys(%{$mput->{$key}}));

        if ($put_count >= 1) {
            #          push(@nonunique, $key);
            #          print STDERR "Nonunique substring '$key' removed\n";
            #          delete $mput->{$key};
            #} elsif ($put_count == 1) {
            #($mem->{$key}) =
            foreach my $k(keys(%{$mput->{$key}})) {
                $mem->{$key}->{$k} = 1;
            }
        } elsif ($put_count < 1) {
            push(@nonmember, $key);
        }
    }
    
#    print STDERR scalar(@nonmember)." nonmembers\n";
    
    ## it might be necessary to loop one extra time to pick up mappings
    ## for the last set of subsequence members added to the hash
    ## so that's what's going on here
    if (scalar(@nonmember) == 0) {
        $zero_count++;
    }
    
} while ($zero_count <= 1);


my $puts = {};
foreach my $gi(keys(%{$mem})) {
#    if (defined($mem->{$gi}->{'PUT-157a-Solanum_tuberosum-75773149'})) {
#        print `xdget -n Solanum_tuberosum.mRNA.Subsequence.fasta $gi 2>/dev/null`;
#    }
    foreach my $put_id(keys(%{$mem->{$gi}})) {
        push (@{$puts->{$put_id}}, $gi);
    }
}

#unless (-e "$members_fasta_file.xni") {
#    system("xdformat -I -n $members_fasta_file");
#}
#unless (-e "$subsequence_fasta_file.xni") {
#    system("xdformat -I -n $subsequence_fasta_file");
#}

print STDERR "DONE\n";


print STDERR "Preparing sequence database files...\n";

my $mem_db = new Bio::DB::Fasta($members_fasta_file);
my $sub_db = new Bio::DB::Fasta($subsequence_fasta_file);
my $puts_db = new Bio::DB::Fasta($puts_fasta_file);

tie my %members_seq, "DB_File", "$prefix.members.db" ;
tie my %subsequence_seq, "DB_File", "$prefix.subsequence.db" ;

my (undef, $con_temp) = tempfile(UNLINK => 1, OPEN => 0);
my (undef, $sub_temp) = tempfile(UNLINK => 1, OPEN => 0);

foreach my $put_id(keys(%{$puts})) {
#    print "## $put_id\n";
    
    my $out_fh;
    
    my $consensus = $puts_db->get_Seq_by_id($put_id);
    my $consensus_seq = $consensus->seq();
    $consensus_seq =~ s/(\S{1,60})/$1\n/g;
    
    open(my $con_fh, ">$con_temp") || die "Failed to open consensus sequence file '$con_temp' for writing: $!";
    
    print $con_fh ">consensus\n$consensus_seq";
    close $con_fh;
   
    
    open(my $sub_fh, ">$sub_temp") || die "Failed to open subsequence file '$sub_temp' for writing: $!";
    
    my $mem_fsa = '';
    my $sub_seq = {};
    my $sub_orient = {};
    my $sub_count = 0;
    foreach my $gi(@{$puts->{$put_id}}) {
        my $member = (defined($in_members_file->{$gi})) ? 1 : 0;
#        print join("\t", ($put_id, $gi, $member))."\n";
        if ($member) {    
            my $seq_obj = $mem_db->get_Seq_by_id($gi);
            my $seq = $seq_obj->seq();
            $seq   =~ s/(\S{1,60})/$1\n/g;
            $mem_fsa .= ">$gi\n$seq";
        } else {
            ## in the subsequence file the gi numbers are prefixed with "gi|"
            my $seq_obj = $sub_db->get_Seq_by_id("gi|$gi");
            my $seq = $seq_obj->seq();
            $sub_seq->{$gi} = $seq;
            $seq   =~ s/(\S{1,60})/$1\n/g;
            print $sub_fh ">$gi\n$seq";
            $sub_count++;
        }
    }
    $members_seq{$put_id} = $mem_fsa;
   
    if ($sub_count > 0) { 
        ## deal with subsequences
        my (undef, $nucmer_prefix) = tempfile(UNLINK => 0, OPEN => 0);

        my $err = system("nucmer $con_temp $sub_temp -p $nucmer_prefix 2>/dev/null");    
        if ($err) { die "Failed running nucmer";}

        ## show-coords output 
        my (undef, $show_coords_temp) = tempfile(UNLINK => 0, OPEN => 0);
        $err = system("show-coords -cdTH $nucmer_prefix.delta 1>$show_coords_temp 2>/tmp/show-coords.log");
        if ($err) {
            die "Failed running show-coords on '$nucmer_prefix.delta': $err";
        }

        open(my $show_coords_fh, $show_coords_temp) || die "Failed opening show-coords output '$show_coords_temp' for reading: $!";
        my $gi_counter = {};
        
        ## CHANGE THIS PART OF THE CODE TO VERIFY THAT WE GOT ALL OF THE GIs THAT WE EXPECTED IN THE OUTPUT
        ## AND THAT THERE WERE NONE THAT HAVE 2 OR MORE ROWS IN THE TABLE
        ##
        ## JUST EXCLUDE ANY SUBSEQUENCE THAT CAN'T BE MAPPED BY NUCMER TO THE CONSENSUS
        
        while (<$show_coords_fh>) {
            
            #print STDERR $_;
            chomp;

            my @col = split("\t", $_);
            
            my $gi = $col[12];
            
            ## for debugging
##            $gi_counter->{$gi}++;
##            if ($gi_counter->{$gi} > 1) {
##                print ">1 alignment for gi '$gi' in $nucmer_prefix.delta\n";
##                die();
##            }
            ## ^^^
            my $orient = '';
            if ($col[10] == -1) { 
                $orient = '-';
                #$sub_orient->{$gi} .= '-';
###                $sub_seq->{$gi} = reverse_complement_dna($sub_seq->{$gi});
                #$sub_seq->{$gi} =~ s/(\S{1,60})/$1\n/g;
                #$sub_seq->{$gi} = ">gi|$gi-\n".$sub_seq->{$gi};
            } elsif ($col[10] == 1) {
                $orient = '+';
                #$sub_orient->{$gi} .= '+';
                #$sub_seq->{$gi} =~ s/(\S{1,60})/$1\n/g;
                #$sub_seq->{$gi} = ">gi|$gi+\n".$sub_seq->{$gi};
            } else {
                die "Something unexpected happened with column 10 in the show-coords output:\n $_";
            }
            if (! defined($sub_orient->{$gi}->{$orient})) {
                $sub_orient->{$gi}->{$orient} = $col[8];
            } elsif (defined($sub_orient->{$gi}->{$orient}) && $col[8] > $sub_orient->{$gi}->{$orient}) {
                $sub_orient->{$gi}->{$orient} = $col[8];
            }
            
        }
#        if ($line_counter != $sub_count) {
#            print "There was a mismatch in subsequence count and show-coords output for '$put_id'\n";
#            print $members_seq{$put_id};
#            my $sub_fsa;
#            foreach my $gi(keys(%{$sub_seq})) {
#                $sub_fsa .= $sub_seq->{$gi};
#            }
#            print $sub_fsa;
#            die();
#        }
    
        unlink("$nucmer_prefix.cluster");
        unlink("$nucmer_prefix.delta");
        unlink("$nucmer_prefix.coords");
        unlink("$show_coords_temp");
    }
    my $sub_fsa;
    foreach my $gi(keys(%{$sub_seq})) {
        my $orient = '';
        if (! defined($sub_orient->{$gi})) {
            print STDERR "WARNING: $put_id - Could not map subsequence '$gi' to consensus\n";
            next;
        } elsif ( scalar(keys(%{$sub_orient->{$gi}})) > 1 ) {
            #both
            $orient = ($sub_orient->{$gi}->{'+'} > $sub_orient->{$gi}->{'-'}) ? '+' : '-';
        } elsif ( $sub_orient->{$gi}->{'-'} ) {
            #reverse
            $orient = '-';
            $sub_seq->{$gi} = reverse_complement_dna($sub_seq->{$gi});
        } elsif ( $sub_orient->{$gi}->{'+'} ) {
            #forward
            $orient = '+';
        }
        $sub_seq->{$gi} =~ s/(\S{1,60})/$1\n/g;
        $sub_seq->{$gi} = ">gi|".$gi.$orient."\n".$sub_seq->{$gi};
        $sub_fsa .= $sub_seq->{$gi};
    }
    $subsequence_seq{$put_id} = $sub_fsa;
    
}

print STDERR "DONE\n";

#print "".(total_size($puts)/1000)."KB\n";
#use Data::Dumper;
#print Dumper $puts;

#print join("\n", values(%{$mem}));

#print scalar(@nonunique)."\n";
#print scalar(@nonmember)."\n";


print STDERR "Converting alignments...\n";


## alignment file conversion

my $alignment_file = "$prefix.alignment.txt";

my $msa_acc;
my $align;
my $block_counter;

open (my $in_fh, $alignment_file) || die "Failed to open alignment file '$alignment_file' for reading: $!";

tie my %alignment_db, "DB_File", "$prefix.alignment.db" ;

while (<$in_fh>) {
    chomp;
    if (/^>(\S+)/) {
        if ($align) {
            $alignment_db{$msa_acc} = process_alignment($msa_acc, $align);
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
## process the last alignment
if ($align) {
    $alignment_db{$msa_acc} = process_alignment($msa_acc, $align);
}    

print STDERR "DONE\n";

sub process_alignment {
    my ($msa_acc, $align) = @_;

    my @consensus = (@{$align->{'consensus'}});

    ## how many alignment blocks
    my $block_count = scalar(@consensus);
    
    ## commenting out removing consensus
    ##  delete($align->{'consensus'});

#    print "## $msa_acc\n";
   
    my $align_fsa = '';
    
    foreach my $acc(keys(%{$align})) {
        
        ## print the FASTA header
#        print ">$acc\n";
        $align_fsa .= ">$acc\n";        
        
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
#            print $align->{$acc}->[$i]."\n";
            $align_fsa .= $align->{$acc}->[$i]."\n";
        }
    }

#    print Dumper $align;
    return $align_fsa;
    
}


sub reverse_complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
    $r_seq = reverse($r_seq);

    return $r_seq;
}
