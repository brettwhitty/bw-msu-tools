#!/usr/bin/perl

use strict;
use warnings;

use DB_File;
use File::Temp qw{ tempfile };
use Bio::DB::Fasta;

my $prefix = shift @ARGV || die "Must provide a prefix";

my $puts_fasta = "$prefix.mRNA.PUT.fasta";

my $alignment_db = "$prefix.alignment.db";
my $members_db = "$prefix.members.db";
my $subsequence_db = "$prefix.subsequence.db";
my $full_alignment_db = "$prefix.full_alignment.db";

my $put_db = new Bio::DB::Fasta($puts_fasta);

tie my %aln_fsa, 'DB_File', $alignment_db, O_RDONLY;
tie my %members_fsa, 'DB_File', $members_db, O_RDONLY;
tie my %subsequence_fsa, 'DB_File', $subsequence_db, O_RDONLY;

## remove database if it exists
unlink($full_alignment_db);
tie my %full_aln_fsa, 'DB_File', $full_alignment_db;

## temp outfile for subsequences
(undef, my $temp_sub) = tempfile(DIR => "/scratch/tmp", OPEN => 0, UNLINK => 1);

## temp outfile for subsequence alignment
(undef, my $sub_aln) = tempfile(DIR => "/scratch/tmp", OPEN => 0, UNLINK => 1);

## temp infile for alignment
(undef, my $mem_aln) = tempfile(DIR => "/scratch/tmp", OPEN => 0, UNLINK => 1);        

## temp outfile for subsequence alignment
(undef, my $full_aln) = tempfile(DIR => "/scratch/tmp", OPEN => 0, UNLINK => 1);

foreach my $put_id(keys(%members_fsa)) {
    ## if there's no subsequences, then we're not going to realign
    if (length($subsequence_fsa{$put_id}) > 0) {
        print STDERR "## $put_id:\n";
        
        ## fetch the PUT consensus
        my $put = $put_db->get_Seq_by_id($put_id);
        my $consensus_seq = $put->seq();
        $consensus_seq =~ s/(\S{1,60})/$1\n/g;
        my $consensus = ">consensus\n".$consensus_seq;
     
        print STDERR "# consensus:\n";
        print STDERR $consensus;
        
        ## fetch the excluded subsequences
        my $subsequence = $subsequence_fsa{$put_id};
        
        print STDERR "# subsequence:\n";
        print STDERR $subsequence;
        
        
        ## temp infile for subsequence alignment
        open(my $temp_subfh, ">$temp_sub") || die "Failed opening subsequence alignment file '$temp_sub' for writing: $!";

        ## prepare the input file for subsequence alignment
        print $temp_subfh $consensus;
        print $temp_subfh $subsequence;
        close $temp_subfh;
        
#        my $muscle_pipe;
        my $err = system("muscle -in $temp_sub -out $sub_aln 2>/dev/null");
        if ($err) {
            print "There was a problem running muscle on '$temp_sub', see '/tmp/__debug__': $@\n";
            use File::Copy;
            copy($temp_sub, "/tmp/__debug__");
            die();
        }

        open(my $sub_fh, $sub_aln) || die "Failed to open subsequence alignment file '$sub_aln' for reading: $!";
        
        my $subsequence_alignment = '';
        my $skip_consensus = 0;
        while (<$sub_fh>) {
            if (/^>consensus$/) {
                $skip_consensus = 1;
            } elsif (/^>/) { 
                $skip_consensus = 0;
            }
            if ($skip_consensus) {
                next;
            }
            $subsequence_alignment .= $_;
        }
        
        ## rewrite subsequence alignment with consensus sequence removed
        open($sub_fh, ">$sub_aln") || die "Failed to open subsequence alignment file '$sub_aln' for reading: $!";
        print $sub_fh $subsequence_alignment;
        close $sub_fh;
        
        print STDERR "# subsequence alignment:\n";
        print STDERR $subsequence_alignment;
        
        ## write subsequence alignment (minus consensus sequence) to temp file
#        print $sub_alnfh $subsequence_alignment;
#        close $sub_alnfh;

#        my $subsequence_alignment_file = $temp_out;
        
        ## fetch the aligned PUT member sequences
        my $alignment = (defined($aln_fsa{$put_id})) ? $aln_fsa{$put_id} : singleton_member($put_id) . $consensus;        

        print STDERR "# alignment:\n";
        print STDERR $alignment;
        
        open(my $mem_alnfh, ">$mem_aln") || die "Failed to open alignment file '$mem_aln' for reading: $!";
        
        print $mem_alnfh $alignment;
        close $mem_alnfh;

        ## do profile-profile alignment on member seqs and subsequence alignments
        $err = system("muscle -profile -in1 $mem_aln -in2 $sub_aln -out $full_aln 2>/dev/null");
        if ($err) {
            print "There was a problem running muscle profile-profile alignment on '$mem_aln' & '$sub_aln', see '/tmp/__debug1__' & '/tmp/__debug2__': $@\n";
            use File::Copy;
            copy($mem_aln, "/tmp/__debug1__");
            copy($sub_aln, "/tmp/__debug2__");
            die();
        }
            
        open(my $full_fh, $full_aln) || die "Failed to open full alignment file '$full_aln' for reading: $!";
        #my $skip_consensus = 0;
        ## should we remove consensus here?
        my $full_alignment = '';
        #my $skip_consensus = 0;
        while (<$full_fh>) {
            #    if (/^>consensus/) {
            #    $skip_consensus = 1;
            #} elsif (/^>/) { 
            #    $skip_consensus = 0;
            #}
            #if ($skip_consensus) {
            #    next;
            #}
            $full_alignment .= $_;
        }
        print STDERR "# full alignment:\n";
        print STDERR $full_alignment;
        $full_aln_fsa{$put_id} = $full_alignment;
    } else {
        if (defined($aln_fsa{$put_id})) {
            print STDERR "## $put_id: has no subsequence, storing existing alignment\n";
            $full_aln_fsa{$put_id} = $aln_fsa{$put_id}
        } else {
            print STDERR "## $put_id: is a singleton and has no subsequence, skipping\n";
        }
    }
}

## reads PUT member fasta sequence from the members fasta hash
## modifies the fasta header to include orientation (which we will assume is forward, for singletons)
sub singleton_member {
    my ($put_id) = @_;

    my $mem_seq = $members_fsa{$put_id};
    
    $mem_seq =~ s/>(\d+)/>gi|$1\+/;

    return $mem_seq;
}
