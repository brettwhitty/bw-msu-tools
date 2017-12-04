#!/opt/rocks/bin/perl

use strict;
use warnings;

## This script takes a mirror of a PlantGDB PUT file set for a given species
## and uses vmatch to map the excluded near-exact subsequence files to an
## optimal PUT member sequence aligment, then adds them into the MSAs present
## in the *.alignment.txt file
##
## It requires the following files from PlantGDB:
##
## *.mRNA.PUTmemberSequence.fasta   a fasta file of genbank EST sequences for members
##                                  where the accession is the genbank GI number
##
## *.mRNA.Subsequence.fasta         a fasta file of genbank EST sequences for the
##                                  excluded near-exact subsequences
##
## *.PUT_member.txt                 table of PUT members used to determine
##                                  which GI #s correspond to which PUTs, and to
##                                  identify singletons
##
## *.alignment.txt                  PUT member CAP3 multiple sequence alignment
##                                  that was used to determine consensus sequence
##
## It produces the following file:
## 
## *.alignment.with_subseq.txt      modified alignment file that will contain
##                                  the additional member subsequence files

use Cwd qw{ abs_path getcwd };
use File::Temp qw{ tempfile };
use File::Basename qw{ basename dirname };
use Bio::DB::Fasta;

my $put_fasta = shift @ARGV || die "Please provide full path to a PUTs fasta file";
#my $prefix = shift @ARGV || die "Please provide PUT file prefix";

$put_fasta = abs_path($put_fasta);
my $base = basename($put_fasta);
my $dir = dirname($put_fasta);

chdir($dir);

$base =~ /^([^.]+)/;
my $prefix = "$dir/$1";

#$prefix =~ s/\.$//; ## for tab completion laziness

print STDERR $put_fasta."\n"; 

## these file are expected to be available from PlantGDB for each PUT set
my $members = "$prefix.mRNA.PUTmemberSequence.fasta";
my $subseqs = "$prefix.mRNA.Subsequence.fasta";
my $putmap  = "$prefix.PUT_member.txt";
my $alignments = "$prefix.alignment.txt";

## these are new files that are created during the run
my $new_alignments = "$prefix.alignment.with_subseq.txt";
my $alignment_fasta = "$prefix.alignment.fasta";
my $vmatch_fasta = "$prefix.vmatch.fasta";

## fasta db
my $db_m = new Bio::DB::Fasta($members, -reindex => 1);
my $db_s = new Bio::DB::Fasta($subseqs, -reindex => 1);

## store put membership
my $puts = {};
## use to indentify singletons
my $member_counts = {};
open (my $in_fh, $putmap) || die "Failed to open '$putmap' for reading: $!";
while (<$in_fh>) {
    chomp;
    
    my @t = split("\t", $_);
    $puts->{$t[2]} = $t[1];
    $member_counts->{$t[1]}++;
} 

my $new_members = {};

## prepare a fasta format file of the aligned member sequences
## which will include gaps present in the MSA
## this will also have singletons
unless (-e $vmatch_fasta) { prepare_alignment_fasta(); }

#my $db_p = new Bio::DB::Fasta($puts_db);
my $db_aln = new Bio::DB::Fasta($alignment_fasta, -reindex => 1);

(undef, my $temp) = tempfile(OPEN => 0, UNLINK => 1);

## get seq stream
my $ss = $db_s->get_PrimarySeq_stream() or die;

while (my $seq = $ss->next_seq) {

    ## for storing the maximal alignment of the subsequence to a PUT member
    my $max_aln = {};

    ## write temp fasta file
    open(my $out_fh, ">$temp") || die "Failed to open '$temp' for writing: $!";
    my $id = $seq->id();
    $id =~ s/^gi\|//;
    my $query_seq = $seq->seq();
    print $out_fh ">$id\n$query_seq\n";
    close $out_fh;
    
    ## run vmatch looking for matches at 99%+ identity --- remove duplicate matches
    my @results = `/share/apps/vmatch/vmatch -showdesc 10 -d -p -l 50 -exdrop 1 -identity 99 -q $temp $vmatch_fasta | sort | uniq`;

    ## discard comment line (at the bottom from the sort)
    pop @results;

    foreach my $line(@results) {
        chomp $line;
        
        $line =~ s/^\s+|\s+$//g;
        
        my @cols = split(/\s+/, $line);
    
        ## subject_match_len  subject_id  subject_offset  D/P  query_match_len query_id  query_offset  edit_dist  eval  score  pid
        my $query_gi = $cols[5];
        my $subject_gi = $cols[1];

        my $q_aln_len = $cols[4];
        my $s_aln_len = $cols[0];

        my $q_len = $db_s->length("gi|$query_gi") or die "Unable to retrieve query length for '$query_gi' from '$subseqs'";
        my $s_len = $db_m->length($subject_gi) or die "Unabled to retrieve subject length for '$subject_gi' from '$members'";
        
        my $q_aln_pct = sprintf("%.2f", $q_aln_len / $q_len * 100);
        
        
        ## 'D' for direct or 'P' for palindrome
        my $revcomp = ($cols[3] eq 'D') ? '0' : '1';

        ## seq query fmin/fmax
        my $q_fmin = $cols[6];
        my $q_fmax = $q_fmin + $q_aln_len;

        ## set subject fmin/fmax
        my $s_fmin = $cols[2];
        my $s_fmax = $s_fmin + $s_aln_len;


        ## if query is shorter that subject, and at least 99% of it is aligned
        if ($q_len < $s_len && $q_aln_pct >= 99) {

            ## and subject sequence is a member of a PUT
            if (defined($puts->{$subject_gi})) {

                ## get put id
                my $s_put_id = $puts->{$subject_gi};

                ## get subject aligned sequence from MSA
                my $s_aln_seq_obj = $db_aln->get_Seq_by_id($subject_gi);
                my $s_aln_header  = $db_aln->header($subject_gi);
                
                my $s_aln_seq    = $s_aln_seq_obj->seq();
                my $s_aln_orient = ($s_aln_header =~ /^\S+ ([+-]+)/) ? $1 : '+';

                ## pad the query sequence so that it's position is relative
                ## to 0 on the subject sequence
                my $left_padding = $s_fmin;

                ## get the aligned region of the query seq
                my $q_aln_seq = substr($query_seq, $q_fmin, $q_fmax - $q_fmin);
                
                ## skip alignments < 50bp long
                if (length($q_aln_seq) < 50) {
                    print STDERR "*** TRIMMED TO NOTHING ***\n";
                    next; ## skip short alignments
                }
              
                ## reverse complement if we need to 
                if ($revcomp) {
                   $q_aln_seq = reverse_complement_dna($q_aln_seq);
                }

                ## add on gap positions if fmin of query > fmin of subject
                ## extra gaps will be added later if the subject sequence 
                ## has additional gaps in this region
                $q_aln_seq = '-' x $left_padding . $q_aln_seq;

                my @s_positions = split("", $s_aln_seq);
                my @q_positions = split("", $q_aln_seq);

                ## iterate through subject sequence from MSA which includes gaps
                ## and add gap positions to the query sequence
                my $q_aligned_seq = '';
                foreach my $s_position(@s_positions) {
                    if ($s_position eq '-') {
                        $q_aligned_seq .= '-';
                    } else {
                        my $q_position = shift(@q_positions) || '-';
                        $q_aligned_seq .= $q_position;
                    }
                }
                
                ## for building the sequence accession
                my $orient = ($revcomp) ? '-' : '+';

                ## find maximal length alignment among PUT members for each PUT
                if (! defined($max_aln->{$s_put_id})) {
                    $max_aln->{$s_put_id} = [ $q_aln_len, int($cols[9]), "gi|$query_gi$orient", $q_aligned_seq ]; 
                } else {
                    if ($max_aln->{$s_put_id}->[0] < $q_aln_len) {
                        $max_aln->{$s_put_id} = [ $q_aln_len, int($cols[9]), "gi|$query_gi$orient", $q_aligned_seq ]; 
                    } elsif ($max_aln->{$s_put_id}->[0] == $q_aln_len
                      &&     $max_aln->{$s_put_id}->[1] < $cols[9]) {
                        $max_aln->{$s_put_id} = [ $q_aln_len, int($cols[9]), "gi|$query_gi$orient", $q_aligned_seq ]; 
                    }
                }

            }
       }
    }
    
    ## replace -'s with space and store the new sequences
    foreach my $k(keys %{$max_aln}) {
        my $seq = $max_aln->{$k}->[3];

        $seq =~ s/^([-]*)//;
        my $left_pad = $1;
        $seq =~ s/([-]*)$//;
        my $right_pad = $1;
        
        $left_pad =~ tr/-/ /;
        $right_pad =~ tr/-/ /;
        push (@{$new_members->{$k}}, [$max_aln->{$k}->[2], $left_pad.$seq.$right_pad]);
    }

}

## add the new members to the alignments
add_new_members();
cleanup();

## remove working files
sub cleanup {
    my @files;
    @files = glob("$vmatch_fasta*");
    foreach my $file(@files) {
        unlink $file;
    }
    @files = glob("$alignment_fasta*");
    foreach my $file(@files) {
        unlink $file;
    }
}

## adds the lines for the new member sequences to the MSA
sub add_new_members {
    my $in_fh;
    my $out_fh;

    open ($in_fh, $alignments) || die "Can't open file '$alignments' for reading: $!";

    open ($out_fh, ">$new_alignments") || die "Failed to open '$new_alignments' for writing: $!";

    my $put_id = '';
    while (<$in_fh>) {
        if (/^>(\S+)/) {
            $put_id = $1;
            print $out_fh $_;
        } elsif (/^\s+\.    :    \./) {
            print $out_fh $_;
            if (defined($new_members->{$put_id})) {
                foreach my $member(@{$new_members->{$put_id}}) {
                    if (length($member->[1]) > 0) {
                        ## each line is 60 bases (or padding)
                        $member->[1] =~ s/^(.{1,60})//;
                        my $seq = $1;

                        ## if there's only whitespace (previously -'s) the line won't be shown
                        ## N.B: there's no chance of long internal gaps here because they 
                        ## won't be present in the vmatch output
                        $seq =~ s/^\s+$//g;

                        ## if there's some sequence print the line
                        if (length($seq) > 0) {
                            print $out_fh $member->[0] . " " x (22 - length($member->[0]));
                            print $out_fh $seq."\n";
                        }
                    }
                }
            }
        } else {
            print $out_fh $_;
        }
    }
    ## add in the singletons

    ## $puts hash isn't needed anymore, so it's being reversed here
    ## to use it to look up PUT ID -> GI mapping for singletons
    %{$puts} = reverse %{$puts};

    foreach my $put_id(keys(%{$member_counts})) {
        if ($member_counts->{$put_id} == 1) {
            if (defined($new_members->{$put_id})) {
                my $gi = $puts->{$put_id};
                my $singleton_header = $db_m->header($gi);
                $singleton_header = "gi|$singleton_header+";
                my $singleton_seq    = $db_m->seq($gi);
                my $consensus_seq = $singleton_seq;
                print $out_fh ">$put_id\n";
                while ($consensus_seq =~ /(.{1,60})/g) {
                    my $con_seq = $1;
                
                    print $out_fh dot_line();

                    foreach my $member(@{$new_members->{$put_id}}) {
                        if (length($member->[1]) > 0) {
                            $member->[1] =~ s/^(.{1,60})//;
                            my $seq = $1;
                            $seq =~ s/^\s+$//g;

                            ## if there's some sequence print the line
                            if (length($seq) > 0) {
                                print $out_fh $member->[0] . " " x (22 - length($member->[0]));
                                print $out_fh $seq."\n";
                            }
                        }
                    }
                    $singleton_seq =~ s/^(.{1,60})//;
                    my $sngl_seq = $1;
                    $sngl_seq =~ s/^\s+$//g;
                    print $out_fh $singleton_header . " " x (22 - length($singleton_header));
                    print $out_fh $sngl_seq."\n";
                    print $out_fh bottom_line();
                    print $out_fh 'consensus             ';
                    print $out_fh $con_seq."\n\n";
                }
            }
        }
    }
}


sub reverse_complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
    $r_seq = reverse($r_seq);

    return $r_seq;
}

## uses vmatch to get the real coordinates of the region of the member sequence
## that is involved in the MSA, as the number in *.PUT_member.txt is incorrect
sub get_member_msa_coords {
    my ($gi) = @_;

    (undef, my $member_temp) = tempfile(OPEN => 0);
    (undef, my $aln_temp)    = tempfile(OPEN => 0);

    ## get subject aligned sequence from MSA
    my $aln_header = $db_aln->header($gi);
    my $aln_seq    = $db_aln->seq($gi);
    unless (defined($aln_header)) {
        $aln_header = $db_m->header($gi);
        $aln_seq    = $db_m->seq($gi);
    }
    $aln_seq =~ s/-//g;
    
    ## get query member sequence
    my $member_header = $db_m->header($gi);
    my $member_seq    = $db_m->seq($gi);

    open(my $outfh, ">$member_temp") || die "Failed to write temporary sequence file '$member_temp': $!";
    print $outfh ">$member_header\n$member_seq\n";
    close $outfh;
    open($outfh, ">$aln_temp") || die "Failed to write temporary sequence file '$aln_temp': $!";
    print $outfh ">$aln_header\n$aln_seq\n";
    close $outfh;

    my $temp_dir = dirname($aln_temp);
    my $cwd = getcwd();

    chdir($temp_dir);
    system("/share/apps/vmatch/mkvtree -allout -dna -pl 3 -db $aln_temp");
    chdir($cwd);

    my @results = `/share/apps/vmatch/vmatch -showdesc 10 -d -p -l 50 -exdrop 1 -identity 99 -q $member_temp $aln_temp | sort | uniq`;

    ## discard comment line (at the bottom from the sort)
    pop @results;

    my $max = [0, 0];
    foreach my $line(@results) {
        chomp $line;
        my @t = split(/\s+/, $line);
        if ($t[4] > $max->[0]) {
            $max = [$t[4], $t[6], $t[4] + $t[6]];
        }
    }
    
    ## remove temp files
    my @temp_files = glob("$aln_temp*");
    foreach my $temp_file(@temp_files) {
        unlink $temp_file;
    }
    unlink($member_temp);
    
    return [$max->[1], $max->[2]];
}


## prepares a fasta file to use for vmatch searches
sub prepare_alignment_fasta {

    ## generate fasta format MSA file
    convert_alignment($alignments, $alignment_fasta);

    ## generate a hash for put membership coords
    my %in_file = ();

    open(my $in_fh, "$alignment_fasta") || die "Failed to open '$alignment_fasta' for reading: $!";
    while (<$in_fh>) {
        if (/^>(\d+)/) {
            $in_file{$1} = 1;
        }
    }

    open(my $out_fh, ">>$alignment_fasta") || die "Failed to open '$alignment_fasta' for append: $!";

    open ($in_fh, $putmap) || die "Failed to open '$putmap' for reading: $!";
    while (<$in_fh>) {
        chomp;

        my @t = split("\t", $_);
        my $gi = $t[2];

        ## if the gi is not in the alignment fasta, fetch it and add it
        unless ($in_file{$gi}) {
            my $header = $db_m->header($gi) or die "Failed to get header for $gi";
            my $seq    = $db_m->seq($gi) or die "Failed to get seq for $gi";
            print $out_fh ">$header +\n$seq\n";
        }

    }
    close $out_fh;

    open ($out_fh, ">$vmatch_fasta") || die "Failed to open '$vmatch_fasta' for writing: $!";
   
    my $db = new Bio::DB::Fasta($alignment_fasta, -reindex => 1);
    my $ss = $db->get_PrimarySeq_stream() or die;
    while (my $seq = $ss->next_seq) {
            my $id = $seq->id();
            my $nt = $seq->seq();
            
            $nt =~ s/[-]+//g;

            print $out_fh ">$id\n$nt\n";
    }
    close $out_fh;

    system("/share/apps/vmatch/mkvtree -allout -dna -pl 6 -db $vmatch_fasta");

}

## converts alignment.txt to fasta
sub convert_alignment {
    my ($infile, $outfile) = @_;

    my $msa_acc;
    my $align;
    my $block_counter;

    open (my $in_fh, $infile) || die "Failed to open '$infile' for reading: $!";
    open (my $out_fh, ">$outfile") || die "Failed to open '$outfile' for writing: $!";
    while (<$in_fh>) {
        chomp;
        if (/^>(\S+)/) {
            if ($align) {
                process_alignment($msa_acc, $align, $out_fh);
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
        process_alignment($msa_acc, $align, $out_fh);
    }    
}


sub process_alignment {
    my ($msa_acc, $align, $out_fh) = @_;

    my @consensus = (@{$align->{'consensus'}});

    ## how many alignment blocks
    my $block_count = scalar(@consensus);
    
    delete($align->{'consensus'});

    foreach my $acc(keys(%{$align})) {
        
        my $gi = $acc;
        $gi =~ s/gi\|(\d+).*([+-])/$1 $2/;

        ## print the FASTA header
        print $out_fh ">$gi\n";
        
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
            print $out_fh $align->{$acc}->[$i]."\n";
        }
    }

}


sub dot_line {
    return <<END;
                          .    :    .    :    .    :    .    :    .    :    .    :
END
}

sub bottom_line {
    return <<END;
                      ____________________________________________________________
END
}
