#!/usr/bin/perl

package CigarAligner;

use strict;
use warnings;

$| = 1;

use Carp;
use Bio::DB::Fasta;
use File::Temp qw{ tempfile };
use String::CRC32;
#use IPC::System::Simple qw{ capture };

my $genome_db;
my $est_db;
my $extend;
my $model;


sub new {
    my ($class, %args) = @_;

    my $genome_fasta    = $args{'genome_fasta'} || undef;
    my $est_fasta       = $args{'est_fasta'}    || undef;
    $extend             = (defined($args{'extend'}))    ?   $args{'extend'} :   5000;
    $model              = $args{'model'}        || 'coding2genome';

    $genome_db   = (defined($genome_fasta))
        ? new Bio::DB::Fasta($genome_fasta, -reindex => 1)
        : confess "Constructor requires 'genome_fasta' argument";

    $est_db      = (defined($est_fasta))
        ? new Bio::DB::Fasta($est_fasta, -reindex => 1)
        : confess "Constructor requires 'est_fasta' argument";


    return bless {}, $class;
}

sub do_alignment {
    my ($self, $est_id, $genome_id, $genome_start, $genome_end) = @_;

    ## start and end coordinates are 1-based coordinates like GFF3

    my $genome_len = $genome_end - $genome_start + 1;

    my ($genome_fh, $genome_fsa )   = tempfile( UNLINK => 1 ); 
    my ($est_fh,    $est_fsa    )   = tempfile( UNLINK => 1 ); 

    ## write genome subseq file
    my $genome_seq = _get_genome_subseq_fasta($genome_id, $genome_start, $genome_len, $extend);
    print $genome_fh $genome_seq;
    close $genome_fh;

    ## write est file
    my $est_seq = _get_est_fasta($est_id);
    print $est_fh $est_seq;
    close $est_fh;

    ## set up the exonerate command line 
    my $exonerate_command = 'exonerate'
                          . ' --verbose 0'
                          . ' --model est2genome'
#                          . ' --exhaustive true'
                          . ' --exhaustive false'
                          . ' --subopt false'
                          . ' --bestn 1'
#                          . ' --showtargetgff true'                          
                          . ' --showalignment false'
                          . ' --showcigar true'
#                          . ' --showcigar false'
                          . ' --showtargetgff false'                          
                          . ' --showvulgar false'
#                          . ' --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n"'
                          . ' --query ' . $est_fsa 
                          . ' --querytype dna'
                          . ' --target ' . $genome_fsa 
                          . ' --targettype dna';

    ## run exonerate
##    my $exonerate_stdout = capture($exonerate_command);
    my $exonerate_stdout = `$exonerate_command`;
  
    chomp $exonerate_stdout;

    $exonerate_stdout =~ /^(cigar: \S+ \S+ \S+ \S+ )(\S+) (\S+) (\S+)(.*)$/;
    my ($left, $genome_seg, $seg_start, $seg_end, $right) = ($1, $2, $3, $4, $5);

    $genome_seg =~ /^[^:]+:subseq\((\d+),(\d+)\)/ || confess "$genome_seg";
    my ($offset, $len) = ($1, $2);
    my $new_genome_start = $offset + $seg_start - 1;
    my $new_genome_end   = $offset + $seg_end   - 1;
#    my $genome_len = $genome_end - $genome_start + 1;

    my $cigar = $left . "$genome_id $new_genome_start $new_genome_end" . $right; 

    ## remove temp files
    File::Temp::cleanup();

    unlink($genome_fsa);
    unlink($est_fsa);

    #my $exon_coords = _process_exonerate(\$exonerate_stdout);
    my ($genome_fasta, $cigar_line) = _process_exonerate_cigar($cigar);

#    my $exon_seqs = '';
#    foreach my $exon(@{$exon_coords}) {
#        my $seq = _get_genome_subseq_fasta($genome_id, $exon->[0], $exon->[1], undef, 1);
#        $exon_seqs .= $seq;
#    }

#    return $exon_seqs;

    return ($genome_fasta, $cigar_line);
}

sub _get_est_fasta {
    my ($est_id) = @_;
    
    my $seq = $est_db->seq($est_id);
    
    $seq =~ s/(.{1,60})/$1\n/g;
    
    return ">$est_id\n$seq";
}

sub _get_genome_subseq_fasta {
    my ($genome_id, $start, $subseq_len, $extension, $digest) = @_;

    my $seq_len = $genome_db->length($genome_id);   ## length of sequence

    my $end = $start + $subseq_len - 1;             ## subsequence end position
    
    ## extend the subsequence region if $extension is != 0
    if (defined($extension)) {
        ## adjust start
        $start -= $extension;
        if ($start < 1) {
            $start = 1;
        }
        ## adjust end
        $end += $extension;
        if ($end > $seq_len) {
            $end = $seq_len;
        }
    }
    $subseq_len = $end - $start + 1;  ## recalculate adjusted subseq length

    my $seq = $genome_db->seq($genome_id, $start, $end);

    $seq =~ s/(.{1,60})/$1\n/g;

    my $seq_id = "$genome_id:subseq($start,$subseq_len)";

    if ($digest) {
        my $short_id = 'g'.uc(sprintf("%x", crc32($seq_id)));
        return ">$short_id $seq_id\n$seq";
    } else {
        return ">$seq_id\n$seq";
    }
}

sub _process_exonerate_cigar {
    my ($cigar_line) = @_;

    $cigar_line =~ /^cigar: (\S+) (\d+) (\d+) ([+-]) (\S+) (\d+) (\d+) ([+-]) (\d+)\s+(.*)/
        || die "Cigar line seems to be in unexpected format\n$cigar_line\n";

    my ($q_id, $q_start, $q_end, $q_orient, $s_id, $s_start, $s_end, $s_orient, $unknown, $rest)
     = ($1,    $2,       $3,     $4,        $5,    $6,       $7,     $8,        $9,       $10);

    ## flip subject start end if end < start
    if ($s_start > $s_end) {
        ($s_end, $s_start) = ($s_start, $s_end);
    }
    ## flip query start end if end < start
    if ($q_start > $q_end) {
        ($q_end, $q_start) = ($q_start, $q_end);
    }

    ## adjust start coords from interbase
    $q_start += 1;
    $s_start += 1;

    my $query_len        = $est_db->length($q_id);
    my $query_full_seq   = $est_db->seq($q_id);
    my $query_seq        = $est_db->seq($q_id, $q_start, $q_end);
    if ($q_orient eq '-') {
        $query_seq = _reverse_complement_dna($query_seq);
    }
    my $subject_len = $genome_db->length($s_id);
    my $subject_seq = $genome_db->seq($s_id, $s_start, $s_end);
    if ($s_orient eq '-') {
        $subject_seq = _reverse_complement_dna($subject_seq);
    }

    ## calculate padding
    my $start_gap_len = $q_start - 1;
    my $start_gap   = '-' x $start_gap_len;
    my $end_gap_len   = $q_end - $query_len;
    my $end_gap     = '-' x $end_gap_len;

    my @cigar = split(" ", $rest);

    my $q_aln;
    my $s_aln;
    while (@cigar) {
        my $type = shift @cigar;
        my $len = shift @cigar;

        if ($type eq 'M') {
            $q_aln .= substr($query_seq, 0, $len);
            $query_seq = substr($query_seq, $len);
            $s_aln .= substr($subject_seq, 0, $len);
            $subject_seq = substr($subject_seq, $len);
        } elsif ($type eq 'D') {
            $subject_seq = substr($subject_seq, $len);
        } elsif ($type eq 'I') {
            $q_aln .= substr($query_seq, 0, $len);
            $query_seq = substr($query_seq, $len);
            $s_aln .= '-' x $len;
        } else {
            die "Unexpected cigar character '$type'";
        }
    }
    $s_aln = $start_gap . $s_aln . $end_gap;
    $q_aln = $start_gap . $q_aln . $end_gap;

    return ($s_aln, $cigar_line);
}


sub _process_exonerate_coords {
    my ($out_ref) = @_;

    my $exon_coords = [];
    my $ryo;
    while (${$out_ref} =~ /(.*)\n/g) {
        my $line = $1;
      
        if ($line =~ /\texon\t/) {
            my @gff = split("\t", $line);
            my ($genome_seg, $seg_start, $seg_end) = ($gff[0], $gff[3], $gff[4]);

            $genome_seg =~ /^[^:]+:subseq\((\d+),(\d+)\)/ || confess "$genome_seg";
            my ($offset, $len) = ($1, $2);

            ## map segment coordinates back to genome coords
            my $genome_start = $offset + $seg_start - 1;
            my $genome_end   = $offset + $seg_end   - 1;
            my $genome_len = $genome_end - $genome_start + 1;

            if ($genome_len < 0) {
                confess "genome_len calculated to be < 0, this should not happen";
            }
            if ($genome_len == 0) {
                confess "genome_len calculated to be == 0, this should not happen";
            }

            push(@{$exon_coords}, [$genome_start, $genome_len]);
        }
        ## last line is ryo, so will exit loop with this value
    }
   
    return $exon_coords; 
}

sub _reverse_complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb-/TtGgCcAaKkYyWwSsRrMmBbDdHhVv-/;
    $r_seq = reverse($r_seq);

    return $r_seq;
}

1;
