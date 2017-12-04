#!/usr/bin/perl

## This script takes a query multifasta database (eg: PUTs, uniref100) and a target
## multifasta file (eg: BACs, split chromosomes) and uses a pre-filtering strategy
## to reduce the size of the query database set before running exonerate
## on each query database and target sequence
##
## Specific exonerate flags are used for generating output consumed by downstream scripts
##
## Both the pre-filtering and exonerate steps support softmasked target sequence
## 

use strict;
use warnings;

use File::Basename;
use File::Temp qw/ :POSIX tempfile tempdir /;
use File::Path;
use Cwd qw/ getcwd abs_path /;

use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::SeqIO;

use Carp;
use File::Basename;
use POSIX qw(ceil);

use Getopt::Long;

my $max_intron_length = 2000;
my $percent_score = 50;

my $target_name;
my $query_name;
my $query_type = 'n';
my $output;
my $work_dir = getcwd();
my $filter_type = 'blast';
## TEMP ## my $nucmer_id_cutoff = 70;
## TEMP ## my $nucmer_cov_cutoff = 50;

my $result = GetOptions(
                           'target|t=s'             =>  \$target_name,
                           'query|q=s'              =>  \$query_name,
                           'query_type=s'           =>  \$query_type,
                           'output_dir|o=s'         =>  \$output,
                           'work_dir|w=s'           =>  \$work_dir,
                           'filter_type|f=s'        =>  \$filter_type,
                           'percent_score|p=i'      =>  \$percent_score,
##                           'nucmer_id_cutoff=i'     =>  \$nucmer_id_cutoff,
##                           'nucmer_cov_cutoff=i'    =>  \$nucmer_cov_cutoff,
                       );

unless ($query_type =~ /[encpENCP]/) {
    confess "Unrecognized query type";
}

if ($output) {
    $output =~ s/\/$//;
} else {
    $output = getcwd();
}

unless (-d $output) {
    mkpath($output);
}

my $model = {
            'e' =>  'est2genome',
            'c' =>  'cdna2genome',
            't' =>  'coding2genome',
            'p' =>  'protein2genome',
         };
   
my $query_type_string = {
            'e' =>  'dna',
            'c' =>  'dna',
            't' =>  'dna',
            'p' =>  'protein',
                        };    
                        
my $temp_dir = tempdir( DIR => $work_dir, CLEANUP => 1);                    
my (undef, $temp_file) = tempfile(DIR => $temp_dir, OPEN => 0);
my $prefix = basename($temp_file);

## split the target file
my @target_files = file_split($target_name, $prefix, $temp_dir, 1);
#`/home/whitty/bin/split_fasta.pl -v -i $target_name -p $prefix -c 1 -o $temp_dir`; 

my ($name, $dir) = fileparse($query_name);
$name =~ /^([^\.]+)/;
my $query_prefix = $1;

foreach my $target_file(@target_files) {
    chomp $target_file;

    open (IN, $target_file) || confess "Can't open '$target_file' for reading: $!";
    my $defline = <IN>;
    close IN;
    
    $defline =~ /\>(\S+)/ || confess "Didn't match defline in '$target_file'";
    my $target_acc = $1;
    
    system("formatdb -p F -o F -i $target_file");
                       
    $target_acc =~ s/\W/_/g;
    
    my $out_raw = "$output/$target_acc.$query_prefix.exonerate.raw";
    
    open (OUT, ">$out_raw");
    close OUT;
   
    my $tmp_db;

    if ($filter_type eq 'nucmer' || $filter_type eq 'p') {
    
        $tmp_db = do_nucmer_filter($target_file, $query_name, $query_type, $work_dir); 

    } else {
    
        $tmp_db = do_blast_filter($target_file, $query_name, $query_type, $work_dir); 

    }
    
    chomp $tmp_db;

    unless (-e $tmp_db) {
        confess "Exonerate prefilter failed to create a filtered database file: '$tmp_db'";
    }
    
    print STDERR "Running:\nexonerate -Q ".$query_type_string->{$query_type}.' -T dna --model '.$model->{$query_type}.' --verbose 0 --showalignment no --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no -q '.$tmp_db.' -t '.$target_file.' --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent '.$percent_score.' --softmasktarget 1 --maxintron '.$max_intron_length.' --subopt 0 --fsmmemory 1000 --dpmemory 1000 >'.$out_raw."\n";
    system('exonerate -Q '.$query_type_string->{$query_type}.' -T dna --model '.$model->{$query_type}.' --verbose 0 --showalignment no --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no -q '.$tmp_db.' -t '.$target_file.' --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent '.$percent_score.' --softmasktarget 1 --maxintron '.$max_intron_length.' --subopt 0 --fsmmemory 1000 --dpmemory 1000 >'.$out_raw);
    
    unlink($tmp_db);
        
}

##
## This sub pre-filters an exonerate query database
## and excludes query sequences with no significant sequence similarity to the target
## to reduce the size of the query database for the exonerate run and reduce overall run-time
##
## The BLAST parameters used also support using softmasked target sequence
##
## CAUTION:
## The first time this script is run, Bio::SeqIO will index the query database.
## Concurrent execution of the script before indexing is complete can cause corruption
## of the index file, and the script won't be able to properly fetch sequences from the
## query database when generating the output database.
sub do_blast_filter {

    my ($db_file, $query_file, $type, $work_dir) = @_;
    
    
    my $cpus = 1;
    my $cutoff = 0.01; ## this shouldn't be changed from 0.01 which was found to be lossless in test runs


    my $prog = {
            'n' =>  'tblastx',
            'c' =>  'tblastx',
            't' =>  'tblastx',
            'p' =>  'tblastn',
               };

    unless (-d $work_dir) {
        mkpath($work_dir);
    }

    my (undef, $tmp_out) = tempfile(DIR => $work_dir, OPEN => 0);

    $tmp_out = abs_path($tmp_out);
    $query_file = abs_path($query_file);

    system("blastall -p $prog->{$type} -d $db_file -i $query_file -U T -F \"m S\" -a $cpus -e $cutoff -o $tmp_out");

    my $fasta_db = new Bio::DB::Fasta($query_file);
    
    my $blast_report = new Bio::SearchIO(
                                            -format =>  'blast',
                                            -file   =>  $tmp_out,
                                        );
    
    my ($out_fh, $out_name) = tempfile( DIR => $work_dir, CLEANUP => 1 );
    
    my $fasta_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $out_fh );
    
    while( my $result = $blast_report->next_result ) {

        my $query_name = $result->query_name;

        while( my $hit = $result->next_hit ) {

            my $hit_name = $hit->name;
            my $rank = $hit->rank;
            my $score = $hit->score;
            my $sig = $hit->significance;

            my $seq_obj = $fasta_db->get_Seq_by_id($query_name);

            $fasta_out->write_seq($seq_obj);
            
        }

    }   

    unlink($tmp_out);
    
    return "$work_dir/$out_name";
}


## as an alternative to the BLAST filtering, nucmer (mummer) can be used
## for reducing the size of the query database without losing alignable sequences
##
sub do_nucmer_filter {
    
    my ($db_file, $query_file, $type, $work_dir) = @_;
    
##    print STDERR "Filtering query database using nucmer, id_cutoff = $nucmer_id_cutoff / cov_cutoff = $nucmer_cov_cutoff\n";
    print STDERR "Filtering query database using nucmer...\n";
    
    unless (-d $work_dir) {
        mkpath($work_dir);
    }
    
    my (undef, $tmp_out) = tempfile(DIR => $work_dir, OPEN => 0);
    my $nucmer_cluster = $tmp_out.".cluster";
    my $nucmer_delta = $tmp_out.".delta";

#    my ($short_seq_fh, $short_seq) = tempfile(DIR => $work_dir, OPEN => 0);
#    my ($long_seq_fh, $long_seq) = tempfile(DIR => $work_dir, OPEN => 0);
    
    $tmp_out = abs_path($tmp_out);
    $query_file = abs_path($query_file);
    
#    system("promer -p $tmp_out --maxmatch -c 9 --nooptimize $db_file $query_file 2>/dev/null"); 
#    print "Running:\nnucmer -p $tmp_out --maxmatch -l 10 -c 20 --nooptimize $db_file $query_file 2>/dev/null\n"; 
    print "Running:\nnucmer -p $tmp_out --maxmatch -l 12 -c 20 --nooptimize $db_file $query_file 2>/dev/null\n"; 
    
#    system("nucmer -p $tmp_out --maxmatch -l 10 -c 20 --nooptimize $db_file $query_file 2>/dev/null"); 
    system("nucmer -p $tmp_out --maxmatch -l 12 -c 20 --nooptimize $db_file $query_file 2>/dev/null"); 
    
    my $fasta_db = new Bio::DB::Fasta($query_file);

    my ($out_fh, $out_name) = tempfile( DIR => $work_dir, CLEANUP => 1 );
    
    my $fasta_out = Bio::SeqIO->new( -format => 'Fasta', -fh => $out_fh );
    
    #open(my $infh, $nucmer_cluster) || confess "Failed to open nucmer cluster file '$nucmer_cluster' for reading: $!";
    #open($showcoords_fh, "show-coords -ckHT $tmp_out.delta |") || die "Failed to open filehandle to show-coords on '$tmp_out.delta': $!";
    open(my $showcoords_fh, "show-coords -cHT $tmp_out.delta |") || die "Failed to open filehandle to show-coords on '$tmp_out.delta': $!";

    my %passed_seqs = ();
    
    while (<$showcoords_fh>) {
        chomp $_;
       
        my @cols = split(/\t/, $_);

## TEMP ##
####        if ($cols[6] < $nucmer_id_cutoff && $cols[8] < $nucmer_cov_cutoff) { next; }
## TEMP ##        
        
        $passed_seqs{$cols[10]} = 1;
    }
    
    foreach my $seq_id(keys(%passed_seqs)) { 
        my $seq_obj = $fasta_db->get_Seq_by_id($seq_id);

        $fasta_out->write_seq($seq_obj);
    }
    
    unlink($nucmer_delta);
    unlink($nucmer_cluster);
    
    print STDERR "Created nucmer filtered database '$out_name', ".scalar(keys(%passed_seqs))." sequences.\n";
    
    return $out_name;
}


sub file_split {

    my ($input, $prefix, $output_dir, $chunk_size) = @_;

    my $file_count = 0;
    
    unless ($chunk_size || $file_count) {
        confess "Must provide value to either --chunk_size|c or --file_count|n flags";
    }
    if ($output_dir && ! -d $output_dir) {
        confess "Specified output dir '$output_dir' doesn't exist";
    }
    unless ($output_dir) {
        $output_dir = dirname($input);
    }

    my $suffix = '';
    if (! $prefix) {
        if ($input) {
            my $basename = basename($input);
            $basename =~ s/(\.[^\.]+(\.masked)?)$//; ## remove file suffix
            ## save the suffix for later
            $suffix = $1;
            unless ($basename) {    ## in case nothing's left
                $prefix = "split_fasta";    
            }
            $prefix = $basename;
        } else {
            $prefix = "split_fasta";
        }
    }

    if ($file_count) {
        my $infh = get_input_fh($input);
        my $header_count = header_count($infh);
        $chunk_size = ceil($header_count / $file_count);
    }

    my $infh = get_input_fh($input);

    my $file_id = 0;
    my $header_count = 0;
    my $outfile;
    my $outfh;
    my @files = ();
    while (<$infh>) {
        if (/^>/) {
            $header_count++;
            if (($header_count - 1) % $chunk_size == 0) {
                $outfile = $output_dir."/".$prefix."_".++$file_id.$suffix;
                push(@files, abs_path($outfile));
                $outfh = get_output_fh($outfile);
            }
        }
        print $outfh $_;
    }

    return @files;
}

sub get_output_fh {
    my ($filename) = @_;

    open (my $outfh, ">".$filename) || confess "Can't open file '$filename' for writing: $!";

    return $outfh;
}

sub header_count {
    my ($infh) = @_;

    my $count = 0;
    while (<$infh>) {
        if (/^>/) {
            $count++;
        } else {
            next;
        }
    }

    return $count;
}


sub get_input_fh {
    
    my ($input) = @_;

    my $infh;

    if ($input) {
        open ($infh, $input) || confess "Failed opening '$input' for reading: $!";
    } else {
        open ($infh, \*STDIN) || confess "Failed opening STDIN for reading: $!";
    }

    return $infh;
}
