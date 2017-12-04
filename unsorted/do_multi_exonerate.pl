#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use File::Temp qw/ :POSIX tempfile tempdir /;
use File::Path;
use Cwd;

use Getopt::Long;

my $max_intron_length = 1000;
my $percent_score = 50;

my $target_name;
my $query_name;
my $query_type = 'n';
my $output;
my $work_dir = getcwd();

my $result = GetOptions(
                           'target|t=s'     =>  \$target_name,
                           'query|q=s'      =>  \$query_name,
                           'query_type=s'   =>  \$query_type,
                           'output|o=s'     =>  \$output,
                           'work_dir|w=s'   =>  \$work_dir,
                       );

unless ($query_type eq 'n' || $query_type eq 'p') {
    die "Unrecognized query type";
}

my $model = {
#            'n' =>  'coding2genome',
            'n' =>  'cdna2genome',
            'p' =>  'protein2genome',
         };
   
my $query_type_string = {
            'n' =>  'dna',
            'p' =>  'protein',
                        };    
                        
my $temp_dir = tempdir( DIR => $work_dir, CLEANUP => 1);                    
my (undef, $temp_file) = tempfile(DIR => $temp_dir, OPEN => 0);
my $prefix = basename($temp_file);

## split the target file
my @target_files = `split_fasta.pl -v -i $target_name -p $prefix -c 1 -o $temp_dir`; 

my ($name, $dir) = fileparse($query_name);
$name =~ /^([^\.]+)/;
my $query_prefix = $1;

foreach my $target_file(@target_files) {
    chomp $target_file;

    open (IN, $target_file) || die "Can't open '$target_file' for reading: $!";
    my $defline = <IN>;
    close IN;
    
    $defline =~ /\>(\S+)/ || die "Didn't match defline in '$target_file'";
    my $target_acc = $1;
    
    system("formatdb -p F -o F -i $target_file");
                        
    my $out_raw = "$target_acc.$query_prefix.exonerate.raw";
   
    my $tmp_db = `exonerate_prefilter.pl -t $target_file -q $query_name --query_type $query_type -w $work_dir`;
    
    chomp $tmp_db;
    
    system('exonerate -Q '.$query_type_string->{$query_type}.' -T dna --model '.$model->{$query_type}.' --verbose 0 --showalignment no --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no -q '.$tmp_db.' -t '.$target_file.' --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent '.$percent_score.' --softmasktarget 1 --maxintron '.$max_intron_length.' --subopt 0 --fsmmemory 1000 --dpmemory 1000 >'.$out_raw);
    
    unlink($tmp_db);
        
}
