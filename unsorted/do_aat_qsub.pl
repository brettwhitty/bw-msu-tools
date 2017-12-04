#!/usr/bin/perl

## script to do qsub for AAT searches and check that they finish

use lib '/mnt/lustre/whitty/perl5lib/lib/perl5/site_perl/5.8.8';
use lib '/mnt/home/whitty/lib';

use strict;
use warnings;
use Carp;
use Getopt::Long;

use Cwd qw( abs_path );
use File::Basename;
use File::Path;
#use Data::GUID;
#my $guid = new Data::GUID();

## user name
my $username = 'whitty';
## maximum number of queued jobs allowed
my $max_jobs = 256;
## queue name
my $queue_name = 'brody_main';
## number of retries when doing qsub
my $retries = 5;
## time to sleep between qsubs and retries
my $sleep_time = 2;
## number of seconds qsub will block if it can't contact submission server
## this is an attempt to stop the qsub hangs, but it doesn't seem to work
my $qsub_block_seconds = 1;

my $input_file;
my $work_dir;
my $output_dir;
my $database;
my $database2;
my $protein = 0;
my $force_rerun = 0;
my $force_split = 0;
my $identity = 0;
my $similarity = 0;
my $coverage = 0;
my $chains = 10;
my $type = 'n';
my $help;

my $results = GetOptions(
                            "input|i=s"             =>  \$input_file,
                            "work_dir|w=s"          =>  \$work_dir,
                            "output_dir|o=s"        =>  \$output_dir,
                            "database|d=s"          =>  \$database,
                            "database2|d2=s"        =>  \$database2,
                            "type|t=s"              =>  \$type,
                            "max_jobs|m=i"          =>  \$max_jobs,
                            "force_rerun!"          =>  \$force_rerun,
                            "force_split!"          =>  \$force_split,
                            "identity|p=i"          =>  \$identity,
                            "similarity|s=i"        =>  \$similarity,
                            "chains|c=i"            =>  \$chains,
                            "coverage=i"            =>  \$coverage,
                            "help|h!"               =>  \$help,
                        );
if ($help) {
    print "Usage:\ndo_aat_qsub.pl -i /path/to/input_multifasta_assemblies.fsa -o /path/to/output_dir -w /path/to/work_dir -d /path/to/ext_database.fsa \\\n [-d2 /path/to/gap2_or_nap_database] [--protein] [--max_jobs 100] [--force_rerun]\n";
    exit(1);
}

## set protein flag based on type flag
$protein = ($type eq 'p') ? 1 : 0;

## this string will be used later for checking expected output files
my $search_type = ($protein) ? 'nap' : 'gap2';

## set shell scripts to use
my $shell_scripts = ($protein) ?
                    { masked =>  '/home/whitty/aat_aa_masked.sh', unmasked =>  '/home/whitty/aat_aa_unmasked.sh' }
                              :
                    { masked =>  '/home/whitty/aat_na_masked.sh', unmasked =>  '/home/whitty/aat_na_unmasked.sh' };
                    
unless ($output_dir) {
    confess "You must specify an output dir with --output_dir";
}
if (! -d $output_dir) {
    #confess "Specified output dir '$output_dir' doesn't exist";
    mkpath([$output_dir]);
}
unless ($work_dir && -d $work_dir) {
    confess "Must provide full path to a working directory with --work_dir";
}    
unless ($input_file) {
    confess "Must provide an input file with --input";
}
unless (-f $input_file) {
    confess "Specified input file '$input_file' doesn't exist";
}
unless ($database) {
    confess "Must provide a database to search with --database";
}
unless (-f $database) {
    confess "Specified database file '$database' doesn't exist";
}
if ($database2 && ! -f $database2) {
    confess "File provided with --database2 flag doesn't exist";
} 
unless ($database2) {
    $database2 = $database;
}

## set shell script to use based on masked/unmasked input file
my $shell_script = ($input_file =~ /\.masked$/) ? $shell_scripts->{'masked'} : $shell_scripts->{'unmasked'};

## strip trailing slash from work_dir
$work_dir =~ s/\/$//;

## expand output dir to its abs_path (this also strips trailing slashes)
$output_dir = abs_path($output_dir);

## we will store qsub related files here
my $qsub_dir = "$output_dir/.qsub";

## create the qsub_dir directory if it doesn't exist
unless (-d $qsub_dir) {
    mkpath([$qsub_dir]) || confess "Failed to make .qsub directory '$qsub_dir': $!";
}

my $input_file_dirname = dirname($input_file);
my $input_file_basename = basename($input_file);

## the directory to store single sequence fasta files generated from the input fasta file
my $split_input_dir = "$work_dir/$input_file_basename";

## get the list of split input files
my @input_file_list = get_input_file_list($input_file, $split_input_dir);

## set the number of jobs we can submit
my $free_jobs = free_jobs();

my $qsub_count = 0;
foreach my $iter_input_file(@input_file_list) {

    ## if we've reached the limit for free slots then we're done for now
    if ($qsub_count >= $free_jobs) {
        print STDERR "Maximum number of jobs reached, rerun script later to continue\n";
        last;
    }
    
    my $iter_file_prefix = get_input_file_prefix($iter_input_file);
    
    ## strip .masked off iter_input_file as the shell script doesn't want it
    $iter_input_file =~ s/\.masked$//;

    ## check if the job has already completed by looking for expected output file
    my $expected_output_file = "$output_dir/$iter_file_prefix.$search_type.btab";
    ## can over-ride with force_rerun flag
    if (! $force_rerun && -e $expected_output_file) {
        ## if so, skip submitting this job
        print STDERR "Output file '$expected_output_file' already exists, skipping qsub\n";
        next;
    }

    ## so the job hasn't completed yet, but has it been submitted before?
    my $job_id_file = get_job_id_file($iter_file_prefix);
    my $job_id = get_job_id_from_file($job_id_file);

    ## if it's in the queue, we're not gonna rerun it
    if (! in_queue($job_id)) {
        for (my $i = 0; $i < $retries; $i++) {
            if ($i > 0) {
                print STDERR "qsub retry #$i\n";
                sleep($sleep_time);
            }
            ## do the qsub
            my $job_id = do_qsub($iter_input_file, $iter_file_prefix, $database, $database2, $output_dir, $queue_name);
            
            ## if we didn't get a job id, we're going to try and see if we can track it down, and if not we'll retry the qsub
            unless ($job_id) {
                ## try to get the job id from the file if it was written
                my $job_id_file = get_job_id_file($iter_file_prefix);
                $job_id = get_job_id_from_file($job_id_file);
                ## unless we were able to get the job id, do a retry on the qsub
                unless ($job_id) {
                    next;
                }
            }
            if (-e $expected_output_file || in_queue($job_id)) {
                print "Job '$job_id' was submitted\n";
                ## increment the qsub count unless the job completed already
                unless (! $force_rerun && -e $expected_output_file) {
                    $qsub_count++;
                }
                last;
            }   
        }
    } else {
        print STDERR "Job for '$iter_input_file' is already in the queue, job id '$job_id'\n";
    }
        
}

## check if an input file has already been qsub'd, and if so is it still in the queue
## will return true (1) only if the job has been qsubbed and is still in the queue
sub in_queue {
    my ($job_id) = @_;

    unless ($job_id) {
        return 0;
    }
    
    my $status = `qstat_retry -f1 $job_id 2>&1`;

    if ($status =~ /Unknown Job Id/) {
        return 0;
    } else {
        return 1;
    }
}   

## fetch job id from file
sub get_job_id_from_file {
    my ($job_id_file) = @_;

    my $job_id = '';
    
    unless (-f $job_id_file) { 
        return '';
    }
    
    ## read the previous job id from the $job_id_file
    open (IN, $job_id_file) || confess "Couldn't open job id file '$job_id_file' for reading: $!";
    $job_id = <IN>;
    close IN;
   
    if (defined($job_id)) { ## so we don't get any complaints 
        chomp $job_id;
    }
    
    return $job_id;
}

## does the qsub
sub do_qsub {
    my ($input_file, $input_file_prefix, $database, $database2, $output_dir, $queue_name) = @_;
 
    ## this will be used to track the job id
    my $job_id_file = get_job_id_file($input_file_prefix); 
  
    ## job stdout will go here
    my $job_stdout = get_job_stdout($input_file_prefix);
    
    ## job stderr will go here
    my $job_stderr = get_job_stderr($input_file_prefix);
  
    my $job_id = '';
    
    while ($job_id eq '') {
    
        eval {
            local $SIG{ALRM} = sub {die "qsub isn't working"};
           
            alarm 2;
            ## do the qsub and save the job_id 
            $job_id = `qsub -b $qsub_block_seconds -o $job_stdout -e $job_stderr -q "$queue_name" -v "INPUT_FILE=$input_file,INPUT_FILE_PREFIX=$input_file_prefix,DATABASE=$database,OUTPUT_DIR=$output_dir,CHAINS=$chains,PCT_SIMILARITY=$similarity,PCT_IDENTITY=$identity" $shell_script | tee $job_id_file`;
            alarm 0;        
        };
        unless($@) {
            last;
        }
    } 
    chomp $job_id;
    
    return $job_id; 
}

## returns a count of a user's jobs in any state currently in the queue
sub get_job_count {
    my $count = `qstat_retry 2>&1 | grep $username | grep $queue_name | wc -l`;
    return $count;
}

## returns the number of open slots a user has avaiable for qsubbing
sub free_jobs {
    my $free_jobs = $max_jobs - get_job_count();

    if ($free_jobs < 0) {
        $free_jobs = 0;
    }

    return $free_jobs;
}

## splits an input fasta file and returns a list of files created
sub split_input {
    my ($input_file, $output_dir) = @_;

    my @file_list = `split_fasta.pl -i $input_file -o $output_dir -c 1 -v | tee $output_dir/.list`;
    chomp @file_list;

    ## also split the unmasked input file if input was a masked input file
    if ($input_file =~ /^(.*)\.masked$/) {
        my $unmasked_input_file = $1;
        unless (-f $unmasked_input_file) {
            confess "Unmasked input file '$unmasked_input_file' doesn't exist";
        }
        system("split_fasta.pl -i $unmasked_input_file -o $output_dir -c 1");
        if ($?) {
            confess "Failed to split unmasked input file '$unmasked_input_file'";
        }
    } 
    
    return @file_list;
}

## parses out the input file prefix
sub get_input_file_prefix {
    my ($input_file) = @_;

    my $prefix = basename($input_file);
    $prefix =~ s/(\.[^\.]+(\.masked)?)$//; ## remove file suffix
                
    return $prefix; 
}

## gets the list of input files
sub get_input_file_list {
    my ($input_file, $fasta_split_dir) = @_;

    my @list;
    
    my $list_file = "$fasta_split_dir/.list";
    
    if (! $force_split && -f $list_file) {
        open (IN, $list_file) || confess "Failed to open list file '$list_file' for reading: $!";
        while (<IN>) {
            chomp;
            if (-f $_) {
                push (@list, $_);
            } else {
                confess "File '$_' specified in list file '$list_file' doesn't actually exist";
            }
        }
    } else {
        unless (-d $fasta_split_dir) {
            mkpath([$fasta_split_dir]) || confess "Failed to make directory '$fasta_split_dir': $!";
        }
        @list = split_input($input_file, $fasta_split_dir);
    }
   
    return @list; 
}

sub get_job_id_file {
    my ($input_file_prefix) = @_;
    
    my $job_id_file = "$qsub_dir/$input_file_prefix.job_id";

    return $job_id_file;
}

sub get_job_stdout {
    my ($input_file_prefix) = @_;
    
    my $out = "$qsub_dir/$input_file_prefix.stdout";
    
    return $out;
}

sub get_job_stderr {
    my ($input_file_prefix) = @_;
    
    my $err = "$qsub_dir/$input_file_prefix.stderr";

    return $err;
}
