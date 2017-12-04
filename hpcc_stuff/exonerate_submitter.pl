#!/usr/bin/perl

##
## This is a quick script to submit exonerate jobs to the MSU HPCC grid
## from any host
##
## It uses Net::SCP and Net::SSH::Perl to copy itself and run itself 
## across a stack of hosts in order to get to the submission host
##
## Brett Whitty
## whitty@msu.edu

use strict;
use warnings;

use Cwd qw{ abs_path };
use File::Basename;
use File::Path;
use Net::SCP;
use Net::SSH::Perl;
use Getopt::Long;
use DateTime;

#my $queue = 'bulk_brody_4s';
my $queue = 'brody_main';
my $array = '';   ## if changed, make sure to change array total also
my $array_total = 100; ## exonerate requires this
my $host = 'buell-general.plantbiology.msu.edu';
my $output_dir = '';
my $user = 'whitty';
my $walltime = '12:00:00';

my $log_dir = '/home/whitty/exonerate/log';

my $ssh_stack = 'hpcc.msu.edu,brody';

my ($query, $target, $command);

my $result = GetOptions(
                        'queue|q=s'         =>  \$queue,
                        'array_indices|t=s' =>  \$array,
                        'array_total=i'     =>  \$array_total,
                        'host|h=s'          =>  \$host,
                        'query=s'           =>  \$query,
                        'target=s'          =>  \$target,
                        'output_dir|o=s'    =>  \$output_dir,
                        'ssh_stack=s'       =>  \$ssh_stack,
                        'command=s'         =>  \$command,
                        'walltime=s'        =>  \$walltime,
                       );

#$queue = ($queue) ? $queue : $default_queue;
#$array_total = ($array_total) ? $array_total : $default_array_total;
#$host = ($host) ? $host : $default_host;
#$output_dir = ($output_dir) ? $output_dir : $default_output_dir;
#$ssh_stack = ($ssh_stack) ? $ssh_stack : $default_ssh_stack;

if ($array_total && $array eq '') {
    $array = "1-$array_total";
}
unless ($query) {
    die "Must provide path to a query file!";
}
unless ($target) {
    die "Must provide path to a target file!";
}

my $script = abs_path($0);
my $script_base = basename($script);
print STDERR $script."\n";

my @ssh_hosts = split(",", $ssh_stack);

if ($ssh_stack ne '') {

    my $ssh_host = shift @ssh_hosts;
    $ssh_stack = join(",", @ssh_hosts);
    my $scp = new Net::SCP();    
    $scp->scp($script, "$ssh_host:/tmp/$script_base");
   
    print STDERR "Connecting to '$ssh_host'\n"; 
    my $ssh = new Net::SSH::Perl($ssh_host, protocol=> '2', debug => 0);
    $ssh->login($user);
    
    $ssh->cmd("chmod +x /tmp/$script_base");
    my ($out, $err, $exit) = $ssh->cmd("/tmp/$script_base --walltime '$walltime' --queue '$queue' --array '$array' --array_total '$array_total' --host '$host' --query '$query' --target '$target' --output_dir '$output_dir' --ssh_stack '$ssh_stack' 1>out.log 2>/home/whitty/error.log");
    print "$out $err $exit\n";

    ## remove script
    $ssh->cmd("unlink /tmp/$script_base");
    
} else {

    my $dt = DateTime->now(time_zone => 'America/Detroit');
    my $timestamp = $dt->strftime("%G%m%d%H%M%S");

    my $stderr_log_dir = "$log_dir/$timestamp/stderr";
    my $stdout_log_dir = "$log_dir/$timestamp/stdout";

    mkpath("$log_dir/$timestamp", { verbose => 0 });
    mkpath("$stderr_log_dir", { verbose => 0 });
    mkpath("$stdout_log_dir", { verbose => 0 });
    
    print STDERR "$timestamp\n";
   
    my $qsub_sh = "/tmp/run_$timestamp.sh";

    open(my $sub_fh, ">$qsub_sh") || die "Failed to write shell script '$qsub_sh': $!";

    while (<DATA>) {
        print $sub_fh $_;
    }

    my $command_line = "qsub -q \"$queue\" -t \"$array\" -l \"walltime=$walltime\" -v \"REMOTE_HOSTNAME=$host,REMOTE_QUERY=$query,REMOTE_TARGET=$target,REMOTE_OUTDIR=$output_dir,PBS_ARRAYTOTAL=$array_total\" -e \"$stderr_log_dir\" -o \"$stdout_log_dir\" $qsub_sh >$log_dir/$timestamp/qsub.log";
   
    open(OUT, ">$log_dir/$timestamp.log") || die "Failed to write to '$log_dir/$timestamp.log': $!";
    print OUT $command_line."\n";
    close OUT;
    
    system("$command_line");
    
}

__DATA__
#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=2000mb
#PBS -N exonerate

## PBS -M whitty@msu.edu
## PBS -m e

if [ "$REMOTE_HOSTNAME" == '' ]
then
    echo "Must specify a REMOTE_HOSTNAME"
    exit 1
fi
if [ "$REMOTE_QUERY" == '' ]
then
    echo "Must specify a REMOTE_QUERY"
    exit 1
fi

if [ "$REMOTE_TARGET" == '' ]
then
    echo "Must specify a REMOTE_TARGET"
    exit 1
fi

if [ "$REMOTE_OUTDIR" == '' ]
then
    echo "Must specify a REMOTE_OUTDIR"
    exit 1
fi

if [ "$PBS_ARRAYTOTAL" == '' ]
then
    echo "Must specify a PBS_ARRAYTOTAL"
    exit 1
fi

if [ "$PERCENT" == '' ]
then
    PERCENT=50
fi
if [ "$MAX_INTRON" == '' ]
then
    MAX_INTRON=2000
fi

WORK_DIR="/mnt/local/$PBS_JOBID"

mkdir -p $WORK_DIR

if [ ! -d $WORK_DIR ]
then
    echo "Failed to create '$WORK_DIR'"
    exit 1
fi

REMOTE_QUERY="$REMOTE_HOSTNAME:$REMOTE_QUERY"
REMOTE_TARGET="$REMOTE_HOSTNAME:$REMOTE_TARGET"
REMOTE_OUTDIR="$REMOTE_HOSTNAME:$REMOTE_OUTDIR"

QUERY_BASE=`basename $REMOTE_QUERY`
TARGET_BASE=`basename $REMOTE_TARGET`
QUERY_FILE="$WORK_DIR/$QUERY_BASE"
TARGET_FILE="$WORK_DIR/$TARGET_BASE"

OUTPUT_FILENAME="$QUERY_BASE.$TARGET_BASE.$PBS_JOBID.$PBS_ARRAYID.exonerate.raw"
OUTPUT_FILE="$WORK_DIR/$OUTPUT_FILENAME"
REMOTE_OUTPUT_FILE="$REMOTE_OUTDIR/$OUTPUT_FILENAME"

scp $REMOTE_QUERY $QUERY_FILE
if [ ! -e $QUERY_FILE ]
then
    echo "Failed to scp '$REMOTE_QUERY' to '$QUERY_FILE'"
    exit 1
fi

scp $REMOTE_TARGET $TARGET_FILE
if [ ! -e $TARGET_FILE ]
then
    echo "Failed to scp '$REMOTE_TARGET' to '$TARGET_FILE'"
    exit 1
fi

exonerate -Q dna -T dna --model est2genome --verbose 0 --showalignment no --showsugar no --showquerygff no --showtargetgff yes --showcigar yes --showvulgar no -q $QUERY_FILE -t $TARGET_FILE --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent 50 --softmasktarget 1 --maxintron $MAX_INTRON --subopt 0 --fsmmemory 1000 --dpmemory 1000 --querychunkid $PBS_ARRAYID --querychunktotal $PBS_ARRAYTOTAL >$OUTPUT_FILE

scp $OUTPUT_FILE $REMOTE_OUTPUT_FILE

rm -rf $WORK_DIR
