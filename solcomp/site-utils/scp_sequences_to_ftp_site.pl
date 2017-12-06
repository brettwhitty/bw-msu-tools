#!/usr/bin/perl

use strict;
use warnings;

## script to copy sequences to the ftp/blast server

use Cwd qw/ abs_path /;

use Net::SSH::Perl;
use Net::SCP;
use Carp;

my $scp = new Net::SCP;

my @subdirs = ( 'bacs', 'bac_ends', 'puts' );

my $ftp_globs   = {
                    'bac_ends'  =>  [
                                        "*.fsa",
                                        "*.gbk",
                                        ".version",
                                    ],
                    'bacs'      =>  [
                                        "*.fsa",
                                        "*.gbk",
                                        ".version",
                                    ],
                    'puts'      =>  [
                                        "*.PUT.fasta",
                                        "*.PUT.gff",
                                        ".version",
                                    ],
                };
my $blast_db_globs = {
                    'bac_ends'  =>  [
                                        "*.fsa",
                                        "*.fsa.*",
                                    ],
                    'bacs'      =>  [
                                        "*.fsa",
                                        "*.fsa.*",
                                    ],
                    'puts'      =>  [
                                        "*.PUT.fasta",
                                        "*.PUT.fasta.*",
                                    ],
                     };

## root directory for the FTP site                     
my $ftp_root = '/data/ftp/pub/sgr/.by_taxon_id';
## root directory for the BLAST databases
my $blast_root = '/data/blast_dbs/sol-blast/.by_taxon_id';

## local test file
my $local_test = '/tmp/scp_test_'.`date +%s`;
chomp $local_test;
## remote test file
my $remote_test = 'blast.plantbiology.msu.edu:'.$local_test;

## full remote path to 
my $remote_ftp_root = 'blast.plantbiology.msu.edu:'.$ftp_root;
my $remote_blast_root = 'blast.plantbiology.msu.edu:'.$blast_root;

my $repository = shift @ARGV || die "Provide path to a sequence repository";
my $compress = shift @ARGV; ## flag to compress files on FTP site

$repository = abs_path($repository);

## check that the repository has all the directories we expect
foreach my $subdir(@subdirs) {
    if (! (-d "$repository/$subdir")) {
        die "Subdirectory '$subdir' is missing from the repository '$repository'";
    }
}

## connect to remote server
print STDERR "Connecting to remote server...\n";
my $ssh = new Net::SSH::Perl('blast.plantbiology.msu.edu', protocol=> '2', debug => 0);
$ssh->login('whitty');

## test that remote commands can be executed
my ($out, $err, $exit) = $ssh->cmd('hostname');
chomp $out;
unless ($out =~ 'blast') {
    die "Remote command execution test failed!";
}

## test that we can SCP a file to the remote machine
`touch $local_test`;
$scp->scp("$local_test", "$remote_test") or croak "Failed to copy test file to remote machine: ".$scp->{errstr};
($out, $err, $exit) = $ssh->cmd("test -e $local_test");
if ($exit) {
    die "Failed to copy test file to remote machine";
}

## remove existing FTP directories
print STDERR "Removing existing FTP directories...\n";
$ssh->cmd("rm -rf $ftp_root/*");

## remove existing BLAST directories
print STDERR "Removing existing BLAST db directories...\n";
$ssh->cmd("rm -rf $blast_root/*");

## process each directory specified in @subdirs
foreach my $subdir(@subdirs) {
    opendir(my $dir, "$repository/$subdir") || die "Failed to open dirhandle for '$repository/$subdir': $!";

    print STDERR $subdir.":\n";
    
    my @records = readdir($dir);
    foreach my $record(@records) {
        my $local_path = "$repository/$subdir/$record";
        if ($record =~ /^(\d+)\./ && -d $local_path) {
            my $taxon_id = $1;
            
            ## create the remote FTP path
            my $remote_path = "$remote_ftp_root/$taxon_id/$subdir";
            my $destination_path = "$ftp_root/$taxon_id/$subdir";
            $ssh->cmd("mkdir -p $destination_path");

            ## mirror to FTP site
            foreach my $glob(@{$ftp_globs->{$subdir}}) {
                foreach my $file(glob("$local_path/$glob")) {
                    print STDERR "Copying '$file' -> '$remote_path'\n";
                    $scp->scp("$file", "$remote_path/") or warn $scp->{errstr};
                }
            }
            
            ## create the remote BLAST db path
            my $db_remote_path = "$remote_blast_root/$subdir";
            my $db_destination_path = "$blast_root/$subdir";
            $ssh->cmd("mkdir -p $db_destination_path");
            
            ## mirror to BLAST db dir
            foreach my $glob(@{$blast_db_globs->{$subdir}}) {
                foreach my $file(glob("$local_path/$glob")) {
                    print STDERR "Copying '$file' -> '$db_remote_path'\n";
                    $scp->scp("$file", "$db_remote_path/") or warn $scp->{errstr};
                }
            }
            
        }
    }
}

## compress the remote files if the option has been specified
if ($compress) {
    print STDERR "Compressing remote files...\n";
    $ssh->cmd('find '.$ftp_root.' -name "*" -exec bzip2 {} \;');
    ## this is a hack because I don't know how to make find not match the .version files right now
    $ssh->cmd('find '.$ftp_root.' -name ".version.bz2" -exec bunzip2 {} \;');
}

## open read permissions on the FTP site files
print STDERR "Adding read permissions to FTP data...\n";
$ssh->cmd("chmod -R a+rX $ftp_root/");

## open read permisions on the BLAST databases
print STDERR "Adding read permissions to BLAST databases...\n";
$ssh->cmd("chmod -R a+rX $blast_root/");
