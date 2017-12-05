#!/opt/rocks/bin/perl

##
## This script will populate/update the SGR FTP site directory
## using the current contents of the Sol::SeqDB database
## 
## It is intended to be used with rsync on the existing directory at:
## blast:/data/ftp/pub/sgr
##
## Brett Whitty
## whitty@msu.edu
##

$| = 1;

use strict;
use warnings;

use Carp;
use IO::Compress::Bzip2 qw{ bzip2 $Bzip2Error };
use File::Path qw{ mkpath };
use FileCache;
use File::Find::Rule;

## use Bio::SeqIO for parsing genbank records
use Bio::SeqIO;

## development install of Sol::SeqDB module
use lib '/home/whitty/work/testseqdb';
use Sol::SeqDB;

## term to use for dir holding sequences of type 'contig'
my $contig_dir = 'bacs';

## base_dir is the base of where the directory structure will be created
my $base_dir = '/share/scratch/whitty/ftp';

## to prevent damage from accidental run
my $safety_flag = shift @ARGV || die "Accidentally running this can cause damage";

## for accessing the contents of the Sol SeqDB
my $db = new Sol::SeqDB;

## populate hash of  all taxa that are in the database
my $taxonomy_hash = $db->get_taxonomy_hash();

## get today's date
my $today = $db->get_today_date();

## get a list of GI numbers for sequences stored in the DB as type 'contig'
my @contigs = $db->get_gi_list('contig');

## get a list of GI numbers for sequences stored in the DB as type 'read'
my @reads = $db->get_gi_list('read');

## read information for each BAC in the list from the database
print STDERR "Reading BACs info from database ... ";
my $bacs = {};
foreach my $contig_gi(@contigs) {
    my $atts = $db->get_gi_atts($contig_gi);
    push (@{$bacs->{$atts->{'taxon_id'}}->{$atts->{'replicon'}}}, $contig_gi);
}
print STDERR "done.\n";

## read information for each read in the list from the database
my $reads = {};
print STDERR "Reading reads info from database ... ";
foreach my $read_gi(@reads) {
    my $atts = $db->get_gi_atts($read_gi);
    push (@{$reads->{$atts->{'taxon_id'}}}, $read_gi);
}
print STDERR "done.\n";

print STDERR "Updating BACs...\n";
foreach my $taxon_id(keys %{$bacs}) {
    ## path to store contig sequences
    my $path = $base_dir.'/.by_taxon_id/'.$taxon_id.'/'.$contig_dir;

    ## created the directory if it doesn't exist
    unless (-d $path) {
        mkpath($path);
    }

    ## prepare an organism string for using in filenames
    my $organism = $taxonomy_hash->{$taxon_id}->{'scientific_name'};
    my $org = $organism;
    print STDERR "$organism:\n";
    $org =~ s/ /_/g;

    ## output file for all sequences
    my $all_fsa = "$path/$taxon_id.$org.all.fsa";

    ## output file for all genbank records
    my $all_gbk = "$path/$taxon_id.$org.all.gbk";

    ## file to store the version information
    my $version = "$path/.version";

    open(my $all_fsa_fh, ">".$all_fsa) || die "Failed to open '$all_fsa' for writing: $!";
    open(my $all_gbk_fh, ">".$all_gbk) || die "Failed to open '$all_gbk' for writing: $!";

    foreach my $replicon(sort keys %{$bacs->{$taxon_id}}) {
        ## file to store subset of sequences by replicon
        my $replicon_fsa = "$path/$taxon_id.$org.$replicon.fsa";

        ## file to store subset of genbank records by replicon
        my $replicon_gbk = "$path/$taxon_id.$org.$replicon.gbk";

        open (my $replicon_fsa_fh, ">".$replicon_fsa) || die "Failed to open '$replicon_fsa' for writing: $!";
        open (my $replicon_gbk_fh, ">".$replicon_gbk) || die "Failed to open '$replicon_gbk' for writing: $!";
        print STDERR "Writing $replicon_gbk ... ";
        foreach my $gi(sort {$a <=> $b} @{$bacs->{$taxon_id}->{$replicon}}) {

            ## write genbank
            my $gbk = $db->get_genbank($gi);
            print $all_gbk_fh $gbk;
            print $replicon_gbk_fh $gbk;

            ## write fasta
            my $fsa = $db->get_gb_fasta($gi);
            print $all_fsa_fh $fsa;
            print $replicon_fsa_fh $fsa;

            ## check for RHPOTKEY
            if ($gbk =~ /RHPOTKEY/i) {
                ## write RHPOTKEY subset sequence file by replicon
                my $rhpotkey_fsa = "$path/4113.Solanum_tuberosum.$replicon.RHPOTKEY.fsa";

                ## write RHPOTKEY subset genbank file by replicon
                my $rhpotkey_gbk = "$path/4113.Solanum_tuberosum.$replicon.RHPOTKEY.gbk";

                ## write RHPOTKEY subset sequence file
                my $rhpotkey_all_fsa = "$path/4113.Solanum_tuberosum.all.RHPOTKEY.fsa";

                ## write RHPOTKEY subset genbank file
                my $rhpotkey_all_gbk = "$path/4113.Solanum_tuberosum.all.RHPOTKEY.gbk";

                ## use FileCache for writing the files
                no strict;
                cacheout $rhpotkey_fsa;
                cacheout $rhpotkey_gbk;
                cacheout $rhpotkey_all_fsa;
                cacheout $rhpotkey_all_gbk;
                print $rhpotkey_fsa $fsa;
                print $rhpotkey_gbk $gbk;
                print $rhpotkey_all_fsa $fsa;
                print $rhpotkey_all_gbk $gbk;
                use strict;
            }
        }


        close $replicon_fsa_fh;

        close $replicon_gbk_fh;

        print STDERR "done.\n";
    }
    
    close $all_fsa_fh;

    close $all_gbk_fh;
    
    ## write version file
    open(my $ver_fh, ">$version") || die "Failed to write '$version': $!";
    print $ver_fh "$organism\t$taxon_id\t$today";

    print STDERR "$organism complete.\n";
}
print STDERR "BACs updated.\n";

print STDERR "Updating reads...\n";
foreach my $taxon_id(keys %{$reads}) {
    my $read_base_dir = $base_dir.'/.by_taxon_id/'.$taxon_id;

    ## prepare an organism string for using in filenames
    my $organism = $taxonomy_hash->{$taxon_id}->{'scientific_name'};
    my $org = $organism;
    print STDERR "$organism:\n";
    $org =~ s/ /_/g;

#    open(my $all_fsa_fh, ">".$all_fsa) || die "Failed to open '$all_fsa' for writing: $!";
#    open(my $all_gbk_fh, ">".$all_gbk) || die "Failed to open '$all_gbk' for writing: $!";
    foreach my $gi(sort {$a <=> $b} @{$reads->{$taxon_id}}) {

        ## fetch genbank record
        my $gbk = $db->get_genbank($gi);

        open(my $gbk_sfh, '<', \$gbk) || croak "Failed to open file handle on string";

        ## get information from the gbk record
        my $seqio_obj = Bio::SeqIO->new(-fh => $gbk_sfh);
        my $seq = $seqio_obj->next_seq();
        my $annot = $seq->annotation();
        my ($comment) = $annot->get_Annotations('comment');
        my $comment_text = $comment->{'text'};
        $comment_text =~ /Class: ([A-Za-z]+\s*[a-z]*)/ || carp "Failed to find Class in comment '$comment_text'";
        my $read_class = $1 || 'gss';
        my $read_type = lc($read_class);
        $read_type =~ s/\s/_/g;

        ## do rhpotkey check
        my $desc = $seq->desc();
        my $is_rhpotkey = 0;
        if ($desc =~ /RHPOTKEY/) {
            $is_rhpotkey = 1;
        }

        ## extract record type from genbank record
        my $path = $read_base_dir.'/'.$read_type;

        unless (-d $path) {
            mkpath($path);
        }

        ## fetch fasta record
        my $fsa = $db->get_gb_fasta($gi);

        my $all_fsa = "$path/$taxon_id.$org.fsa";
        my $all_gbk = "$path/$taxon_id.$org.gbk";
###        my $version = "$path/.version";


        ## write fasta and genbank to 'all' files
        no strict;
###        my $version = "$path/.version";
#        my $all_fsa_fh = cacheout $all_fsa;
#        my $all_gbk_fh = cacheout $all_gbk;
        cacheout $all_fsa;
        cacheout $all_gbk;
#        print $all_fsa_fh $fsa;
#        print $all_gbk_fh $gbk;
        print $all_fsa $fsa;
        print $all_gbk $gbk;
        use strict;


        ## check for RHPOTKEY 
        if ($is_rhpotkey) {
            
            my $rhpotkey_fsa = "$path/4113.Solanum_tuberosum.RHPOTKEY.fsa";
            my $rhpotkey_gbk = "$path/4113.Solanum_tuberosum.RHPOTKEY.gbk";
            no strict;
#            my $rh_fsa_fh = cacheout $rhpotkey_fsa;
#            my $rh_gbk_fh = cacheout $rhpotkey_gbk;
            cacheout $rhpotkey_fsa;
            cacheout $rhpotkey_gbk;
#            print $rh_fsa_fh $fsa;
#            print $rh_gbk_fh $gbk;
            print $rhpotkey_fsa $fsa;
            print $rhpotkey_gbk $gbk;
            use strict;
        }


    }
#    print STDERR "Compressing $all_gbk ... ";
    
#    close $all_fsa_fh;
#    my $status = bzip2 $all_fsa => $all_fsa.'.bz2'
#        or die "bzip2 failed: $Bzip2Error";
#    unlink $all_fsa;

#    close $all_gbk_fh;
#    $status = bzip2 $all_gbk => $all_gbk.'.bz2'
#        or die "bzip2 failed: $Bzip2Error";
#    unlink $all_gbk;
    
#    print STDERR "done.\n";

    ## write version file
    #open(my $ver_fh, ">$version") || die "Failed to write '$version': $!";
    #print $ver_fh "$organism\t$taxon_id\t$today";

    print STDERR "$organism complete.\n";
}
print STDERR "Reads updated.\n";

## compress gbk files
print STDERR "Compressing gbk files ... \n";
my @gbk_files = File::Find::Rule->file()->name( '*.gbk' )->extras({ follow => 0 })->in( $base_dir );
foreach my $gbk_file(@gbk_files) {
    print $gbk_file."\n";
    my $status = bzip2 $gbk_file => $gbk_file.'.bz2'
        or die "bzip2 failed: $Bzip2Error";
    unlink $gbk_file;
}
print STDERR "done.\n";

## compress fsa files
print STDERR "Compressing fsa files ... \n";
my @fsa_files = File::Find::Rule->file()->name( '*.fsa' )->extras({ follow => 0 })->in( $base_dir );
foreach my $fsa_file(@fsa_files) {
    print $fsa_file."\n";
    my $status = bzip2 $fsa_file => $fsa_file.'.bz2'
        or die "bzip2 failed: $Bzip2Error";
    unlink $fsa_file;
}
print STDERR "done.\n";

print STDERR "All done.\n";
