#! /usr/bin/perl -w

# This script will submit a series of muscle alignment jobs to the local grid.
# All results will be written to a single directory that is indicated by an input
# parameter.  Each qsub job will involve a single multifasta file that contains
# sequences that have already been clustered by orthomcl.  Each qsub job will 
# submit all muscle alignments for a single multifasta file.

# 4 August 2009

# Kevin Childs

use strict;
use Getopt::Std;
use Digest::MD5 qw(md5_base64);

my $usage = "\n$0  -d root_of_working_directories  -c new_file_all_command_line_calls  -s new_file_for_shell_script  -D path_to_decoder_file\n\n";

our ($opt_c, $opt_s, $opt_d, $opt_D, $opt_h);
getopts("c:s:d:D:h") or die "X  $usage";

if (defined($opt_h)) {
    die "$usage";
}

if (!defined($opt_D)) {
    die "opt_D  $usage\n";
}
if (!defined($opt_c)) {
    die "opt_c  $usage";
}
if (!defined($opt_d)) {
    die "opt_d  $usage";
}
if (!defined($opt_s)) {
    die "opt_s  $usage";
}

my $working_dir_root = $opt_d;
my $command_line_file = $opt_c;
my $shell_out = $opt_s;
my $path_to_decoder_file = $opt_D;

if ((-e $command_line_file)) {
    print "A $usage";
    print "$command_line_file\n";
    exit;
}
if (!(-e $working_dir_root)) {
    print "B $usage";
    print "$working_dir_root\n";
    exit;
}
if ((-e $shell_out)) {
    print "D $usage";
    print "$shell_out\n";
    exit;
}
if (!(-e $path_to_decoder_file)) {
    print "D $usage";
    print "$path_to_decoder_file\n";
    exit;
}

my $cwd = `pwd`;
chomp $cwd;
$cwd .= "/";
$command_line_file = $cwd . $command_line_file;

write_shell();

system ("chmod a+x $shell_out");  

open OUT, ">$command_line_file" || die "\nUnable to open $command_line_file for writing.\n\n";
my $num_jobs = 0;
my $count = 0;
my $total_count = 0;

# Get a list of all fasta files in the working directory tree.
# Just regenerate the paths after the root directory.  
# Traversing the tree with File::Find::Rule is too slow.
# But first, we need to divine the $path_depth.
my $path_depth = divine_path_depth($working_dir_root);
## hard coding because it's wrong on rerun
$path_depth = 2;
open IN, $path_to_decoder_file || die "\nUnable to open $path_to_decoder_file for reading.\n\n";
my (@files, %num_seqs, $num_in_cluster);
while (my $line = <IN>) {
    chomp $line;
  NEW_CLUSTER:
    if ($line =~ /#(.+)/) {
        my $fasta_filename = $1 . ".fasta";
        my $multifasta_file = $working_dir_root . "/" . generate_file_path($fasta_filename, $path_depth) . $fasta_filename;
        push @files, $multifasta_file;
        $num_in_cluster = 0;
        while ($line = <IN>) {
            # While we're going through the decoder file, figure out how many seqs in each cluster.
            if ($line !~ /#(.+)/) {
                ++$num_in_cluster;
            }
            else {
                $num_seqs{$multifasta_file} = $num_in_cluster;
                goto NEW_CLUSTER;
            }
        }
        $num_seqs{$multifasta_file} = $num_in_cluster;
    }
}

# Now loop through fasta file array and print command line calls to file.
foreach my $file (@files) {

    ++$count;
    ++$total_count;

    # Prepare a mess of file and directory names.
    $file =~ /(.+)\.fasta$/;
    my $file_minus_suffix = $1;
    my $phylip_interleave_output_file = $file_minus_suffix . ".phyi";
    my $gblocks_phylip_interleave_output_file = $file_minus_suffix . ".phyi.trim";
    my $gblocks_fasta_input_file = $file_minus_suffix . ".gfsa";
    my $gblocks_fasta_output_file = $file_minus_suffix . ".gfsa.trim";
    my $msa_output_file = $file_minus_suffix . ".msa";

    $file =~ /(.+\/)\w+\.fasta/;
    my $path_for_current_file = $1;

    my $temp_dir = $path_for_current_file . "/temp_dir_" . $total_count;  # Use $total_count so that there will be no collisions.
    my $parameter_file = $temp_dir . "/parameter_file.txt";
    my $path_to_outtree = $temp_dir . "/outtree";
    my $path_to_outfile = $temp_dir . "/outfile";
    my $proml_output_file = $file_minus_suffix . ".newick";
    my $phyloxml_output_file = $file_minus_suffix . ".xml";

    print OUT "/share/apps/bin/muscle -in $file -quiet -phyiout $phylip_interleave_output_file -maxmb 4000\n";
    print OUT "java -jar /share/apps/readseq.jar -inform=12 -f=8 -o $gblocks_fasta_input_file $phylip_interleave_output_file\n";
    print OUT "/share/apps/bin/Gblocks $gblocks_fasta_input_file -t=p -e=.trim -p=t -b5=h\n";
    print OUT "java -jar /share/apps/readseq.jar -inform=8 -f=12 -o $gblocks_phylip_interleave_output_file $gblocks_fasta_output_file\n";
    print OUT "mkdir $temp_dir\n";
    print OUT "cd $temp_dir\n";
    print OUT "cp $gblocks_phylip_interleave_output_file infile\n";
    print OUT "echo \"Y\" | cat > $parameter_file\n";  # $parameter_file needed because proml is so 1970's.
    if (!defined($num_seqs{$file})) {
        print "$file\n";
    }
    if ($num_seqs{$file} > 2) {
        # Clusters with more than two seqs are handled correctly by proml.
        print OUT "/opt/bio/phylip/exe/proml < $parameter_file\n";  
        print OUT "#placeholder line for splitting the file\n";
    }
    else {
        # This cluster just has two sequences.  Proml will seg fault with just two sequences.
        # Use protdist to calculate a distance.  Then make newick file for a fake tree.
        print OUT "/opt/bio/phylip/exe/protdist < $parameter_file\n"; 
        # Write a newick file. (A:0.1,B:0.2)
        print OUT "/run/whitty/potato/v3_2/orthomcl/bin/create_newick_file.pl  -i $path_to_outfile  -o $path_to_outtree  \n";
    }
    # Convert the newick file to xml.
    # Keep a copy of the newick tree.    
    # Give the sequences in the xml file proper names and make a few other changes to the xml file.
    # Clean up the junky temp directory.
    print OUT "java -cp /share/apps/forester_dev.jar org.forester.application.phyloxml_converter -f=gn  $path_to_outtree $phyloxml_output_file\n";  
    print OUT "cp $path_to_outtree $proml_output_file\n";
    print OUT "/run/whitty/potato/v3_2/orthomcl/bin/decode_sequence_ids.pl  -i $phyloxml_output_file  -d $path_to_decoder_file\n"; 
    print OUT "/run/whitty/potato/v3_2/orthomcl/bin/convert_and_decode_phyi_file.pl  -i $phylip_interleave_output_file  -o $msa_output_file  -d $path_to_decoder_file\n"; 
    print OUT "cd $cwd\n";
    print OUT "rm -rf $temp_dir\n";  
    if ($count == 100) {
        # The grid jobs will be run in batches of 1000.
        $count = 0;
        ++$num_jobs;
       print OUT "^\n";
    }
}
close OUT;

print "qsub -R y -V  -t 1-$num_jobs -v \"INFILE=$command_line_file\" $shell_out\n";
system("qsub -R y -V  -t 1-$num_jobs -v \"INFILE=$command_line_file\" $shell_out");

exit 0;

sub write_shell {
    open SHELL_OUT, ">$shell_out" || die "\nUnable to open $shell_out for writing.\n\n";
    my $shell_script = <<END_DOC;
#!/bin/bash

#\$ -cwd   -l virtual_free=4G
#\$ -o /dev/null
#\$ -e /dev/null

date

OLD_IFS=IFS
export OLD_IFS

COUNTER=0
while read -d ^ LINE
do

   COUNTER=\$((\$COUNTER+1))
   if [ \$COUNTER == \$SGE_TASK_ID ]; then
       IFS="\n"
       set -- \$LINE
       array=(\$LINE)

       for COMMAND in \${array[\@]}; do
           echo \$COMMAND
           eval \$COMMAND
       done
   fi
done < \$INFILE

date

END_DOC

    print SHELL_OUT $shell_script;
    close SHELL_OUT;

}

sub generate_file_path {

    my ($fasta_filename, $path_depth) = @_;

    if ($path_depth > 22) {
        $path_depth = 22;
    }
    my $md5_sum = md5_base64($fasta_filename);

    my $md5_sum_short = $md5_sum;
    $md5_sum_short =~ s/\W//g;
    $md5_sum_short = lc($md5_sum_short);

    my @elems = split "", $md5_sum_short;

    my $path;
    for (my $i = 0; $i < $path_depth; ++$i) {
        $path .= $elems[$i] . "/";
    }

    return $path;
}

sub divine_path_depth {
    my ($dir) = @_;

    my $depth = 0;

    opendir(DIR, $dir);
    my @records = readdir(DIR);
    closedir(DIR);

    foreach my $record (@records) {
        if ($record eq "." || $record eq "..") {
            next;
        }

        my $full_path_to_record = $dir . "/" . $record;
        if (-d $full_path_to_record) {
            $depth = divine_path_depth($full_path_to_record);
            ++$depth;
            last;
        }       
    }
    return $depth;
}

