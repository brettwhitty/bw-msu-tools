#!/usr/bin/perl

# parse_orthomcl_results.pl

# 3 August 2009
# Kevin Childs

# This script will parse an orthomcl output file and will output multifasta files 
# for each cluster.  Because of the large numbers of files that may be created, 
# a directory tree will be used to store the files.
# Additionally, the cluster fastas will be renamed and a decoder
# file will be written with the original names of the fastas.

use Getopt::Std;
use Bio::SeqIO;
use Digest::MD5 qw(md5_base64);
use File::Path;

use strict;

my $usage = "\n$0 -i orthomcl_results_file  [-f file_with_paths_to_fasta_files || -F multifasta_file]  -d path_to_root_working_dir  -P depth_multifasta_directory_tree  -p prefix_for_multifasta_files  -I index_start_number  [-D cluster_division_factor]  -o output_decoder_file\n\n";

our ( $opt_i, $opt_o, $opt_p, $opt_P, $opt_D, $opt_d, $opt_I, $opt_f, $opt_F, $opt_h );
getopts("i:o:p:P:D:d:f:F:I:h") or die usage();

if ($opt_h) {
    print $usage;
    exit;
}

my $input_file = $opt_i;
my $output_directory = $opt_d;
my $output_prefix = $opt_p;
my $decoder_file = $opt_o;
my $fasta_path_file = $opt_f;
my $multifasta_file = $opt_F;
my $index_start_number = $opt_I;
my $path_depth = $opt_P;
my $division_factor = $opt_D;

if (   !defined($input_file)
       || !defined($output_directory)
       || !( -e $input_file )
       || !( -d $output_directory )
       || !( -w $output_directory ) 
       || !defined($output_prefix)
       || !defined($path_depth)
       || $path_depth =~ /\D/
       || !defined($decoder_file)
       || (!defined($fasta_path_file) && !defined($multifasta_file)) 
       || !defined($index_start_number)) {
    die "\nB Missing or invalid input values.\n" . $usage;
}

if (!defined($division_factor)) {
    $division_factor = 101;
}

# Prepare $output_directory value for use.
if ( $output_directory !~ /\/$/ ) {
    $output_directory .= "/";
}

my %seqs;

# Read in all sequences that could be needed here.
if (defined($fasta_path_file)) {
    open IN, $fasta_path_file || die "\nUnable to open $fasta_path_file for reading.\n\n";
    while (my $line = <IN>) {
        chomp $line;
        my $multifasta = Bio::SeqIO->new( -file => $line, -format => "fasta" );

        while ( my $seq_obj = $multifasta->next_seq() ) {
            if (exists($seqs{$seq_obj->display_id()})) {
                die "The sequence, " . $seq_obj->display_id() . ", already exists\n";
            }
            $seqs{$seq_obj->display_id()} = $seq_obj->seq();
        }
    }
}
else {
    my $multifasta = Bio::SeqIO->new( -file => $multifasta_file, -format => "fasta" );

    while ( my $seq_obj = $multifasta->next_seq() ) {
        if (exists($seqs{$seq_obj->display_id()})) {
            die "The sequence, " . $seq_obj->display_id() . ", already exists\n";
        }
        $seqs{$seq_obj->display_id()} = $seq_obj->seq();
    }
}

# Read orthomcl results file.  For each cluster, rename each sequence, write the
# new name/old name pair to the decoder file, write the sequence with its new name
# to the output file.
open IN, $input_file || die "\nUnable to open input file, $input_file, for reading.\n\n";

my $cluster_number = $index_start_number;
my $seqs_in_cluster = 0;
my $index = 0;
my $output_file;
my @decoder_data;
# Print each individual fasta into its own file.
while ( my $line = <IN> ) {
    chomp $line;

    if ($line !~ /^ORTHOMCL/) {
        next;
    }

    $seqs_in_cluster = 0;


    ### HACKING FOR HAINING'S PARSED OUTPUT

    # Strip out the group information from $line.
    #$line =~ s/^ORTHOMCL\d+\(\d+ genes,\d+ taxa\):\s+//;

    my @t = split("\t", $line);
    my @elems = split " ", $t[4];

    # Get the name of each sequence from the @elems array.
    # Recode the names and print the sequences in $multifasta_out.
    # Print the decoding information in DECODER.
    my $seq_name;
    my $fasta_filename = $output_prefix . "_$cluster_number" . ".fasta";
    my $cluster_name = "#" . $output_prefix . "_$cluster_number";
    push  @{$decoder_data[$index]}, $cluster_name;
    my $coded_output_directory = $output_directory . generate_file_path($fasta_filename, $path_depth);
    if (!(-d $coded_output_directory)) {
        mkpath($coded_output_directory);
    }
    my $multifasta_out = $coded_output_directory . $fasta_filename;
    my $out = Bio::SeqIO->new(-file => ">$multifasta_out" , '-format' => 'fasta');
    foreach my $elem (@elems) {

#        if ($elem =~ /(.+)\(/) {
#            $seq_name = $1;
#        }

        $seq_name = $elem;

        if (!exists($seqs{$seq_name})) {
            die "\nUnable to find $seq_name in the fasta sequences\n\n";
        }
        ++$seqs_in_cluster;
        push @{$decoder_data[$index]}, "$seqs_in_cluster\t$seq_name";

        my $seq_obj = Bio::Seq->new(-display_id => "$seqs_in_cluster",
                                    -seq => $seqs{$seq_name});

        $out->write_seq($seq_obj);

    }
    ++$cluster_number;
    ++$index;
    if ($index == $division_factor) {
        $index = 0;
    }
}
close IN;

# Because the data are in a hash, they should come out in a rather random order.
open DECODER, ">$decoder_file" || die "\nUnable to open decoder file, $decoder_file, for writing.\n\n";
for (my $i = 0; $i < $division_factor; ++$i) {
    if (defined(@{$decoder_data[$i]})) {
        my $cluster_set = join "\n", @{$decoder_data[$i]};
        print DECODER "$cluster_set\n";
    }
}
close DECODER;

exit;

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

