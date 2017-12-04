#!/opt/rocks/bin/perl

# Created August 2008

# Three input files are needed.  One is the fasta file with the sequences, one is the 
# gff3 output file from the ssr_finder, and the other is the output file name
#
# Original author Morgan Chaires
# Bug fixes and refactoring by Brett Whitty

use warnings;
use strict;

use lib '/home/whitty/SVN/lib/';
use GFFTextEncoder;
use Bio::DB::Fasta;
use Data::Dumper;
use File::Temp qw{ tempfile };
use Getopt::Long;
use File::Basename;
use Carp;
use Cwd qw{ abs_path };

my ($fasta, $gff, $output_dir, $out_gff, $out_txt, $inplace);

my $usage = "\n$0   -f fasta_file -g ssr_gff -o output_dir\n\n";


GetOptions(
    'gff|g=s'           =>  \$gff,
    'fasta|f=s'         =>  \$fasta,
    'output_dir|o=s'    =>  \$output_dir,
    'out_gff=s'         =>  \$out_gff,
    'out_txt=s'         =>  \$out_txt,
    'inplace|i!'        =>  \$inplace,
);

if (! defined($gff)) {
    carp "WARNING: You haven't provided a GFF file name, so one will be automatically provided!";
    my $fullpath = abs_path($fasta);
    my $dir = dirname($fullpath);
    my $base = basename($fullpath);
    $base =~ /^(.*)(\.fsa|\.fna|\.fa|\.fasta)/ || croak "Unexpected fasta file suffix, GFF file name magic failed!";
    my $prefix = $1;

    $gff = "$dir/$prefix.ssr.gff";
    if (! -f $gff) {
        croak "Best guess at GFF file name '$gff' does not exist, please provide it with --gff flag";
    }
}

## provide the option to write the output to where the fasta files are
if ($inplace) {
    my $fullpath = abs_path($fasta);
    my $dir = dirname($fullpath);
    $output_dir = $dir;
    carp "WARNING: You have used the --inplace flag, output will be written to '$output_dir'"; 
}

## output dir must exist
if (! defined($output_dir) || ! -d $output_dir) {
    croak "You must provide an output directory with the --output_dir flag, or use the --inplace flag";
}

## autonaming in effect
if (! defined($out_gff) || ! defined($out_txt)) {
    if ($gff !~ /\.ssr\.gff$/) {
        croak "If you are not specifying output file names with --out_gff and --out_txt you must have an input SSR GFF filename that matches '*.ssr.gff'";
    } else {
        my $base = basename($gff, '.ssr.gff');
        $out_gff = "$output_dir/$base.ssr.primer3.gff";
        $out_txt = "$output_dir/$base.ssr.primer3.txt";
    }
} else {
    ## strip any paths that the user might have provided with the output file names
    ## and use the required output dir
    $out_gff = $output_dir.'/'.basename($out_gff);
    $out_txt = $output_dir.'/'.basename($out_txt);
}

my $flanking_sequence_size = 25;
my $ssr_ds = {};

# handle for the fasta file
my $db = Bio::DB::Fasta->new($fasta);

## strip any paths that the user might have provided with the output file names
## and use the required output dir
# first step for primer3
my $primer3_input = create_primer_global_record();

open (my $infh, '<', $gff) or croak "$!";

carp "Parsing SSR GFF...";
# Goes through each SSR in a file
while (my $line = <$infh>) {

    # gets SSR information
    my @ssr_info = split ('\t', $line);  
    my ($ssr_id, $motif, $num_repeats) = ($ssr_info[8] =~ /ID=([\w|\-|\d|\.]+);Name=\((\w+)\)(\d+)/);
    my $motif_length = length($motif);
    my $ssr_length = $motif_length * $num_repeats;
   
    # retrieves sequence for current SSR
    my $obj = $db->get_Seq_by_id($ssr_info[0]);
    my $seq = $obj->seq;
    my $seq_length = $obj->length;
        
    # if the sequence starts before the flanking sequence size, or is shorter than it, skip this SSR
    if ( ( $ssr_info[3] < $flanking_sequence_size ) 
         || ( $obj->length < $flanking_sequence_size ) ) {
        next;
    } 
    else {
        $ssr_ds->{$ssr_id} = {};
        $ssr_ds->{$ssr_id}->{id} = $ssr_info[0];
        $ssr_ds->{$ssr_id}->{seq_len} = $seq_length;
        $ssr_ds->{$ssr_id}->{seq} = $obj->seq;
        $ssr_ds->{$ssr_id}->{motif} = $motif;
        $ssr_ds->{$ssr_id}->{motif_len} = $motif_length;
        $ssr_ds->{$ssr_id}->{ssr_id} = $ssr_id; 
        $ssr_ds->{$ssr_id}->{num_repeats} = $num_repeats; 
        $ssr_ds->{$ssr_id}->{ssr_start} = $ssr_info[3];
        $ssr_ds->{$ssr_id}->{ssr_end} = $ssr_info[4];

        $primer3_input = create_primer3_seq_record($primer3_input, $seq, $ssr_length,  $ssr_info[3], $ssr_id);  
    }
}
close $infh;

my (undef, $primer3_in_tmp) = tempfile(OPEN => 0);
my (undef, $primer3_out_tmp) = tempfile(OPEN => 0);

open (my $primer3_in_fh, '>', $primer3_in_tmp) || die "$!";
print $primer3_in_fh $primer3_input;
close ($primer3_in_fh) if $primer3_in_fh;

# run primer3
carp "Running primer3...";
system("/share/apps/primer3-1.1.4/bin/primer3 <$primer3_in_tmp >$primer3_out_tmp");

carp "Generating GFF...";
generate_gff($primer3_out_tmp, $out_gff);
carp "Generating text format table...";
generate_text($primer3_out_tmp, $out_txt);
carp "Done.";

#############################################################################

sub generate_gff {

    my ($primer3_out_tmp, $out_gff) = @_;

    open (my $primer3_out_fh, '<', $primer3_out_tmp) || die "$!";
    
    my $primers = parse_primer3($primer3_out_fh);
    close ($primer3_out_fh) if $primer3_in_fh;

    my %primer_seq;

    #print Dumper $primers;
    open (my $outfh, '>', $out_gff) || die "$!";

    for my $ssr_id(sort keys %{$ssr_ds}) {
    
        # if primers exist, prints primer information
        if (defined %{$primers->{$ssr_id}}) {
    
            for my $primer(keys %{$primers->{$ssr_id}}) {

                # current primer number .... primer_1
                $primer =~ /\w+(\d+)/; 
        
                #prints PCR product to the gff3 file
                print $outfh $ssr_ds->{$ssr_id}->{id}."\tprimer3\tPCR_product\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_start}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_start}."\t.\t.\t.\t";
                print $outfh "ID=".$ssr_id.".pcr_product_$1";
                print $outfh ";Note=".encode(
                    'F: GC='. $primers->{$ssr_id}->{$primer}->{primer_left_GC}.'%'
                     .' Tm='. $primers->{$ssr_id}->{$primer}->{primer_left_TM}.'C'
                  .' R: GC='. $primers->{$ssr_id}->{$primer}->{primer_right_GC}.'%'
                     .' Tm='. $primers->{$ssr_id}->{$primer}->{primer_right_TM}.'C'
                     .' / product='.$primers->{$ssr_id}->{$primer}->{product_size}.'bp'
                );  
                print $outfh ";derives_from=$ssr_id\n";
       
                # prints left primer to the gff3 file
                my $primer_left_end = ($primers->{$ssr_id}->{$primer}->{primer_left_start} + $primers->{$ssr_id}->{$primer}->{primer_left_length}) - 1;
        
                print $outfh $ssr_ds->{$ssr_id}->{id}."\tprimer3\tprimer\t";

                if ($primers->{$ssr_id}->{$primer}->{primer_left_start} < $primer_left_end) {
                    print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_start}."\t";
                    print $outfh $primer_left_end."\t.\t+\t.\t";
                } else {          
                    print $outfh $primer_left_end."\t";
                    print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_start}."\t.\t+\t.\t";
                }

                print $outfh "ID=".$ssr_id.".primer_$1f;";
                print $outfh "Name=".$ssr_id.".primer_$1f;";
                print $outfh "Parent=".$ssr_id.".pcr_product_$1;";
                print $outfh "tm=".$primers->{$ssr_id}->{$primer}->{primer_left_TM}.";";
                print $outfh "gc=".$primers->{$ssr_id}->{$primer}->{primer_left_GC}."\n";
  
                # stores primer ID and sequence to be printed at the end of the gff file
                $primer_seq{$ssr_id.'.primer_'.$1.'f'} = $primers->{$ssr_id}->{$primer}->{primer_left_seq};
        
                # prints right primer to the gff3 file
                my $primer_right_end = ($primers->{$ssr_id}->{$primer}->{primer_right_start} - $primers->{$ssr_id}->{$primer}->{primer_right_length}) + 1;
        
                print $outfh $ssr_ds->{$ssr_id}->{id}."\tprimer3\tprimer\t";
      
                if ($primers->{$ssr_id}->{$primer}->{primer_right_start} <  $primer_right_end ) {
                    print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_start}."\t";
                    print $outfh $primer_right_end."\t.\t-\t.\t";    
                } else {
                    print $outfh $primer_right_end."\t";
                    print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_start} ."\t.\t-\t.\t";
                }
        
                print $outfh "ID=".$ssr_id.".primer_$1r;";
                print $outfh "Name=".$ssr_id.".primer_$1r;";
                print $outfh "Parent=".$ssr_id.".pcr_product_$1;";
                print $outfh "tm=".$primers->{$ssr_id}->{$primer}->{primer_right_TM}.";";
                print $outfh "gc=".$primers->{$ssr_id}->{$primer}->{primer_right_GC}."\n";

                # stores primer ID and sequence to be printed at the end of the gff file
                $primer_seq{$ssr_id.'.primer_'.$1.'r'} = $primers->{$ssr_id}->{$primer}->{primer_right_seq};
            }
        }
    }

    if (keys %primer_seq) {
        print $outfh "##FASTA\n";

        for my $primer(keys %primer_seq) {
            print $outfh ">$primer\n".$primer_seq{$primer}."\n";    
        }
    }

}


sub generate_text {

    my ($primer3_out_tmp, $out_txt) = @_;

    open (my $primer3_out_fh, '<', $primer3_out_tmp) || die "$!";
    my $primers = parse_primer3($primer3_out_fh);
    close ($primer3_out_fh) if $primer3_in_fh;

    # Print primers with ssr's to tab delimited file
    #print Dumper $primers;
    open (my $outfh, '>', $out_txt) || die "$!";

    # prints header to file
    print $outfh "GB_ACC\tSSR_ID\tSEQ_LEN\tSSR_START\tSSR_END\t(MOTIF)#REPEATS\t";
    print $outfh "L_PRIMER_1\tR_PRIMER_1\tL_START_1\tL_LEN_1\tL_TM_1\tL_%GC_1\tR_START_1\tR_LEN_1\tR_TM_1\tR_%GC_1\tPRODUCT_SIZE_1\t";
    print $outfh "L_PRIMER_2\tR_PRIMER_2\tL_START_2\tL_LEN_2\tL_TM_2\tL_%GC_2\tR_START_2\tR_LEN_2\tR_TM_2\tR_%GC_2\tPRODUCT_SIZE_2\n";

    for my $ssr_id(sort keys %{$ssr_ds}) {

        # prints SSR information
        print $outfh $ssr_ds->{$ssr_id}->{id}."\t";
        print $outfh $ssr_id."\t";
        print $outfh $ssr_ds->{$ssr_id}->{seq_len}."\t";
        print $outfh $ssr_ds->{$ssr_id}->{ssr_start}."\t";
        print $outfh $ssr_ds->{$ssr_id}->{ssr_end}."\t";
        print $outfh "(".$ssr_ds->{$ssr_id}->{motif}.")".$ssr_ds->{$ssr_id}->{num_repeats}."\t";
        
        # if primers exist, prints primer information
        if (defined %{$primers->{$ssr_id}}) {
            for my $primer(keys %{$primers->{$ssr_id}}) {

                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_seq}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_seq}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_start}."\t";      
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_length}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_TM}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_left_GC}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_start}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_length}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_TM}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{primer_right_GC}."\t";
                print $outfh $primers->{$ssr_id}->{$primer}->{product_size}."\t";
            }
        }   
        print $outfh "\n";    
    }
}


sub create_primer_global_record{
    my $primer3_global_input = '';
    
    #set up the intial 'Global' primer3 settings
    my $primer3_settings = get_primer3_settings();

    foreach my $key(keys(%{$primer3_settings})) {
        my $value = $primer3_settings->{$key};
        $primer3_global_input .= join("=", ($key, $value))."\n";
    }

    return $primer3_global_input;
}

sub get_primer3_settings{
    my $primer3_settings = {
                            PRIMER_EXPLAIN_FLAG => 0,
                            PRIMER_PRODUCT_SIZE_RANGE => "100-150 150-200 200-250 250-300 350-400 400-450 450-500 500-550 550-600 700-800 800-1000",
                            PRIMER_NUM_RETURN => 3,
                            PRIMER_MAX_END_STABILITY => 9.0,
                            PRIMER_MAX_MISPRIMING => 12.00,
                            PRIMER_PAIR_MAX_MISPRIMING => 24.00,
                            PRIMER_MIN_SIZE => 18,
                            PRIMER_OPT_SIZE => 20,
                            PRIMER_MAX_SIZE => 27,
                            PRIMER_MIN_TM => 57.0,
                            PRIMER_OPT_TM => 60.0,
                            PRIMER_MAX_TM => 63.0,
                            PRIMER_MAX_DIFF_TM => 100.0,
                            PRIMER_MIN_GC => 20.0,
                            PRIMER_MAX_GC => 80.0,
                            PRIMER_SELF_ANY => 8.00,
                            PRIMER_SELF_END => 3.00,
                            PRIMER_NUM_NS_ACCEPTED => 0,
                            PRIMER_MAX_POLY_X => 5,
                            PRIMER_OUTSIDE_PENALTY => 0,
                            PRIMER_GC_CLAMP => 0,
                            PRIMER_SALT_CONC => 50.0,
                            PRIMER_DNA_CONC => 50.0,
                            PRIMER_LIBERAL_BASE => 1,
                            PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS => 0,
                            PRIMER_MIN_QUALITY => 0,
                            PRIMER_MIN_END_QUALITY => 0,
                            PRIMER_QUALITY_RANGE_MIN => 0,
                            PRIMER_QUALITY_RANGE_MAX => 100,
                            PRIMER_WT_TM_LT => 1.0,
                            PRIMER_WT_TM_GT => 1.0,
                            PRIMER_WT_SIZE_LT => 1.0,
                            PRIMER_WT_SIZE_GT => 1.0,
                            PRIMER_WT_GC_PERCENT_LT => 0.0,
                            PRIMER_WT_GC_PERCENT_GT => 0.0,
                            PRIMER_WT_COMPL_ANY => 0.0,
                            PRIMER_WT_COMPL_END => 0.0,
                            PRIMER_WT_NUM_NS => 0.0,
                            PRIMER_WT_REP_SIM => 0.0,
                            PRIMER_WT_SEQ_QUAL => 0.0,
                            PRIMER_WT_END_QUAL => 0.0,
                            PRIMER_WT_POS_PENALTY => 0.0,
                            PRIMER_WT_END_STABILITY => 0.0,
                            PRIMER_PAIR_WT_PRODUCT_SIZE_LT => 0.0,
                            PRIMER_PAIR_WT_PRODUCT_SIZE_GT => 0.0,
                            PRIMER_PAIR_WT_PRODUCT_TM_LT => 0.0,
                            PRIMER_PAIR_WT_PRODUCT_TM_GT => 0.0,
                            PRIMER_PAIR_WT_DIFF_TM => 0.0,
                            PRIMER_PAIR_WT_COMPL_ANY => 0.0,
                            PRIMER_PAIR_WT_COMPL_END => 0.0,
                            PRIMER_PAIR_WT_REP_SIM => 0.0,
                            PRIMER_PAIR_WT_PR_PENALTY => 1.0,
                            PRIMER_PAIR_WT_IO_PENALTY => 0.0,
                            PRIMER_IO_WT_TM_LT => 1.0,
                            PRIMER_IO_WT_TM_GT => 1.0,
                            PRIMER_IO_WT_SIZE_LT => 1.0,
                            PRIMER_IO_WT_SIZE_GT => 1.0,
                            PRIMER_IO_WT_GC_PERCENT_LT => 0.0,
                            PRIMER_IO_WT_GC_PERCENT_GT => 0.0,
                            PRIMER_IO_WT_COMPL_ANY => 0.0,
                            PRIMER_IO_WT_NUM_NS => 0.0,
                            PRIMER_IO_WT_REP_SIM => 0.0,
                            PRIMER_IO_WT_SEQ_QUAL => 0.0,
                            PRIMER_TASK => "pick_pcr_primers",
                            PRIMER_FIRST_BASE_INDEX => 1,
                            PRIMER_PICK_ANYWAY => 0,
                           };
    return $primer3_settings;
}

sub create_primer3_seq_record{
    #this sub prints the "sequence input tags" portion of the primer3 input
    my ($primer3_input, $seq, $ssr_length, $ssr_start, $ssr_id) = @_;
        
    #$primer3_input .= "PRIMER_SEQUENCE_ID=".$seq_obj->display_id."|".$ssr_count."\n";
    $primer3_input .= "PRIMER_SEQUENCE_ID=$ssr_id\n";
    $primer3_input .= "SEQUENCE=".$seq."\n";
    $primer3_input .= "TARGET=$ssr_start,$ssr_length\n"; 
    $primer3_input .= "=\n";
    return $primer3_input;
}


sub parse_primer3 {
    my ( $primer3_out_fh ) = @_;
    seek( $primer3_out_fh, 0, 0 );
    #now parse the primer file
    my $primers = {};
    my %seen_primers = ();
    my $locus;
    while (<$primer3_out_fh>) {
    chomp;
    if (/^=/) {
            $locus = undef;
        } elsif (/PRIMER_SEQUENCE_ID/) {
            (undef, $locus) = split(/=/);
            $primers->{$locus}= {};
    } elsif (/PRIMER_LEFT_SEQUENCE/ || /PRIMER_LEFT_(\d+)_SEQUENCE/) {
            my $primer_name;
            #since seeing a primer seq is guaranteed, populated seen hash here
            if (/PRIMER_LEFT_SEQUENCE/) {
                $primer_name = "primer_1";
                $seen_primers{$locus."_".$primer_name} = [ $locus, $primer_name ];
            } else {
                $primer_name = "primer_".$1;
                $seen_primers{$locus."_".$primer_name} = [ $locus, $primer_name ];
            }
            my (undef, $primer_left_seq) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_left_seq'} = $primer_left_seq;
    } elsif (/PRIMER_RIGHT_SEQUENCE/ ||  /PRIMER_RIGHT_(\d+)_SEQUENCE/) {
            my $primer_name;
            if (/PRIMER_RIGHT_SEQUENCE/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $primer_right_seq) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_right_seq'} = $primer_right_seq;
    } elsif (/^PRIMER_LEFT=/ || /^PRIMER_LEFT_(\d+)=/) {
            my $primer_name;
            if (/^PRIMER_LEFT=/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my ($primer_start, $primer_length) = $_ =~ /=(\d+),(\d+)$/;
            $primers->{$locus}->{$primer_name}->{'primer_left_start'} = $primer_start;
            $primers->{$locus}->{$primer_name}->{'primer_left_length'} = $primer_length;
    } elsif (/^PRIMER_RIGHT=/ || /^PRIMER_RIGHT_(\d+)=/) {
            my $primer_name;
            if (/^PRIMER_RIGHT=/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my ($primer_start, $primer_length) = $_ =~ /=(\d+),(\d+)$/;
            $primers->{$locus}->{$primer_name}->{'primer_right_start'} = $primer_start;
            $primers->{$locus}->{$primer_name}->{'primer_right_length'} = $primer_length;
    } elsif (/PRIMER_LEFT_TM/ || /PRIMER_LEFT_(\d+)_TM/) {
            my $primer_name;
            if (/PRIMER_LEFT_TM/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $primer_left_TM) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_left_TM'} = $primer_left_TM;
    } elsif (/PRIMER_RIGHT_TM/ || /PRIMER_RIGHT_(\d+)_TM/) {
            my $primer_name;
            if (/PRIMER_RIGHT_TM/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $primer_right_TM) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_right_TM'} = $primer_right_TM;
    } elsif (/PRIMER_LEFT_GC_PERCENT/ || /PRIMER_LEFT_(\d+)_GC_PERCENT/) {
            my $primer_name;
            if (/PRIMER_LEFT_GC_PERCENT/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $primer_left_GC) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_left_GC'} = $primer_left_GC;
    } elsif (/PRIMER_RIGHT_GC_PERCENT/ || /PRIMER_RIGHT_(\d+)_GC_PERCENT/) {
            my $primer_name;
            if (/PRIMER_RIGHT_GC_PERCENT/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $primer_right_GC) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'primer_right_GC'} = $primer_right_GC;
    } elsif (/PRIMER_PRODUCT_SIZE=/ || /PRIMER_PRODUCT_SIZE_(\d+)/) {
            my $primer_name;
            if (/PRIMER_PRODUCT_SIZE=/) {
                $primer_name = "primer_1";
            } else {
                $primer_name = "primer_".$1;
            }
            my (undef, $product_size) = split(/=/);
            $primers->{$locus}->{$primer_name}->{'product_size'} = $product_size;
    } else {
            next;
    }
    }
    close($primer3_out_fh);
    return $primers;
}



