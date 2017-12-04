#!/usr/bin/perl

=head1  NAME

do_pangenome_analysis.pl - Performs the merging and analysis of BLAST results to create pangenome data.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input_list=/path/to/somefile.dat.list
        --output=/path/to/output.list
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_list,-i>
    List of files containing serialized array of BLAST data.

B<--output_path,-o>
    Path to which output files will be written.

B<--output_file,-o>
    Full path and name of pangenome data output file.

B<--multiplicity,-m>
    The multiplicity value to be used in data sampling.

B<--comparisons,-c>
    The number of comparisons to sample.

B<--log,-d>
    optional. Will create a log file with summaries of all actions performed.

B<--debug>
    optional. Will display verbose status messages to STDERR if logging is disabled.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    The pangenome analysis script creates an array of BLAST results data which is then
    processed to create pangenome data.

=head1 INPUT

    The input should be a list of files containing serialized BLAST results array data.

=head1 OUTPUT

    There is no output unless you use the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use Pod::Usage;
use File::Basename;
use XML::Twig;
use Math::Combinatorics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Storable qw(nstore retrieve);
#use lib "/export/lib/5.8.8"; ## temp
use Data::Random qw(:all);   ## temp
use Benchmark;
#use DBM::Deep;

#use warnings;
use strict;

my @results = ();
my $dups = {};
my %options = ();
my $results = GetOptions (  \%options,
                            'input_file|i=s',
                            'db_list|dl:s',
                            'blast_stored_file|b=s',
                            'output_path|o=s',
                            'output_file|f=s',
                            'write_lists:i',
                            'multiplicity|m:s',
                            'comparisons|c:s',
                            'debug|d=s',
                            'command_id=s',       ## passed by workflow
                            'logconf=s',          ## passed by workflow (not used)
                            'log|l=s',
                            'help|h',
                         ) || pod2usage();


my $comparisons;
my $multiplicity;
my $db_filter = undef;

if($options{'db_list'}) {
    &read_db_list();
}

unless ($options{'comparisons'} || $options{'multiplicity'}) {
#    die "must provide flag --comparisons or --multiplicity";
}

my $output_path = $options{'output_path'};
$output_path =~ s/\/$//;
if ($output_path eq '') {
    $output_path = '.';
}

my $output_file;
if ($options{'output_file'}) {
    $output_file = $options{'output_file'};
} else {
    $output_file = 'pangenome.table.txt';
}


##open (IN, $options{'input_list'}) || die "couldn't open input list";


if (!-e $options{'input_file'}) {
    die "must specify an input file which exists";
}

open (IN, $options{'input_file'}) || die "couldn't open input file: $!";


my $profiles = {};
while (<IN>) {
    chomp;

    ## skip comment lines
    if (/^#/) { next; }

    my ($db, $gene, @profile) = split("\t", $_);

    $profiles->{$db}->{$gene} = join('', @profile); ##$profile;
}
my @genomes = sort {$a cmp $b} keys(%{$profiles});

my $genes={};

print STDERR "Processing results...\n";

foreach my $genome(@genomes) {
    foreach my $gene(keys(%{$profiles->{$genome}})) {
        for (my $i = 0; $i < scalar(@genomes); $i++) {
            my $flag = substr($profiles->{$genome}->{$gene}, $i, 1);
            if ($flag) {
                $genes->{$genome}->{$gene}->{$genomes[$i]} = 1
            }
        }
    }
}

open(RESULT, ">".$output_path."/".$output_file) || die "couldn't open $output_file for writing";

my %genome_index = ();
my $i_counter = 1;
foreach (@genomes) {
    $genome_index{$_} = $i_counter;
    #print RESULT "## $i_counter\t$_\n";
    $i_counter++;
}

# output the gene counts for single genomes
foreach (@genomes) {
    my $gene_count = scalar keys(%{$genes->{$_}});
    ## our output table will take the form: 
    ## genomecount\tcore\tshared\tnew\tcore dup count\tshared dup\tnew dup
    print RESULT "1\t0\t0\t$gene_count\t0\t0\t0\t".$genome_index{$_}."\t$_\n";
}
print STDERR "done.\n";

if($options{'comparisons'} || $options{'multiplicity'}) {
    &do_analysis_with_sampling();
}
else {
    &do_analysis_no_sampling();
}

exit(0);

sub do_analysis_no_sampling {
my $comp_counter = {};
my $max = scalar @genomes;
for (my $i = 1; $i <= $max; $i++) {
    my $combinat = Math::Combinatorics->new(
                                            count => $i,
                                            data  => [@genomes],
                                           );

#    my $cnt = 1;
    while(my @reference = $combinat->next_combination){
        my $ref_string = "(".join(",",@reference).")\n";
        my @comparison_set = @{array_diff(\@genomes, \@reference)};
        foreach my $comp_genome(@comparison_set) {
            my $dup_counts = {};
            $comp_counter->{$comp_genome}->{$i}++;
            my $genes_by_category = {};
            GENE: foreach my $gene(keys(%{$genes->{$comp_genome}})) {
                my $count = 0;

                foreach my $ref_genome(@reference) {
                    ## if the hash value == 1
                    if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
                        $count++;
                        ## we have a hit
                        $genes_by_category->{$comp_genome}->{'shared'}->{$gene}=1;
                    }
                }       
                ## if the number of genomes we have hits against is equal to the total
                if ($count == scalar @reference) {
                    ## then this is a core gene
                    $genes_by_category->{$comp_genome}->{'core'}->{$gene}=1;
                }
                ## if we didn't have any hits at all
                if ($count == 0) {
                    ## then this is a new (strain specific) gene
                    $genes_by_category->{$comp_genome}->{'new'}->{$gene}=1;
                }
            }
            my $rgcount = (scalar @reference) + 1;
            
            ## process lists to see how many duplicated genes are in each category
            foreach my $cat(('core', 'shared', 'new')) {
                $dup_counts->{$comp_genome}->{$cat} = 0;
                ## we'll look at each set of duplicates one at a time
                foreach my $dup_set(keys(%{$dups->{$comp_genome}})) {
                    my @dup_genes = split(" ", $dup_set);
                    my $dup_count = scalar(@dup_genes);
                    my $dup_counter=0;
                    ## for each gene of a duplicate set
                    foreach my $dup(@dup_genes) {
                        ## check if it's in the category
                        if ($genes_by_category->{$comp_genome}->{$cat}->{$dup}) {
                            $dup_counter++;
                        }
                    }
                    ## check if all of the dups of a set weren't found in the same category
                    if ($dup_counter > 0 && $dup_count != $dup_counter) {
                        print STDERR "Only $dup_counter of the following dup set found:\n$dup_set\n***This could be a problem\n";
                    ## and if they are, then add the dup overcount to the total
                    } elsif ($dup_counter == $dup_count) {
                        $dup_counts->{$comp_genome}->{$cat} += ($dup_count - 1);
                    }
                }
            }
            
            my $core_count     = scalar(keys(%{$genes_by_category->{$comp_genome}->{'core'}}));
            my $shared_count   = scalar(keys(%{$genes_by_category->{$comp_genome}->{'shared'}}));
            my $new_count      = scalar(keys(%{$genes_by_category->{$comp_genome}->{'new'}}));
            my $core_dup_count   = $dup_counts->{$comp_genome}->{'core'};
            my $shared_dup_count = $dup_counts->{$comp_genome}->{'shared'};
            my $new_dup_count    = $dup_counts->{$comp_genome}->{'new'};
            
            print RESULT "$rgcount\t$core_count\t$shared_count\t$new_count\t$core_dup_count\t$shared_dup_count\t$new_dup_count\t$genome_index{$comp_genome}\t$comp_genome\n";

            if ($options{'write_lists'}) { 
            
                open (OUT, ">".$output_path."/".$comp_genome."_core_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
                print OUT "#".$comp_genome."\n";
                print OUT "#".join(" ", @reference)."\n";
                foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'core'}})) {
                    print OUT $g."\n";
                }
                close OUT;
                open (OUT, ">".$output_path."/".$comp_genome."_shared_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
                print OUT "#".$comp_genome."\n";
                print OUT "#".join(" ", @reference)."\n";
                foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'shared'}})) {
                    print OUT $g."\n";
                }
                close OUT;
                open (OUT, ">".$output_path."/".$comp_genome."_new_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
                print OUT "#".$comp_genome."\n";
                print OUT "#".join(" ", @reference)."\n";
                foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'new'}})) {
                    print OUT $g."\n";
                }
                close OUT;
            }
##            $cnt++;                
        }

    }
}
}

sub do_analysis_with_sampling {
if ($options{'comparisons'}) {
    my ($est_multiplicity, $t, $t2) = estimate_multiplicity(scalar(@genomes), $options{'comparisons'});
    print STDERR "# of comparisons: $options{comparisons}\nest. multiplicity: $est_multiplicity\n";
    $options{'multiplicity'} = $est_multiplicity;
} else {
    my ($est_comparisons, $t) =  estimate_comparisons(scalar(@genomes), $options{'multiplicity'}); 
    print STDERR "multiplicity: $options{multiplicity}\nest. # comparisons: $est_comparisons\n";
}

my $comp_counter = {};

my $max = scalar @genomes;

print STDERR "Total of $max genomes in results set\n";
print STDERR "Running analysis...\n";

my %seen;

my @i_genomes = (0 .. ($max - 1));

for (my $i = 1; $i < $max; $i++) {

    # start timer
    my $start = new Benchmark;

    print STDERR "Running analysis with $i reference genomes...";
    
    ## new
    my $true_max = int(factorial($max - 1) / (factorial(($i + 1) - 1) * factorial($max - ($i + 1))));
    print STDERR "Maximum number of computations for $i on $max genomes: $true_max\n";
    
    for (my $j = 0; $j < $max; $j++) {
        
        my $comp_genome = $genomes[$j];
        
        my @ref_genomes = @i_genomes;
        splice @ref_genomes, $j, 1;
         
        for (my $point_count=0; ($point_count < $true_max && $point_count < $options{'multiplicity'});){
            my @reference_set = rand_set( set => \@ref_genomes, size => $i, shuffle => 0 );

            my @seen_vec = undef;

            for (my $ii = 0; $ii < scalar(@reference_set); $ii++) {
                $seen_vec[$reference_set[$ii]] = 1;
            }   

            unless ($seen{$j}{"@seen_vec"}){
            
                my $dup_counts = {};
                $comp_counter->{$comp_genome}->{$i}++;
                my $genes_by_category = {};
                GENE: foreach my $gene(keys(%{$genes->{$comp_genome}})) {
                    my $count = 0;

                    foreach my $i_ref_genome(@reference_set) {
                        my $ref_genome = $genomes[$i_ref_genome];
                        ## if the hash value == 1
                        if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
                            $count++;
                            ## we have a hit
                            $genes_by_category->{$comp_genome}->{'shared'}->{$gene}=1;
                        }
                     }       
                     ## if the number of genomes we have hits against is equal to the total
                     if ($count == scalar @reference_set) {
                         ## then this is a core gene
                         $genes_by_category->{$comp_genome}->{'core'}->{$gene}=1;
                     }
                     ## if we didn't have any hits at all
                     if ($count == 0) {
                         ## then this is a new (strain specific) gene
                         $genes_by_category->{$comp_genome}->{'new'}->{$gene}=1;
                     }
                }
                my $rgcount = (scalar @reference_set) + 1;
            
                ## process lists to see how many duplicated genes are in each category
                foreach my $cat(('core', 'shared', 'new')) {
                    $dup_counts->{$comp_genome}->{$cat} = 0;
                    ## we'll look at each set of duplicates one at a time
                    foreach my $dup_set(keys(%{$dups->{$comp_genome}})) {
                        my @dup_genes = split(" ", $dup_set);
                        my $dup_count = scalar(@dup_genes);
                        my $dup_counter=0;
                        ## for each gene of a duplicate set
                        foreach my $dup(@dup_genes) {
                            ## check if it's in the category
                            if ($genes_by_category->{$comp_genome}->{$cat}->{$dup}) {
                                $dup_counter++;
                            }
                        }
                        ## check if all of the dups of a set weren't found in the same category
                        if ($dup_counter > 0 && $dup_count != $dup_counter) {
                            print STDERR "Only $dup_counter of the following dup set found:\n$dup_set\n***This could be a problem\n";
                            ## and if they are, then add the dup overcount to the total
                        } elsif ($dup_counter == $dup_count) {
                            $dup_counts->{$comp_genome}->{$cat} += ($dup_count - 1);
                        }
                    }
                }
            
                my $core_count     = scalar(keys(%{$genes_by_category->{$comp_genome}->{'core'}}));
                my $shared_count   = scalar(keys(%{$genes_by_category->{$comp_genome}->{'shared'}}));
                my $new_count      = scalar(keys(%{$genes_by_category->{$comp_genome}->{'new'}}));
                my $core_dup_count   = $dup_counts->{$comp_genome}->{'core'};
                my $shared_dup_count = $dup_counts->{$comp_genome}->{'shared'};
                my $new_dup_count    = $dup_counts->{$comp_genome}->{'new'};
         
                print RESULT "$rgcount\t$core_count\t$shared_count\t$new_count\t$core_dup_count\t$shared_dup_count\t$new_dup_count\t".$genome_index{$comp_genome}."\t$comp_genome\n";

                if ($options{'write_lists'}) { 
        
                    open (OUT, ">".$output_path."/".$comp_genome."_core_".$rgcount."_".$comp_counter->{$comp_genome}->{$i}.".ids.txt");
                    print OUT "#".$comp_genome."\n";
                    print OUT "#".join(" ", subset(\@genomes, \@reference_set))."\n";
                    foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'core'}})) {
                        print OUT $g."\n";
                    }
                    close OUT;
                    open (OUT, ">".$output_path."/".$comp_genome."_shared_".$rgcount."_".$comp_counter->{$comp_genome}->{$i}.".ids.txt");
                    print OUT "#".$comp_genome."\n";
                    print OUT "#".join(" ", subset(\@genomes, \@reference_set))."\n";
                    foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'shared'}})) {
                        print OUT $g."\n";
                    }
                    close OUT;
                    open (OUT, ">".$output_path."/".$comp_genome."_new_".$rgcount."_".$comp_counter->{$comp_genome}->{$i}.".ids.txt");
                    print OUT "#".$comp_genome."\n";
                    print OUT "#".join(" ", subset(\@genomes, \@reference_set))."\n";
                    foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'new'}})) {
                        print OUT $g."\n";
                    }
                    close OUT;
                }
            
                ## new
                $seen{$j}{"@seen_vec"}=1;
                $point_count++;
        
            } ##end of foreach i_comp_genome
            
        }## end of unless
    }
    
    # end timer
    my $end = new Benchmark;

    # calculate difference
    my $diff = timediff($end, $start);

    print "runtime: ".timestr($diff, 'noc')."\n";
}
}

sub array_diff {
    my ($superset_ref, $subset_ref) = @_;
    my %index;
    my @diff_array;
    foreach (@{$superset_ref}) {
        $index{$_} = 1;
    }
    foreach (@{$subset_ref}) {
        delete $index{$_};
    }
    @diff_array = keys(%index);
    return \@diff_array;
}

no warnings;
sub factorial {
        my $i = shift @_;
        if ($i <= 1) {
                return 1;
        } else {
                return ($i * factorial($i-1));
        }
}
use warnings;

sub subset {
    my ($array_ref, $indices_array_ref) = @_;

    my @result;

    foreach my $i(@{$indices_array_ref}) {
        push (@result, $array_ref->[$i]);
    }

    return @result;
}

sub estimate_multiplicity {
    my ($ngenomes, $req_comparisons) = @_;
    my $i = 0;
    my $lower_mult = 0;
    my $lower_comp = 0;
    my $lower_theor = 0;
    my $upper_mult = 0;
    my $upper_comp = 0;
    my $upper_theor = 0;
    my $ldiff = 0;
    my $udiff = 0;
LOOP:   for ($i=5; $i<=5000; $i++){
        my ($a, $b) = estimate_comparisons($ngenomes, $i); 
        if ($a < $req_comparisons){
            $lower_mult = $i;
            $lower_comp = $a;
            $lower_theor = $b;
        } else {
            $upper_mult = $i;
            $upper_comp = $a;
            $upper_theor = $b;
            last LOOP;
        }
    }
    if ($upper_mult == 5){
        my ($a,$b) = estimate_comparisons($ngenomes, 5);
        return (5,$a,$b);
    }
    elsif ($lower_mult == 5000){
        return (5000,$lower_comp,$lower_theor);
    }
    else {
        $ldiff = $req_comparisons - $lower_comp;
        $udiff = $upper_comp - $req_comparisons;
        if ($ldiff <= $udiff){
            return ($lower_mult, $lower_comp, $lower_theor);
        } else {
            return ($upper_mult, $upper_comp, $upper_theor);
        }
    }   
}

sub estimate_comparisons{
    my ($ngenomes, $multiplex) = @_;
    
    my $i = 0;
    my $tot_comparisons = 0;
    my $theor_comparisons = 0;
    my $theor = 0;
    my $real = 0;
    for ($i=2; $i<=$ngenomes; $i++){
        my $theor = factorial($ngenomes) / ( factorial($i - 1) * factorial($ngenomes - $i) );
        my $real = $multiplex * $ngenomes;
        $theor_comparisons += $theor;
        if ($theor < $real){
            $tot_comparisons += $theor;
        } else {
            $tot_comparisons += $real;
        }
    }
    
    return ($tot_comparisons, $theor_comparisons);
}

sub read_db_list {
    open IN, "<".$options{'db_list'} || die "couldn't open '$options{db_list}' for reading: $!";
    while(<IN>) {
        chomp;
        $db_filter->{$_} = 1;
    }close IN;
    
}
