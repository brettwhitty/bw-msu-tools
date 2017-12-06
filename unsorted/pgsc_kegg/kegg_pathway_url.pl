#!/opt/rocks/bin/perl

use strict;
use warnings;

use Carp;
use Getopt::Long;
use Set::Scalar;

my ($kegg_mapping, $output, $gene_list, $pathway_prefix, $mode, $cuffdiff);

GetOptions(
    'gene_list|g=s'         =>  \$gene_list,
    'kegg_mapping|k=s'      =>  \$kegg_mapping,
    'output|o=s'            =>  \$output,
    'pathway_prefix|p=s'    =>  \$pathway_prefix,
    'mode|m=i'              =>  \$mode,
    'cuffdiff|c=s'          =>  \$cuffdiff,
);

$mode ||= 1;

if (defined($cuffdiff)) {
    $mode = 3;
} else {
    if ($mode == 3) {
        croak "You have specified --mode=3 which requires a cuffdiff output table as input with --cuffdiff flag";
    }
}
my $valid_mode = {
    1   =>  'text output',
    2   =>  'KEGG pathway URL output',
    3   =>  'heatmap KEGG pathway URL output using log fold change',
    4   =>  'text output with gene and ko columns',
};

if (! defined($valid_mode->{$mode})) {
    croak "Invalid mode specified with --mode flag";
}

## default to the KEGG reference pathway prefix
#$pathway_prefix ||= 'ko'; 

if (! defined($gene_list)) {
    croak "Must provide a gene list with --gene_list";
}
if (! defined($kegg_mapping)) {
    croak "Must provide KEGG mapping table with --kegg_mapping";
}

my $outfh;
if (defined($output)) {
    open $outfh, '>', $output or croak "Failed to open '$output' for writing: $!";
} else {
    $outfh = \*STDOUT;
}

my $infh;

my $genes = {};

## initialize hash of genes of interest
open $infh, '<', $gene_list or croak "Failed to open gene list '$gene_list' for reading: $!";
while (<$infh>) {
    chomp;
    $genes->{$_} = 1;
}

my $logs;
if ($mode == 3) {
    open $infh, '<', $cuffdiff or croak "Failed to open cuffdiff output file '$cuffdiff' for reading: $!";
#    my $max = 0;
    while (<$infh>) {
        chomp;

        if (/^test_id/) {
            next;
        }

        my @t = split(/\t/, $_);
        my ($gene_id, $log_change) = ($t[0], $t[6]);

        $logs->{$gene_id} = $log_change;

#        if (abs($log_change) > $max) {
#            $max = abs($log_change);
#        }
    }
#    ## normalize the log values to make coloring easier
#    foreach my $gene_id(keys(%{$logs})) {
#        $logs->{$gene_id} = $heat_map->{$gene_id} / $max;
#    }

}

## store the pathways
my $pathways = {};
my $gene_to_path = {};
my $title = {};

## parse kegg mapping table
open $infh, '<', $kegg_mapping or croak "Failed to open KEGG mapping file '$kegg_mapping' for reading: $!";
while (<$infh>) {
    chomp;
    
    my @t = split("\t", $_);

    my ($gene_id, $ko_id, $path_id, $path_desc) = @t[0 .. 3];

    if (! defined($genes->{$gene_id})) {
        next;
    }
    
    my $title_id = $path_id;
    $title_id =~ s/[a-z]//ig;
    
    if (! defined($title->{$title_id})) {
        $title->{$title_id} = $path_desc;
    }

    ## if user has specified a different reference pathway, then adjust the prefix
    if ($pathway_prefix) {
        $path_id =~ s/^[a-z]+//;
        $path_id = $pathway_prefix . $path_id;
    }

    ## store the members and corresponding KO in arrays under the pathway ID
    ## there will be duplicates in here for KO
    push(@{$pathways->{$path_id}->{'gene'}}, $gene_id);
    push(@{$pathways->{$path_id}->{'ko'}}, $ko_id);

#    push(@{$gene_to_path->{$gene_id}}, $path_id);
}

foreach my $path_id(sort keys(%{$pathways})) {
    
    my $title_id = $path_id;
    $title_id =~ s/[a-z]//ig;
    
    if ($mode == 1) {
        my $pathway_members = new Set::Scalar(@{$pathways->{$path_id}->{'ko'}});
        print $path_id."\t".join(',', $pathway_members->members())."\n";
    } elsif ($mode == 4) {
        my @genes = @{$pathways->{$path_id}->{'gene'}};
        my @kos = @{$pathways->{$path_id}->{'ko'}};
        
        my $gene_count = scalar(@genes);
        
        print join("\t", (
                $path_id,
                $title->{$title_id},
                join(',', @genes),
                join(',', @kos),
            ))."\n";
    } elsif ($mode == 2) {
        my $pathway_members = new Set::Scalar(@{$pathways->{$path_id}->{'ko'}});
        print join("\t", (
            $path_id,
            $title->{$title_id},
            "http://www.genome.jp/kegg-bin/show_pathway?$path_id/"
            .join('/', (
                $pathway_members->members()
            )),
        ))."\n";
    } elsif ($mode == 3) {
        my @genes = @{$pathways->{$path_id}->{'gene'}};
        my @kos = @{$pathways->{$path_id}->{'ko'}};
        my @logs;
        my @url_kos;

        my $gene_count = scalar(@genes);

        for (my $i =0; $i < $gene_count; $i++) {
            ## skip genes we have no log change for in the cuffdiff file
            if (! defined($logs->{$genes[$i]})) {
                    next;
            }
            my $log_change = $logs->{$genes[$i]};
            
            my $r = 255;
            my $g = 255;
            my $b = 0;

            if ($log_change < 0) {
                $g = 0;
                $r = 128 + 128 * abs($log_change / 30);
            } else {
                $r = 0;
                $g = 128 + 128 * abs($log_change / 30);
            }

            my $rgb = dec2hex($r).dec2hex($g).dec2hex($b);

#            push(@logs, $logs->{$genes[$i]});
            push(@url_kos, $kos[$i]
                #."%09%23FFFFFF,
                ."%09%23$rgb,%23f0f0f0"
#                .":"
#                .$logs->{$genes[$i]}
            );
        }
        print join("\t", (
                $path_id,
                $title->{$title_id},
                "http://www.genome.jp/kegg-bin/show_pathway?$path_id/"
                .join('/', (
#            $pathway_members->members(),
                    'default%3dblack,white',
                    @url_kos
                )),
        ))."\n";
    }
}

sub dec2hex {
    my ($dec) = @_;

    my $hex = unpack("H8", pack("N", $dec));

    return substr($hex, 6, 2);
}
