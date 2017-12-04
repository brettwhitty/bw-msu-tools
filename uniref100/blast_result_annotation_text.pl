#!/usr/bin/perl

use strict;
use warnings;

## Input is BLAST results from searching against UniRef100
## Returns a table of query accession -> annotation text
##
## Requires BioPerl is in your PERL5LIB or Perl @INC
##
## Brett Whitty
## whitty@msu.edu
##
## usage:
## blast_result_annotation_text.pl --db /state/partition1/db/uniref100.aasubst.fasta --input blast_output.raw
## or
## zcat blast_output.raw.gz | blast_result_annotation_text.pl --db /state/partition1/db/uniref100.aasubst.fasta

use Bio::SearchIO;
use Bio::DB::Fasta;
use Getopt::Long;
use Carp;

my $uniref_db_path;
my $input_file;
my $conserved_unknown_text  = 'Conserved gene of unknown function';
my $unknown_text            = 'Gene of unknown function';
my $add_taxon;
my $reindex;

GetOptions(
    'db|d=s'                =>  \$uniref_db_path,
    'input|i=s'             =>  \$input_file,
    'unknown|u=s'           =>  \$unknown_text,
    'conserved_unknown|c=s' =>  \$conserved_unknown_text,
    'taxon|t!'              =>  \$add_taxon,
    'reindex|r!'            =>  \$reindex,
);

unless (-f $uniref_db_path) {
    croak "Failed to provide path to UniRef100 database fasta file with --db flag";
}

my $infh;

if (defined($input_file) && -f $input_file) {
    open($infh, '<', $input_file) || die "Failed to open '$input_file' for reading: $!";
} else {
    $infh = \*STDIN;
}

my $uniref_db;
if ($reindex) {
    $uniref_db = new Bio::DB::Fasta($uniref_db_path, -reindex => 1);
} else {
    $uniref_db = new Bio::DB::Fasta($uniref_db_path);
}

my $searchio = Bio::SearchIO->new( 
                -format => 'blast',
                -fh =>  $infh,
);

my ($query_name, $subject_name);
while ( my $result = $searchio->next_result() ) {
    my $query_name = $result->query_name;

    my $annotation = [];
    my $hit_counter = 0;
    while( my $hit = $result->next_hit ) {
        my $subject_name = $hit->name;
        $hit_counter++;

        my $fun_ref = get_uniref_function($subject_name);

        if (defined($fun_ref)) {
            push(@{$annotation}, $fun_ref);
        }        
    }
    
    unless (@{$annotation}) {
        if ($hit_counter) {
            $annotation = [
                        [
                            '',
                            $conserved_unknown_text,
                            '',
                            0,
                        ]
            ];
        } else {
            $annotation = [
                        [
                            '',
                            $unknown_text,
                            '',
                            0,
                        ]
            ];
        }
    }

    print join("\t", (
        $query_name,
        (($add_taxon && $annotation->[0]->[0])
            ?
                '{'
                . $annotation->[0]->[0]
                . '|'
                . $annotation->[0]->[2]
                . '} ' 
            : 
                '') ## append species if there is one
        . ucfirst($annotation->[0]->[1]),
    ))."\n";

}

sub get_uniref_function {
    my ($id) = @_;

    my $header = $uniref_db->header($id);

    $header =~ /^\S+\s+(.*)\s+n=(\d+)\s+Tax=(.*)\s+RepID=(\S+)$/;
    my ($desc, $n, $species, $rep_id) = ($1, $2, $3, $4);

    my $weight = 1;

    my $result = undef;

    #### Reject certain annotations outright
    if (_reject_check($desc)) {
        return undef;
    }

    ## Remove meaningless qualifier terms
    if ($desc =~ / \(At\d+g[\d.]+(\/[^)]+)?\)/i) {
        $desc =~ s/ \(At\d+g[\d.]+(\/[^)]+)?\)//i;
    }
    if ($desc =~ / \(At\d+g[\d.]+\/[^)]+\)/i) {
    # Case 1 gets the first AT id, case 1.5 gets the second parenthetical term.
    # F24P17.21 protein  (At3g06330)  (At3g06330/F24P17_21) 
        $desc =~ s/ \(At\d+g[\d.]+\/[^)]+\)//i;
    }
    if ($desc =~ / At(\d+|[MC])g[\d.]+/i) {
        $desc =~ s/ At(\d+|[MC])g[\d.]+//i;
    }
    ## anopheles gambiae locus IDs
    if ($desc =~ /AGAP\d{6}\S+( )?/i) {
        $desc =~ s/AGAP\d{6}\S+( )?//i;
    }
    ## danio rerio
    if ($desc =~ /^Si:ch\d+-\d+\S+( )?/i) {
        $desc =~ s/^Si:ch\d+-\d+\S+( )?//i;
    }
    if ($desc =~ /^At\d+g[\d.]+ \(([^\)]+)\)$/i) {
        $desc = $1;
    }
    if ($desc =~ /^[a-zA-Z0-9]+v\d \(([^\)]+)\)$/i) {
        $desc = $1;
    }
    if ($desc =~ / \(Os\d+g[\d.]+[^)]+\)$/i) {
        $desc =~ s/ \(Os\d+g[\d.]+[^)]+\)$//i;
    }
    ## strip orf\d+ from the end
    if ($desc =~ / orf\d+( homolog)?$/i) {
        $desc =~ s/ orf\d+( homolog)?$//i;
    }
    if ($desc =~ /^At\d+g[\d.]+ \(([^\)]+)\)$/i) {
        $desc = $1;
    }
    if ($desc =~ /^[a-zA-Z0-9]+v\d \(([^\)]+)\)$/i) {
        $desc = $1;
    }
    if ($desc =~ / \(Os\d+g[\d.]+[^)]+\)$/i) {
        $desc =~ s/ \(Os\d+g[\d.]+[^)]+\)$//i;
    }
    if ($desc =~ / \([a-z]+\d{3,}[a-z]+\d+(\.\d+)?( protein)\)$/i) {
        $desc =~ s/ \([a-z]+\d{3,}[a-z]+\d+(\.\d+)?( protein)\)$//i;
    }
    if ($desc =~ /^Similar(ity)? to /i) {
        $desc =~ s/^Similar(ity)? to //i;
    }
    if ($desc =~ /^Contains similar(ity)? to /i) {
        $desc =~ s/^Contains similar(ity)? to //i;
    }
    if ($desc =~ /^((Highly|Partially|Very|Strong) )?(Gene, \S*, )?(Novel protein )?similar to/i) {
        $desc =~ s/^((Highly|Partially|Very|Strong) )?(Gene, \S*,)?(Novel protein )?similar to//i;
    }
    if ($desc =~ /\s*[- ]\s*like protein/i) {
        $desc =~ s/\s*[- ]\s*like protein//i;
    }
    if ($desc =~ /-like$/i) {
        $desc =~ s/-like$//i;
    }
    if ($desc =~ /, IPR\d{6}$/i) {
        $desc =~ s/, IPR\d{6}$//i;
    }
    if ($desc =~ /^Putative /i) {
        $desc =~ s/^Putative //i;
    }
    if ($desc =~ /^PREDICTED(:)? /i) {
        $desc =~ s/^PREDICTED(:)? //i;
    }
    if ($desc =~ /^Probable /i) {
        $desc =~ s/^Probable //i;
    }
    if ($desc =~ / \(fragment\)$/i) {
        $desc =~ s/ \(fragment\)$//i;
    }
    if ($desc =~ /^Tobacco /i) {      
        $desc =~ s/^Tobacco //i;
    }
    if ($desc =~ /; \d+-\d+( \([^)]+\))?$/) {      
        $desc =~ s/; \d+-\d+( \([^)]+\))?//;
    }
    if ($desc =~ / \(putative uncharacteri[zs]ed protein[^)]*\)/i) {      
        $desc =~ s/ \(putative uncharacteri[zs]ed protein[^)]*\)//ig;
    }
    if ($desc =~ / \([a-z]+\|\S+\)/i) {      
        $desc =~ s/ \([a-z]+\|\S+\)//ig;
    }
    if ($desc =~ /, putat(i)?ve/i) {
        $desc =~ s/, putat(i)?ve//i;
    }
    if ($desc =~ /, [35]'[ -]partial$/) {      
        $desc =~ s/, [35]'[ -]partial$//;
    }
    if ($desc =~ /, expressed$/i) {      
        $desc =~ s/, expressed$//i;
    }
    if ($desc =~ /, identical$/i) {      
        $desc =~ s/, identical$//i;
    }
    ## 'Elongation factor 1B alpha-subunit0like'
    if ($desc =~ /[\d-]+like$/) {
        $desc =~ s/[\d-]+like$//;
    }

## hopefully this will deal with:
## UniRef100_UPI0000ECD818 Dipeptidyl-peptidase 1 precursor (EC 3.4.14.1) (Dipeptidyl-peptidase I) (DPP-I) (DPPI) (Cathepsin C) (Cathepsin J) (Dipeptidyl transferase) [Contains: Dipeptidyl-peptidase 1 exclusion domain chain (Dipeptidyl- peptidase I exclusion domain chain);
    if ($desc =~ / \[Contains:.*$/i) {
        $desc =~ s/ \[Contains:.*$//ig;
    }

##### Kevin's additions

    if ($desc =~ /Jasmonate inducible protein isolog/) {
    # This one crazy case is not easily generalized.
    # Deal with it directly.
        $desc = "Jasmonate inducible protein isolog";
    }
    if ($desc =~ /Ferric reductase-like transmembrane component\)/) {
    # This one crazy case is not easily generalized.
    # Deal with it directly.
        $desc = "Ferric reductase-like transmembrane component";
    }
    if ($desc =~ /At2g39190\/T16B24.17 \(Putative ABC transporter;/) {
    # This one crazy case is not easily generalized.
    # Deal with it directly.
        $desc = "ABC transporter";
    }
    if ($desc =~ /^([^(]+)(\s+\(.+)$/i) {  
    # UDP-N-acetylglucosamine--dolichyl-phosphate N-acetylglucosaminephosphotransferase (UDP-GlcNAc:dolichol phosphate N-acetylglucosamine-1-phosphate transferase) (T26J13.8/T26J13.8) 
        my $first = $1;
        my $second = $2;
        if (length($first) > 30 && length($second) > 30) {
#            $second =~ s/\(/\\\(/g;
#            $second =~ s/\)/\\\)/g;
            $desc =~ s/\Q$second//;
        }
    }
    if ($desc =~ /^(.+)(;.+)$/i) {  
    # Cellular retinaldehyde binding/alpha-tocopherol transport; Cellular retinaldehyde-binding/triple function, N-terminal 
        my $first = $1;
        my $second = $2;
        if (length($first) > 30 && length($second) > 30) {
            $desc =~ s/\Q$second//i;
        }
    }
    if ($desc =~ /^(\[.+\]\s+)/i) {  
    # F24P17.21 protein (At3g06330/F24P17_21) 
        $desc =~ s/\[.+\]\s+//i;
    }
    if ($desc =~ /3\'(\'+)/i) {  ## 3''''-5''''  Duh-oh!
        $desc =~ s/3\'(\'+)/3\'/i;
    }
    if ($desc =~ /5\'(\'+)/i) {  ## 3''''-5''''
        $desc =~ s/5\'(\'+)/5\'/i;
    }
    if ($desc =~ /At\d+g[\d.]+\/[\w]+\s+\(([^)]+)\)/i) {  
    # AT5g20650/T1M15_50 (COPT5)
        my $match = $1;  # In this case = "COPT5"
        $desc =~ s/At\d+g[\d.]+\/[\w]+\s+\(([^)]+)\)/$match/i;
    }
    if ($desc =~ /Gb\|[\w.]+\s+\(([^)]+)\)/i) {  
    # Gb|AAF57656.1  (CDC48-interacting UBX-domain protein)
        my $match = $1;  # In this case = "COPT5"
        $desc =~ s/Gb\|[\w.]+\s+\(([^)]+)\)/$match/i;
    }
    if ($desc =~ /COG\d+\:\s+/i) {  
    # COG3188: P pilus assembly protein, porin PapC
        $desc =~ s/COG\d+\:\s+//i;
    }
    if ($desc =~ / ?\(Putative [^\(]+\)/i) {  
    # Enough with the Putatives.
    # NOI protein (Putative nitrate-induced NOI protein) (Putative NOI protein, nitrate-induced)
        $desc =~ s/ ?\(Putative [^\(]+\)//i;
    }
    if ($desc =~ /\)\s*\(.+/i) {  
    # By the second paranthetical expression, there is usually no additional useful information.
    # hydrolase/acyltransferase (Alpha/beta hydrolase superfamily) (ISS)
        $desc =~ s/\)\s*\(.+/\)/i;
    }
    if ($desc =~ /\s*\(Uncharacterized protein\)/i) {  
    # F12F1.32 protein (Uncharacterized protein)
        $desc =~ s/\s*\(Uncharacterized protein\)//i;
    }
    if ($desc =~ /\s*\([^\(^\)]+whole genome shotgun sequence\)/i) {  
    # Class IV chitinase (Chromosome chr5 scaffold_72, whole genome shotgun sequence)
        $desc =~ s/\s*\([^\(^\)]+whole genome shotgun sequence\)//i;
    }
    if ($desc =~ /\(RAP Annotation release2\)\s+/i) {  
    # Those guys.
        $desc =~ s/\(RAP Annotation release2\)\s+//i;
    }
    if ($desc =~ / \(Lycopersicum esculentum mRNA sequence\)/i) {  
        $desc =~ s/ \(Lycopersicum esculentum mRNA sequence\)//i;
    }
    ## need to hit this one twice apparently
    if ($desc =~ /^Similar to /i) {
        $desc =~ s/^Similar to //i;
    } 
    ## have only seen this once
    if ($desc =~ /; alternative splicing isoform gene prediction data combined with cDNA alignment data to generate this model/i) {
        $desc =~ s/; alternative splicing isoform gene prediction data combined with cDNA alignment data to generate this model//i;
    } 
    ## remove this sort of garbage from end of string: C6orf149
    if ($desc =~ /\s*C\d+orf\d+( homolog)?$/i) {
        $desc =~ s/\s*C\d+orf\d+( homolog)?$//i;
    }

#########################


    ## trim ends
    $desc =~ s/^\s+|\s+$//g;

    ## check again after cleanup
    if (_reject_check($desc)) {
        return undef;
    }

    $result = [
        $id,
        $desc,
        $species,
        $weight,
    ];

    #### Assign a weight to others
#    if ($desc =~ /(putative )?uncharacterized protein/i) {
#        #$result->[2] -= 1;
#        print STDERR "REJECTED: $desc\n";
#        $result = undef;
#    }

    return $result;
}

sub _reject_check {
    my ($desc) = @_;
    
    ## reject empty strings
    if ($desc =~ /^\s*$/) {
        return 1;
    }
    ## Vitis vinifera annotation is useless
    if ($desc =~ /^Chromosome.*whole genome shotgun sequence/i) {
        return 1;
    }
    ## Vitis vinifera annotation is useless
    elsif ($desc =~ /^Whole genome shotgun sequence of/i) {
        return 1;
    }
    ## reject putative uncharacterized protein
    elsif ($desc =~ /^Putative uncharacterized protein/i) {
        return 1;
    }
    ## reject the word Mitochondrion by itself
    elsif ($desc =~ /^Mitochondrion$/i) {
        return 1;
    }
    ## reject the word Protein by itself
    elsif ($desc =~ /^Protein$/i) {
        return 1;
    }
    ## reject the word Binding by itself
    elsif ($desc =~ /^Binding$/i) {
        return 1;
    }
    ## reject the word ORF by itself
    elsif ($desc =~ /^(Nicotiana tabacum )?ORF\d*( protein)?$/i) {
        return 1;
    }
    ## reject one character strings
    elsif ($desc =~ /^[A-Z0-9]$/i) {
        return 1;
    }
    ## Garbage Arabidopsis annotation
    elsif ($desc =~ /Genomic DNA,.*chromosome.*clone/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^CDNA( sequence)? [A-Z0-9]{6,}/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^(\S+ )?CDNA(,)? clone( \S+)?[:.]/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^\S+ CDNA \S+ fis$/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^CDNA: \S+ fis, clone \S+$/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^CDNA, \S+, highly similar to /i) {
        return 1;
    }
    ## RIKEN CDNA clone info
    elsif ($desc =~ /^RIKEN cDNA/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^MRNA, .*, clone:/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^CDNA-\S+-encoded protein$/i) {
        return 1;
    }
    ## no good way to match this in a generalized form 
    elsif ($desc =~ /^CDS localized after complete sequencing of a cognate cDNA$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^AT\d+g\d+(\/[\w.]+)?( protein)?$/i) {
        return 1;
    }
    ## locus-only annotation (unknown, eg: A_IG002N01.4)
    elsif ($desc =~ /^A_\w+\d+\w+\d+\.\d+ protein$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^CM\d+\..* protein$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^[A-Z]\d+[A-Z]\d+(\.\d+)?( protein)?$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^[A-Z]\d+[A-Z]\d+(\.\d+)?( \S+)?$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^\S+ protein \(At\d+g\d+\)$/i) {
        return 1;
    }
    ## locus-only annotation
    elsif ($desc =~ /^LOC\d+ protein$/i) {
        return 1;
    }
    elsif ($desc =~ /^Zgc:\d+ protein$/i) {
        return 1;
    }
    elsif ($desc =~ /^[A-Z]+\d+.\d+ protein$/) {
        return 1;
    }
    elsif ($desc =~ /^CG\d+\S+ isoform \d+$/i) {
        return 1;
    }
    elsif ($desc =~ /^\S+ specific protein [A-Z]+\d+[A-Z]+\S+$/i) {
        return 1;
    }
    elsif ($desc =~ /^isoform \d+ of Uncharacterized protein \S+$/i) {
        return 1;
    }
    ## Rice-locus-only annotation
    elsif ($desc =~ /^Os\S*\d{4,}\S*( protein)?( \(fragment\))?/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^(expressed|hypothetical|predicted|conserved|unknown|uncharacteri[sz]ed) protein( \(fragment\))?/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^unidentified$/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^protein binding$/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^protein binding protein$/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^Down syndrome/i) {
        return 1;
    }
    ## genbank ID only
    elsif ($desc =~ /^[a-z]+\|\S+$/i) {
        return 1;
    }
    ## generic locus ID only
    elsif ($desc =~ /^[a-z0-9]+\.\d+$/i) {
        return 1;
    }
    ## generic locus ID only --- hard to believe, but for tobacco: CYP81B2v2
    elsif ($desc =~ /^[a-zA-Z0-9]+v\d+$/) {
        return 1;
    }
    ## B0103C08-B0602B01.10 protein don't know where this is coming from, but hopefully this will catch all of them
    elsif ($desc =~ /^[a-z]+\d+[a-z]+\d+-[a-z]+\d+[a-z]+\d+(\.\d+)? protein$/i) {
        return 1;
    }
    ## Location of EST ... 
    elsif ($desc =~ /^Location of EST /i) {
        return 1;
    }
    ## UPIxxx-related cluster 
    elsif ($desc =~ /^UPI\S+ related cluster$/i) {
        return 1;
    }
    ## UPF0727 protein WS02710_H03
    elsif ($desc =~ /^UPF\d+ protein \S+$/i) {
        return 1;
    }
    ## Similar to ENSANGP00000021579
    elsif ($desc =~ /^Similar to [A-Z0-9]+$/i) {
        return 1;
    }
    ## Our own annotation hounds us??
    elsif ($desc =~ /Rice Genome Annotation Project/i) {
        return 1;
    }
    ## Sometimes a size measure is all that is left.  25.7 kDa protein
    elsif ($desc =~ /^[\d.]+ kDa protein$/i) {
        return 1;
    }
    ## Transient annotation trap.
    elsif ($desc =~ /^By genscan$/i) {
        return 1;
    }
    ## Arabidopsis locus-only annotation
    elsif ($desc =~ /^Protein F\d+F\d+\.\d+(-)?$/) {
        return 1;
    }
    ## Transient annotation trap.
    elsif ($desc =~ /related/i) {
        return 1;
    }
    # Isoform 2 of Uncharacterized protein  ???????????
    elsif ($desc =~ /^Uncharacteri[sz]ed (conserved )?protein/i) {
        return 1;
    }
    ## transitive domain annotation 
    elsif ($desc =~ /^Polypeptide with a(n)? (.*) domain$/i) {
        return 1;
    }
    ## the world 'catalytic' only 
    elsif ($desc =~ /^catalytic$/i) {
        return 1;
    }
#    ## the world 'polyprotein' only 
#    elsif ($desc =~ /^polyprotein$/i) {
#        return 1;
#    }
    elsif ($desc =~ /^(\[S\] KOG\d+)/i) {  
    ## Some sort of COG protein family.
        return 1;
    } else {
        return 0;
    }
}
