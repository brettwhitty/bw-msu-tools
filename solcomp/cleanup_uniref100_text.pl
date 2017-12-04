#!/usr/bin/perl

use strict;
use warnings;

use lib "/home/whitty/SVN/lib";
use MyIO;

my $uniref_file = shift @ARGV || die "Provide path to uniref100 fasta";

my $infh = get_infh($uniref_file);

while (<$infh>) {

    if (/^>((\S+) .*)/) {
        my ($header, $id) = ($1, $2);

        my $annotation_string = get_uniref_function($header);

        unless( defined($annotation_string) && $annotation_string ne '' ) {
            $annotation_string = 'Gene of unknown function';
        }

        print join("\t", (
            $id,
            ucfirst($annotation_string),
            #$annotation_string,
        ))."\n";

    } else {
        next;
    }
}


sub get_uniref_function {
    my ($header) = @_;

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
    if ($desc =~ / At\d+g[\d.]+/i) {
        $desc =~ s/ At\d+g[\d.]+//i;
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
    if ($desc =~ /\s*[- ]\s*like protein/i) {
        $desc =~ s/\s*[- ]\s*like protein//i;
    }
    if ($desc =~ /-like$/i) {
        $desc =~ s/-like$//i;
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
    if ($desc =~ / \(fragment(s)?\)$/i) {
        $desc =~ s/ \(fragment(s)?\)$//i;
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
    if ($desc =~ / homolog(ue)?((?=\s)|(?=$))/i) {      
        $desc =~ s/ homolog(ue)?//i;
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

#########################


    ## trim ends
    $desc =~ s/\s+$//;
    $desc =~ s/^\s+//;

    ## check again after cleanup
    if (_reject_check($desc)) {
        return undef;
    }

    #### Assign a weight to others
#    if ($desc =~ /(putative )?uncharacterized protein/i) {
#        #$result->[2] -= 1;
#        print STDERR "REJECTED: $desc\n";
#        $result = undef;
#    }

    return $desc;
}

sub _reject_check {
    my ($desc) = @_;

    ## Vitis vinifera annotation is useless
    if ($desc =~ /^Chromosome.*whole genome shotgun sequence/i) {
        return 1;
    }
    ## Garbage Arabidopsis annotation
    elsif ($desc =~ /Genomic DNA,.*chromosome.*clone/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^CDNA(,)? clone:/i) {
        return 1;
    }
    ## CDNA clone info
    elsif ($desc =~ /^MRNA, .* cds, clone:/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^AT\d+g\d+(\/[\w.]+)?$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^[A-Z]\d+[A-Z]\d+(\.\d+)?( protein)?$/i) {
        return 1;
    }
    ## Arabidopsis-locus-only annotation
    elsif ($desc =~ /^\S+ protein \(At\d+g\d+\)$/i) {
        return 1;
    }
    ## Rice-locus-only annotation
    elsif ($desc =~ /^Os\S*\d{4,}\S* protein( \(fragment\))?/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^(expressed|hypothetical|predicted|conserved|unknown|uncharacteri[sz]ed) protein( \(fragment\))?/i) {
        return 1;
    }
    ## predicted protein
    elsif ($desc =~ /^protein binding$/i) {
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
    ## UPIxxx-related cluster 
    elsif ($desc =~ /^UPI\S+ related cluster$/i) {
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
    elsif ($desc =~ /related/i) {
        return 1;
    }
    # Isoform 2 of Uncharacterized protein  ???????????
    elsif ($desc =~ /^Uncharacteri[sz]ed (conserved )?protein/i) {
        return 1;
    }
    elsif ($desc =~ /^(\[S\] KOG\d+)/i) {  
    ## Some sort of COG protein family.
        return 1;
    } else {
        return 0;
    }
}
