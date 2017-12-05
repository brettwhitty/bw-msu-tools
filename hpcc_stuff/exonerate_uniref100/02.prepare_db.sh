#!/bin/bash

#$ -cwd
##$ -e /dev/null
##$ -o /dev/null

# OUTPREFIX
PERL5LIB=/opt/rocks/lib/perl5/site_perl/5.10.1/

## merge the output
cat $OUTPREFIX.*.blastx.raw >$OUTPREFIX.blastx.raw

## remove the unmerged files
#rm $OUTPREFIX.*.blastx.raw

## convert the blast output
/home/whitty/SVN/bp_search2gff.pl -v 3 -i $OUTPREFIX.blastx.raw -method match_part -m -s "blast" 1>$OUTPREFIX.blastx.gff 2>/dev/null

## remove merged output file
#rm $OUTPREFIX.blastx.raw

## tier the blast output
/home/whitty/SVN/gff/tier_gff3.pl -t $HITLIMIT -i $OUTPREFIX.blastx.gff | perl -ne 'chomp; if (/Name=([^;]+)/) { print $1."\n";}' | sort | uniq | sed -r "s/lcl\|//" | /run/whitty/potato/v3/v3_exonerate/blastx/extract_sequence_db_by_list.pl --db $DB 1>$OUTPREFIX.db 2>/dev/null

## remove gff file
#rm $OUTPREFIX.blastx.gff

## get hit IDs and retrieve from database
#cat $OUTPREFIX.*.blastx.raw | egrep -v "^#" | cut -f 2 | sort | uniq | sed -r "s/lcl\|//" | /run/whitty/potato/v3/v3_exonerate/blastx/extract_sequence_db_by_list.pl --db $DB 1>$OUTPREFIX.db 2>/dev/null
