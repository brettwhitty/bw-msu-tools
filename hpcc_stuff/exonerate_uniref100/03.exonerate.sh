#!/bin/bash

#$ -cwd
#$ -e /dev/null
#$ -o /dev/null

export PERL5LIB=/opt/rocks/lib/perl5/site_perl/5.10.1/

## run exonerate
/share/apps/bin/exonerate --query $OUTPREFIX.db --target $QUERY --model $MODEL --verbose false --showalignment false --showsugar false --showquerygff false --showtargetgff true --showcigar true --showvulgar false --ryo "%ti\t%qi\t%ql\t%qal\t%r\t%s\t%pi\t%ps\n" --percent 20 --softmaskquery true --softmasktarget true --maxintron $MAXINTRON --minintron $MININTRON --subopt false --fsmmemory 1024 --dpmemory 1024 --forcefsm normal --wordjump 2 --hspfilter 100 --geneseed 250 --geneseedrepeat 4 1>$OUTPREFIX.exonerate.raw 2>$OUTPREFIX.exonerate.log

## convert raw output to gff
/home/whitty/SVN/gff/convert_to_gff/convert_exonerate_gff_to_gff3.pl -i $OUTPREFIX.exonerate.raw -s "$MODEL" -m "match" 1>$OUTPREFIX.exonerate.gff 2>/dev/null

## filter gff to HITLIMIT tiers
/home/whitty/SVN/gff/tier_gff3.pl -t $HITLIMIT -i $OUTPREFIX.exonerate.gff 1>$OUTPREFIX.exonerate.${HITLIMIT}_tier.gff 2>/dev/null

## remove database
rm $OUTPREFIX.db

touch $OUTPREFIX.END
