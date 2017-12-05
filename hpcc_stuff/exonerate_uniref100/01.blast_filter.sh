#!/bin/bash

#$ -cwd
#$ -e /dev/null
#$ -o /dev/null

BVLIM=$(( $HITLIMIT * 10 ))

## calculate start and end positions
ITER=$(($SGE_TASK_ID - 1))
SEGS=$SGE_TASK_LAST

START=$(($ITER * $SEGLEN + 1))
if [[ $SGE_TASK_ID == $SEGS ]]; then
    END=$QLEN
else
    END=$(( ($ITER + 1) * $SEGLEN ))
fi

## run blast
#/opt/bio/ncbi/bin/blastall -L "$START,$END" -a $THREADS -p blastx -d $DB -i $QUERY -F "m S" -U T -m 8 -K $HITLIMIT -t $MAXINTRON -n T -e $EXPECT -W $WORDSIZE -o $OUTPREFIX.$SGE_TASK_ID.blastx.raw -z 1000000000

#/opt/bio/ncbi/bin/blastall -L "$START,$END" -a $THREADS -p blastx -d $DB -i $QUERY -F "m S" -U T -K $HITLIMIT -t $MAXINTRON -n T -e $EXPECT -W $WORDSIZE -o $OUTPREFIX.$SGE_TASK_ID.blastx.raw -z 1000000000 -b $BVLIM -v $BVLIM

/opt/bio/ncbi/bin/blastall -L "$START,$END" -a $THREADS -p blastx -d $DB -i $QUERY -F "m S" -U T -K $HITLIMIT -t $MAXINTRON -n T -e $EXPECT -W $WORDSIZE -o $OUTPREFIX.$SGE_TASK_ID.blastx.raw -z 1000000000 -b $BVLIM -v $BVLIM -g F -P 1

#echo "/share/apps/ncbi-blast-2.2.23+/bin/blastx -query $QUERY -db $DB -word_size $WORDSIZE -num_threads $THREADS -culling_limit $HITLIMIT -lcase_masking -seg yes -soft_masking true -max_intron_length $MAXINTRON  -query_loc \"$START-$END\" -outfmt 6"
