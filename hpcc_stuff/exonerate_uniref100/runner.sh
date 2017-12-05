#!/bin/bash

#$ -cwd
#$ -N xonr8f
##$ -e /dev/null
##$ -o /dev/null

export TOOLPATH="/home/whitty/tools/exonerate_filter/uniref100"

export THREADS=4

export QUERY=`/bin/sed -n ${SGE_TASK_ID}p ${INPUT_LIST}`
export DB='/db/whitty/uniref100.fasta'

## params
export HITLIMIT=5
export MININTRON=10
export MAXINTRON=30000
export WORDSIZE=3
export EXPECT=0.01
export MODEL='p2g'
export PERCENT=20

## chunk params
export SEGLEN=25000
export SEGMIN=10000

## unique id used for subtasks
export UNAME="e$JOB_ID-$SGE_TASK_ID"

## get query length
FLEN=`/share/apps/bin/fastalength $QUERY`

## file prefix for output files
export OUTPREFIX=`echo "$FLEN" | cut -f 2 -d " "`

### CHUNKER PART
export QLEN=`echo "$FLEN" | cut -f 1 -d " "`
SEGS=$(($QLEN / $SEGLEN))
SMOD=$(($QLEN % $SEGLEN))

## created another segment from remainder if it >SEGMIN
if [[ $SMOD > $SEGMIN ]]; then 
    SEGS=$(( $SEGS + 1))
fi
export SEGS

touch $OUTPREFIX.START

## iterated blast filter step
if [[ $THREADS > 1 ]]; then
    qsub -cwd -pe smp $THREADS -N "$UNAME.1" -t 1-$SEGS -V $TOOLPATH/01.blast_filter.sh
else
    qsub -cwd -N "$UNAME.1" -t 1-$SEGS -V $TOOLPATH/01.blast_filter.sh
fi

## merge results and prepare db
qsub -cwd -N "$UNAME.2" -hold_jid "$UNAME.1" -V $TOOLPATH/02.prepare_db.sh

## run exonerate
qsub -cwd -N "$UNAME.3" -hold_jid "$UNAME.2" -V $TOOLPATH/03.exonerate.sh
