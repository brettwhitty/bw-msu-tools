#!/bin/bash

#$ -N omcl2tree
#$ -cwd
#$ -e /dev/null
#$ -o /dev/null

INPUT_FILE=`/bin/sed -n ${SGE_TASK_ID}p ${INPUT_LIST}`
#SUFFIX=`basename ${INPUT_FILE}`

export PERL5LIB=/opt/rocks/lib/perl5/site_perl/5.8.8

source ${INPUT_FILE}
