#!/bin/bash

#$ -N newick2png 
#$ -cwd
#$ -V
#$ -e /dev/null
#$ -o /dev/null

INPUT_FILE=`/bin/sed -n ${SGE_TASK_ID}p ${INPUT_LIST}`

export PERL5LIB='/opt/rocks/lib/perl5/site_perl/5.10.1/'

BINDIR='/run/whitty/potato/v3_2/orthomcl/bin/'
ANNOT_FILE='/run/whitty/potato/v3_2/orthomcl/all_annot.txt'

${BINDIR}/create_newick_tree_png_svg_imap.pl -t ${INPUT_FILE} -n ${ANNOT_FILE}
