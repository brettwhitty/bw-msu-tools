#!/bin/bash

#$ -cwd
#$ -e /dev/null
#$ -o /dev/null
##$ -l h_vmem=4G
##$ -pe smp 2
#$ -t 1-6946
#$ -N blastp
#$ -R y

export WUBLASTMAT=/share/apps/wublast/matrix
#export WUBLASTDB=/share/db

/share/apps/wublast/blastp /db/whitty/kegg/ko_genes.faa /db/whitty/potato/v3_2/pep/PGSC0003DMP.pep_$SGE_TASK_ID.fa -warnings -notes -shortqueryok -novalidctxok -cpus 1 -e 1e-5 -b 1 -v 1 -hitdist 40 -wordmask "seg" -W 3 -T 999 -matrix BLOSUM80 -Q 11 -R 2 -noseqs 2>/dev/null | gzip -c >PGSC0003DMP.ko.blastp.$SGE_TASK_ID.raw.gz
