#!/bin/bash

TIMESTAMP=`date +%Y%m%d`
WORK_DIR=/projects/compsol/cron/$TIMESTAMP

TODAY=`date +%Y/%m/%d`
CONTIG_MAX_DATE=`~/SVN/solcomp/gb_max_date.pl contig`
READ_MAX_DATE=`~/SVN/solcomp/gb_max_date.pl read`

mkdir -p $WORK_DIR
cd $WORK_DIR

~/SVN/gb/fetch_genbank_gis_by_entrez_query.pl -db nuccore --query "txid4070[Organism:exp] AND 10000:10000000[SLEN] AND (all[filter] NOT \"wgs master\"[Properties]) AND (all[filter] NOT 41343[Genome Project]) AND (\"$CONTIG_MAX_DATE\"[PDAT] : \"$TODAY\"[PDAT])" -o sol.contig.new.list -l contig.list.log
~/SVN/gb/fetch_gbff_from_genbank_by_gi_list.pl --db nucleotide --input_list sol.contig.new.list -o sol.contig.new.gbff -l contig.gbff.log
~/SVN/gb/fetch_genbank_gis_by_entrez_query.pl -db nucgss --query "(bac ends[Library Class] OR cosmid ends[Library Class] OR fosmid ends[Library Class]) AND txid4070[Organism:exp] AND (\"$READ_MAX_DATE\"[PDAT] : \"$TODAY\"[PDAT])" -o sol.read.new.list -l read.list.log
~/SVN/gb/fetch_gbff_from_genbank_by_gi_list.pl --db nucleotide --input_list sol.read.new.list -o sol.read.new.gbff -l read.gbff.log

if [ -s sol.contig.new.gbff ];
then
#    ~/SVN/solcomp/add_gbff_contents_to_sol_seqdb.pl sol.contig.new.gbff
    echo "New contigs added to database."
else
    echo "No new contigs to add to database."
fi
if [ -s sol.contig.new.gbff ];
then
#    ~/SVN/solcomp/add_gbff_contents_to_sol_seqdb.pl sol.read.new.gbff
    echo "New reads added to database."
else
    echo "No new reads to add to database."
fi
