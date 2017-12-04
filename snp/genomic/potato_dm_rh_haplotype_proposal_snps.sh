#!/bin/bash

#$ -R y
#$ -pe smp 8
#$ -cwd
#$ -N bwtsamsnp

#sample qsub script alignment to ref with Bowtie (with PE and SE lanes) and consensus and SNP calling with SAMtools

export PATH=$PATH:/share/apps/bowtie:/share/apps/samtools/bin

#GENOME_GFF=/projects/potato/genome_0428/PGSC0003DMS.gene.tophat_input.gtf
GENOME_BOWTIE=/projects/potato/genome_0428/PGSC0003DMS
GENOME_FILE=/projects/potato/genome_0428/PGSC0003DMS.fasta2

SAMPLE=`basename ${INPUT_FILE} .fastq`

#GENOME_FILE="potato_dm.fa"
#READ_PATH="/share/data1/hamilton/potato_reads/${VAR}"
#GENOME_PATH="/share/data1/hamilton/potato_genomes/"

#cd ${READ_PATH}
#symlink the seq
#ln -s "${GENOME_PATH}${GENOME_FILE}" ${GENOME_FILE} 

echo "aligning reads to reference"

#echo "aligning PE reads"
#/share/apps/bowtie/bowtie -q --phred64-quals -X 400 -k 1 -m 1 -t --un ${SAMPLE}.unmapped.fastq --max ${SAMPLE}.maxhits.fastq --sam --sam-nohead ${GENOME_BOWTIE} -1 s_1_pe_1_filter.fastq,s_4_pe_1_filter.fastq -2 s_1_pe_2_filter.fastq,s_4_pe_2_filter.fastq ${SAMPLE}.sam


if [[ $NSLOTS == 8 ]]
then
echo "aligning SE reads"
/share/apps/bowtie/bowtie -p 8 -q --phred64-quals -X 400 -k 1 -m 1 -t --un ${SAMPLE}.unmapped.fastq --max ${SAMPLE}.maxhits.fastq --sam --sam-nohead ${GENOME_BOWTIE} ${INPUT_FILE} ${SAMPLE}.sam
qalter -pe smp 1 $JOB_ID
exit 99
fi

echo "index reference fasta"
/share/apps/samtools/bin/samtools faidx ${GENOME_FILE}

echo "convert SAM to BAM"
samtools import ${GENOME_FILE}.fai ${SAMPLE}.sam ${SAMPLE}.bam

echo "sort BAM"
samtools sort ${SAMPLE}.bam ${SAMPLE}.bam.sort
mv ${SAMPLE}.bam.sort.bam ${SAMPLE}.bam

echo "index BAM"
/share/apps/samtools/bin/samtools index ${SAMPLE}.bam

#echo "create consensus pileup"
#/share/apps/samtools/bin/samtools pileup -cf ${GENOME_FILE} ${SAMPLE}.bam > ${SAMPLE}.pileup

echo "create SNP call pileup"
/share/apps/samtools/bin/samtools pileup -vcf ${GENOME_FILE} ${SAMPLE}.bam > ${SAMPLE}.snp

#echo "create consensus sequence"
#/share/apps/samtools/bin/samtools.pl pileup2fq ${SAMPLE}.pileup > ${SAMPLE}.cns.fastq

echo "create filtered snps"
/share/apps/samtools/bin/samtools.pl varFilter -d 20 -D 240 -W 100 -N 2 -w 50 ${SAMPLE}.snp > ${SAMPLE}.filtered.snp
