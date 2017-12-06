#!/bin/bash

#$ -N treebuild
#$ -cwd

export PERL5LIB="/opt/rocks/lib/perl5/site_perl/5.10.1/"

INFILE=`/bin/sed -n ${SGE_TASK_ID}p ${INPUT_LIST}`
PREFIX=`basename ${INFILE} .fasta`
DIRNAME=`dirname ${INFILE}`
TMPDIR="/state/partition1/tmp/tmp_`date +%Y%m%d`_$$"

BINDIR='/run/whitty/potato/v3_2/orthomcl/bin'
MAPPINGS='/run/whitty/potato/v3_2/orthomcl/14_plants.mappings.list'

## change to location of input file
cd ${DIRNAME}

echo "## working dir '${DIRNAME}'" 1>&2
echo "##  input base '${PREFIX}'" 1>&2

## cleanup previous output files
rm ${DIRNAME}/${PREFIX}.phyi
rm ${DIRNAME}/${PREFIX}.phyi.trim
rm ${DIRNAME}/${PREFIX}.gfsa
rm ${DIRNAME}/${PREFIX}.gfsa.trim
rm ${DIRNAME}/${PREFIX}.gfsa.trim.txt
rm ${DIRNAME}/${PREFIX}.newick
rm ${DIRNAME}/${PREFIX}.msa
rm ${DIRNAME}/${PREFIX}.xml
rm ${DIRNAME}/${PREFIX}.final.newick
rm ${DIRNAME}/${PREFIX}.codes.txt

## run muscle
echo "## running muscle" 1>&2
/share/apps/bin/muscle -in ${DIRNAME}/${PREFIX}.fasta -quiet -phyiout ${DIRNAME}/${PREFIX}.phyi -maxmb 2000

## check how many sequences are in the alignment
SEQCOUNT=`cut -f 1 -d " " ${DIRNAME}/${PREFIX}.phyi | head -1`
echo "## alignment contains '${SEQCOUNT}' sequences" 1>&2

## run readseq to convert MSA to input format for Gblocks
echo "## converting MSA to Gblocks format input" 1>&2
java -jar /share/apps/readseq.jar -inform=12 -f=8 -o ${DIRNAME}/${PREFIX}.gfsa ${DIRNAME}/${PREFIX}.phyi

## run Gblocks to trim MSA
echo "## running Gblocks" 1>&2
/share/apps/bin/Gblocks ${DIRNAME}/${PREFIX}.gfsa -t=p -e=.trim -p=t -b5=h

## convert Gblocks output to input for phylip
echo "## converting Gblocks output to phylip format" 1>&2
java -jar /share/apps/readseq.jar -inform=8 -f=12 -o ${DIRNAME}/${PREFIX}.phyi.trim ${DIRNAME}/${PREFIX}.gfsa.trim

## make the temp dir
mkdir ${TMPDIR}

## cd into temp dir
cd ${TMPDIR}
echo "## cwd: '`pwd`'"

## copy phylip format alignment into temp dir
cp ${DIRNAME}/${PREFIX}.phyi.trim infile

## create "parameter file" for phylip
echo "Y" | cat > ${TMPDIR}/parameter_file.txt

if [ ${SEQCOUNT} -lt 3 ] ; then
    ## need to handle two seqs differently because we can't make a tree with proml
    echo "## running protdist" 1>&2
    /opt/bio/phylip/exe/protdist <parameter_file.txt
    ${BINDIR}/create_newick_file.pl -i ${TMPDIR}/outfile -o ${TMPDIR}/outtree 
else
    ## run proml
    echo "## running proml" 1>&2
    /opt/bio/phylip/exe/proml <parameter_file.txt
fi

## convert newick format tree to phyloXML
echo "## converting newick to phyloXML" 1>&2
java -cp /share/apps/forester_dev.jar org.forester.application.phyloxml_converter -f=gn ${TMPDIR}/outtree ${DIRNAME}/${PREFIX}.xml

## copy the newick tree to the output dir
cp ${TMPDIR}/outtree ${DIRNAME}/${PREFIX}.newick

## convert the numeric IDs back to gene names in the XML
echo "## converting numeric IDs to gene IDs in XML" 1>&2
${BINDIR}/decode_sequence_ids.pl -i ${DIRNAME}/${PREFIX}.xml -d ${MAPPINGS}

## convert the numeric IDs in the MSA file to gene names
echo "## converting numeric IDs to gene IDs in MSA file" 1>&2
${BINDIR}/convert_and_decode_phyi_file.pl -i ${DIRNAME}/${PREFIX}.phyi  -o ${DIRNAME}/${PREFIX}.msa -d ${MAPPINGS}

## convert the numeric IDs in the newick file to gene names
echo "## converting numeric IDs to gene IDs in newick file" 1>&2
${BINDIR}/orthomcl_tree_id_remapper.pl ${MAPPINGS} ${DIRNAME}/${PREFIX}.newick

## remove the temp dir
rm -rf ${TMPDIR}
