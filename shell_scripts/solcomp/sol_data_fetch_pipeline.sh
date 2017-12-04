#!/bin/bash

## bioperl
export PERL5LIB="/scratch/whitty/gmod/bioperl-live"

## fail on errors
set -e
## set a umask that will allow group rw
umask 006;


## set working directory root
WORK_DIR_ROOT=/projects/solcomp/data_working


## scripts used
FETCH_GENBANK_GIS="/home/whitty/SVN/gb/fetch_genbank_gis_by_entrez_query.pl"
FETCH_GENBANK_FLATFILES="/home/whitty/SVN/gb/fetch_gbff_from_genbank_by_gi_list.pl"
SPLIT_GENBANK_FLATFILES="/home/whitty/SVN/gb/split_genbank_flat_files.pl"
MERGE_GENBANK_FLATFILES="/home/whitty/SVN/gb/merge_gb_by_list_file.pl"
GENBANK_TO_GFF3="/home/whitty/SVN/bp_genbank2gff3_wrapper.pl"
FETCH_PUTS="/home/whitty/SVN/mirror_plantgdb_puts_by_taxon.pl"
PUTS_FASTA_TO_GFF3="/home/whitty/bin/gmod_fasta2gff3.pl" ## copied from GMOD installation
TAXONOMY_TEST="/home/whitty/SVN/taxonomy/taxon_id_is_solanaceae.pl"
GET_SCIENTIFIC_NAME="/home/whitty/SVN/taxonomy/get_scientific_name.pl"
SCP_FILES="/home/whitty/SVN/scp_sequences_to_ftp_site.pl"


## flag for allowing resuming of previous runs
RESUME=false

if [ "$1" != "" ]; then
    TIMESTAMP="$1"
    RESUME=true
else
    ## set a timestamp
    TIMESTAMP=`date +"%Y%m%d"`
fi

## set a prefix for the data set
SOL_BAC_PREFIX="solanaceae_bacs_$TIMESTAMP"
SOL_BAC_ENDS_PREFIX="solanaceae_bac_ends_$TIMESTAMP"

## create a work directory for this timestamp
WORK_DIR="$WORK_DIR_ROOT/$SOL_BAC_PREFIX"

if [ $RESUME == false ]; then
    mkdir -p $WORK_DIR
else
    if [ ! -d $WORK_DIR ]; then
        echo "Can't resume run with work dir '$WORK_DIR' because it doesn't exist"
        exit 1
    fi
fi

## change to the working directory
cd $WORK_DIR

## Download Solanaceae BAC sequences from GenBank
BAC_GI_LIST="$SOL_BAC_PREFIX.gi.list"
BAC_GI_LOG="$SOL_BAC_PREFIX.gi.log"
BAC_GBFF_FILE="$SOL_BAC_PREFIX.gbk.gbff"
BAC_GBFF_LOG="$SOL_BAC_PREFIX.gbk.log"
if [ ! -e ".checkpoint_1" ]; then
    ## fetch GI list
    trap "echo 'FAILED during GI list fetch'; " INT TERM EXIT
    echo "Fetching GI list..."
    $FETCH_GENBANK_GIS --db nuccore --output "$BAC_GI_LIST" --log "$BAC_GI_LOG" --query "txid4070[Organism:exp] AND 20000:10000000[SLEN]"
    trap - INT TERM EXIT
    touch ".checkpoint_1"
else 
    echo "== BAC GI list fetching checkpoint complete" 2>&1
fi
if [ ! -e ".checkpoint_2" ]; then
    ## fetch GenBank records
    trap "echo 'FAILED while fetching GenBank records'; " INT TERM EXIT
    echo "Fetching GenBank records..."
    $FETCH_GENBANK_FLATFILES --db nuccore --input_list "$BAC_GI_LIST" --output "$BAC_GBFF_FILE" --log "$BAC_GBFF_LOG"
    touch ".checkpoint_2"
else 
    echo "== BAC genbank file fetching checkpoint complete" 2>&1   
fi
if [ ! -e ".checkpoint_3" ]; then
    ## split GenBank records by organism
    trap "echo 'FAILED while splitting BAC GenBank records and creating fasta files'; " INT TERM EXIT
    echo "Splitting BAC GenBank records..."
    $SPLIT_GENBANK_FLATFILES $WORK_DIR/bacs <$BAC_GBFF_FILE
    trap - INT TERM EXIT
    touch ".checkpoint_3"
else 
    echo "== BAC genbank file splitting checkpoint complete" 2>&1   
fi

## verify taxonomy of retrieved sequences
if [ ! -e ".checkpoint_4" ]; then
    ## double-check that we haven't fetched any non-solanaceae sequence
    trap "echo 'FAILED while removing non-Solanaceae BAC sequences'; " INT TERM EXIT
    echo "Verifying taxonomy of retrieved BAC sequence files..."
    cd $WORK_DIR/bacs

    for FILE in *
    do
      if [ -d $FILE ]
      then
        echo "Verifying that '$FILE' is in the Solanaceae..."
        RESULT=`$TAXONOMY_TEST $FILE`
        if [ "$RESULT" == "false" ]
        then
            echo "Directory '$FILE' is not in the Solanaceae, removing..."
            if [ -d "$WORK_DIR/bacs/$FILE" ]
            then
                rm -rf "$WORK_DIR/bacs/$FILE"
            else
                echo "Can't remove directory '$WORK_DIR/bacs/$FILE' because it doesn't exist!"
                exit 1
            fi
        fi
      fi
    done
    ## return to working directory
    cd $WORK_DIR
    
    trap - INT TERM EXIT
    touch ".checkpoint_4"
else 
    echo "== BAC taxonomy verification checkpoint complete" 2>&1
fi

## merge files by chromosome/organelle
if [ ! -e ".checkpoint_5" ]; then
    ## double-check that we haven't fetched any non-solanaceae sequence
    trap "echo 'FAILED while merging BAC sequence files'; " INT TERM EXIT
    echo "Merging BAC sequence files..."
    cd $WORK_DIR/bacs

    for FILE in *
    do
      if [ -d $FILE ]
      then
        echo "Merging files in '$FILE'..."
        $MERGE_GENBANK_FLATFILES "$WORK_DIR/bacs/$FILE"
      fi
    done
    ## return to working directory
    cd $WORK_DIR
    
    trap - INT TERM EXIT
    touch ".checkpoint_5"
else 
    echo "== Merging of BAC files by replicon checkpoint complete" 2>&1
fi

## index BACs with xdformat
if [ ! -e ".checkpoint_6" ]; then
    trap "echo 'FAILED while indexing BACs with xdformat'; " INT TERM EXIT
    echo "Indexing BACs with xdformat..."
    find $WORK_DIR/bacs/ -name "*.fsa" -exec xdformat -I -n {} \;
    trap - INT TERM EXIT
    touch ".checkpoint_6"
else 
    echo "== Indexing of BACs checkpoint complete" 2>&1
fi

## convert BACs to GFF using bp_genbank2gff3.pl
##
if [ ! -e ".checkpoint_7" ]; then
    trap "echo 'FAILED while converting BAC GenBank flat files to GFF3'; " INT TERM EXIT
    echo "Converting BAC GenBank flat files to GFF3..."
    TEMP_GBK_LIST=`mktemp`
    find $WORK_DIR/bacs/ -name '*.gbk' 1>$TEMP_GBK_LIST 2>/dev/null
    while read GBK_FILE
    do
        BAC_GFF_DIR=`echo $GBK_FILE | sed -r "s/\.gbk/\.gff/"`
        $GENBANK_TO_GFF3 "--noCDS --nolump --outdir $BAC_GFF_DIR --typesource contig" $GBK_FILE 1>/dev/null 2>/dev/null
    done < $TEMP_GBK_LIST
    rm $TEMP_GBK_LIST
    trap - INT TERM EXIT
    touch ".checkpoint_7"
else 
    echo "== Conversion of BAC GenBank flat files to GFF3 complete" 2>&1
fi


## BAC ENDS VARIABLES
BAC_ENDS_GI_LIST="$SOL_BAC_ENDS_PREFIX.gi.list"
BAC_ENDS_GI_LOG="$SOL_BAC_ENDS_PREFIX.gi.log"
BAC_ENDS_GBFF_FILE="$SOL_BAC_ENDS_PREFIX.gbk.gbff"
BAC_ENDS_GBFF_LOG="$SOL_BAC_ENDS_PREFIX.gbk.log"

if [ ! -e ".checkpoint_11" ]; then
    ## fetch GI list
    trap "echo 'FAILED during GI list fetch'; " INT TERM EXIT
    echo "Fetching GI list..."
    $FETCH_GENBANK_GIS --db nucgss --output "$BAC_ENDS_GI_LIST" --log "$BAC_ENDS_GI_LOG" --query "txid4070[organism:exp] AND bac ends[Library Class]"
    trap - INT TERM EXIT
    touch ".checkpoint_11"
else 
    echo "== BAC ends GI list fetching checkpoint complete" 2>&1
fi
if [ ! -e ".checkpoint_12" ]; then
    ## fetch GenBank records
    trap "echo 'FAILED while fetching GenBank records'; " INT TERM EXIT
    echo "Fetching GenBank records..."
    $FETCH_GENBANK_FLATFILES --db nuccore --input_list "$BAC_ENDS_GI_LIST" --output "$BAC_ENDS_GBFF_FILE" --log "$BAC_ENDS_GBFF_LOG"
    touch ".checkpoint_12"
else 
    echo "== BAC ends genbank file fetching checkpoint complete" 2>&1 
fi
if [ ! -e ".checkpoint_13" ]; then
    ## split GenBank records by organism
    trap "echo 'FAILED while splitting GenBank records and creating fasta files'; " INT TERM EXIT
    echo "Splitting GenBank records..."
    $SPLIT_GENBANK_FLATFILES $WORK_DIR/bac_ends 1 <$BAC_ENDS_GBFF_FILE
    trap - INT TERM EXIT
    touch ".checkpoint_13"
else 
    echo "== BAC ends genbank file splitting checkpoint complete" 2>&1
fi

if [ ! -e ".checkpoint_14" ]; then
    ## double-check that we haven't fetched any non-solanaceae sequence
    trap "echo 'FAILED while removing non-Solanaceae BAC end sequences'; " INT TERM EXIT
    echo "Verifying taxonomy of retrieved BAC end sequence files..."
    cd $WORK_DIR/bac_ends

    for FILE in *
    do
      if [ -d $FILE ]
      then
        echo "Verifying that '$FILE' is in the Solanaceae..."
        RESULT=`$TAXONOMY_TEST $FILE`
        if [ "$RESULT" == "false" ]
        then
            echo "Directory '$FILE' is not in the Solanaceae, removing..."
            if [ -d "$WORK_DIR/bac_ends/$FILE" ]
            then
                rm -rf "$WORK_DIR/bac_ends/$FILE"
            else
                echo "Can't remove directory '$WORK_DIR/bac_ends/$FILE' because it doesn't exist!"
                exit 1
            fi
        fi
      fi
    done
    ## return to working directory
    cd $WORK_DIR
    
    trap - INT TERM EXIT
    touch ".checkpoint_14"
else 
    echo "== BAC ends taxonomy verification checkpoint complete" 2>&1
fi

## index bac ends with xdformat
if [ ! -e ".checkpoint_15" ]; then
    ## split GenBank records by organism
    trap "echo 'FAILED while indexing BAC ends with xdformat'; " INT TERM EXIT
    echo "Indexing BAC ends with xdformat..."
    find $WORK_DIR/bac_ends/ -name "*.fsa" -exec xdformat -I -n {} \;
    trap - INT TERM EXIT
    touch ".checkpoint_15"
else 
    echo "== Indexing of BAC ends checkpoint complete" 2>&1
fi


## convert BAC ends to GFF using bp_genbank2gff3.pl
##
if [ ! -e ".checkpoint_16" ]; then
    trap "echo 'FAILED while converting BAC end GenBank flat files to GFF3'; " INT TERM EXIT
    echo "Converting BAC GenBank end flat files to GFF3..."
    TEMP_GBK_LIST=`mktemp`
    find $WORK_DIR/bac_ends/ -name "*.gbk" 1>$TEMP_GBK_LIST 2>/dev/null
    while read GBK_FILE
    do
        echo $GBK_FILE
        BAC_ENDS_GFF_DIR=`echo $GBK_FILE | sed -r "s/\.gbk/\.gff/"`
        $GENBANK_TO_GFF3 "--noCDS --nolump --outdir $BAC_ENDS_GFF_DIR --typesource read" $GBK_FILE 2>/dev/null 1>/dev/null
    done < $TEMP_GBK_LIST
    rm $TEMP_GBK_LIST
    trap - INT TERM EXIT
    touch ".checkpoint_16"
else 
    echo "== Conversion of BAC ends GenBank flat files to GFF3 complete" 2>&1
fi


## Fetch PlantGDB PUTs for the Solanaceae
PUTS_DIR="$WORK_DIR/puts"

if [ ! -e ".checkpoint_21" ]; then
    if [ ! -d "$PUTS_DIR" ]; then
        mkdir -p $PUTS_DIR
    fi
    ## split GenBank records by organism
    trap "echo 'FAILED while fetching PlantGDB PUTs'; " INT TERM EXIT
    echo "Fetching PlantGDB PUTs.."
    $FETCH_PUTS -r $PUTS_DIR -t
    trap - INT TERM EXIT
    touch ".checkpoint_21"
else 
    echo "== PlantGDB PUTs fetching checkpoint complete" 2>&1
fi

## index PUTs with xdformat
if [ ! -e ".checkpoint_22" ]; then
    trap "echo 'FAILED while indexing PUTs with xdformat'; " INT TERM EXIT
    echo "Indexing PUTs with xdformat..."
    find $WORK_DIR/puts/ -name "*.PUT.fasta" -exec xdformat -I -n {} \;
    trap - INT TERM EXIT
    touch ".checkpoint_22"
else 
    echo "== Indexing of PUTs checkpoint complete" 2>&1
fi

## convert PUTs to GFF3
if [ ! -e ".checkpoint_23" ]; then
    trap "echo 'FAILED while creating GFF3 for PUTs'; " INT TERM EXIT
   
    echo "Converting PUTs fasta to GFF3..."
    
    cd $WORK_DIR/puts
    for FILE in *
    do
        if [ -d $FILE ]; then
            TAXON_ID=`echo $FILE | sed -r "s/\..*//"`
            SCIENTIFIC_NAME=`$GET_SCIENTIFIC_NAME $TAXON_ID`
            PUTS_FASTA_FILE=`find $FILE -name "*.PUT.fasta"`
            PUTS_GFF3_FILE=`echo $PUTS_FASTA_FILE | sed -r "s/\.fasta//"`.gff
            $PUTS_FASTA_TO_GFF3 --fasta_dir $WORK_DIR/puts/$PUTS_FASTA_FILE --source PlantGDB --type sequence_assembly --attributes "Dbxref=taxon:$TAXON_ID;organism=$SCIENTIFIC_NAME" --gfffilename $PUTS_GFF3_FILE
        fi
    done
    cd $WORK_DIR

    trap - INT TERM EXIT
    touch ".checkpoint_23"
else 
    echo "== PUTs GFF3 creation checkpoint complete" 2>&1
fi


if [ ! -e ".checkpoint_31" ]; then
    trap "echo 'FAILED while deploying files to FTP site'; " INT TERM EXIT
    echo "Transferring files to FTP site..."
    $SCP_FILES "$WORK_DIR" 1
    trap - INT TERM EXIT
    touch ".checkpoint_31"
else 
    echo "== Transfer of files to FTP site complete" 2>&1
fi

## populate the project fasta repository

## project fasta repository
FASTA_DIR_ROOT=/projects/solcomp/fasta_repository
FASTA_DIR="/projects/solcomp/fasta_repository/$TIMESTAMP"
if [ ! -e ".checkpoint_41" ]; then
    trap "echo 'FAILED while updating project fasta repository'; " INT TERM EXIT
    echo "Updating project fasta repository..."
    if [ -d $FASTA_DIR ]; then
       rm -rf $FASTA_DIR
    fi
    mkdir -p $FASTA_DIR 
    cp -rf $WORK_DIR/bacs $FASTA_DIR/
    cp -rf $WORK_DIR/bac_ends $FASTA_DIR/
    cp -rf $WORK_DIR/puts $FASTA_DIR/
    rm -rf $FASTA_DIR_ROOT/current
    chmod -R g-w $FASTA_DIR
    ln -s $FASTA_DIR $FASTA_DIR_ROOT/current
    trap - INT TERM EXIT
    touch ".checkpoint_41"
else 
    echo "== Update of project fasta repository complete" 2>&1
fi
