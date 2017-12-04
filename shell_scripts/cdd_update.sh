#!/bin/bash

# This script checks that our installed current version of CDD is up to date.
# If it isn't, it will download and build up-to-date CDD databases for RPSBLAST.
#
# Brett Whitty
# bwhitty@jcvi.org

set -e

export CDD_BUILD_PATH="/home/whitty/projects/db/CDD"
export CDD_VERSION_URL="ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info"
export CDD_URL="ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz"

cd $CDD_BUILD_PATH

echo "Retrieving info on current CDD version"
wget -q $CDD_VERSION_URL
export VERSION=`perl -e '$_ = <>; chomp; /cdd version (.*)/; print $1;' ${CDD_BUILD_PATH}/cdd.info`
rm -f $CDD_BUILD_PATH/cdd.info
echo "Current CDD version is $VERSION"

if ! test -d $CDD_BUILD_PATH/cdd/$VERSION
then
    echo "CDD is out of date, fetching new build..."
    mkdir -p $CDD_BUILD_PATH/cdd/$VERSION
    mkdir -p $CDD_BUILD_PATH/curated/$VERSION
    mkdir $CDD_BUILD_PATH/tmp
    cd $CDD_BUILD_PATH/tmp
    trap "echo '...update aborted!'; cd $CDD_BUILD_PATH; rm -rf $CDD_BUILD_PATH/cdd/$VERSION; rm -rf $CDD_BUILD_PATH/curated/$VERSION; rm -rf $CDD_BUILD_PATH/tmp" INT TERM EXIT
    wget $CDD_URL
    echo "Extracting CDD data files..."
    tar xzf $CDD_BUILD_PATH/tmp/cdd.tar.gz
    rm -f $CDD_BUILD_PATH/tmp/cdd.tar.gz
    ls cd*.smp >Curated.pn
    echo "Building rpsblast cdd database..."
    formatrpsdb -t cdd -i Cdd.pn -o T -f 9.82 -n cdd -S 100.0
    mv -f $CDD_BUILD_PATH/tmp/cdd.* $CDD_BUILD_PATH/cdd/$VERSION
    ln -sf $CDD_BUILD_PATH/cdd/$VERSION $CDD_BUILD_PATH/cdd/current
    echo "Building rpsblast cdd_curated database..."
    formatrpsdb -t cdd_curated -i Curated.pn -o T -f 9.82 -n cdd_curated -S 100.0
    mv -f $CDD_BUILD_PATH/tmp/cdd_curated.* $CDD_BUILD_PATH/curated/$VERSION
    ln -sf $CDD_BUILD_PATH/curated/$VERSION $CDD_BUILD_PATH/curated/current
    cd $CDD_BUILD_PATH
    rm -rf $CDD_BUILD_PATH/tmp
    trap - INT TERM EXIT
    echo "CDD has been updated"
else 
    echo "CDD is up to date"
fi
