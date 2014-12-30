#!/bin/bash
set -beEu -x -o pipefail
source $1
bed=$2
dir=$3
name=$4
if [[ -s $bed ]] ; then
    if [[ -a ${bed%%.bed}.recent.bed ]] ; then
        mkdir -p $WEBROOT/$3
        echo "<td>$4</td>" >> $WEBROOT/$dir/../index.html
        for sp in `echo $SPECIES` ; do mkdir -p $BASE/version$VERSION/$sp/$EXPDIR ; pushd $BASE/version$VERSION/$sp/$EXPDIR ; $SCRIPT/makeHtmlSpecies.sh $sp $bed $dir $name $SCRIPT $WEBROOT $WEBSERVER ; popd ; done
        echo "</TR>" >> $WEBROOT/$dir/../index.html
    else
        echo "NO ANCIENT OR RECENT BED FILES"
        mkdir -p $WEBROOT/$3/
        echo "<td>$4</td>" >> $WEBROOT/$dir/../index.html
        for sp in `echo $SPECIES` ; do mkdir -p $BASE/version$VERSION/$sp/$EXPDIR; pushd ../../$sp/$EXPDIR ; $SCRIPT/makeHtmlSpeciesOne.sh $sp $bed $dir $name $SCRIPT $WEBROOT $WEBSERVER ; popd ; done
        echo "</TR>" >> $WEBROOT/$dir/../index.html
    fi
fi
