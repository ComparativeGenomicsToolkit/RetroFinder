#
set -beEu -o pipefail
source $1
bed=$2
dir=$3
name=$4
if [[ -s $bed ]] ; then
    if [[ -a ${bed%%.bed}.recent.bed ]] ; then
        mkdir -p $ROOTDIR/retro/$3
        echo "<td>$4</td>" >> $ROOTDIR/retro/$dir/../index.html
        for db in `echo $SPECIES` ; do mkdir -p ../../$db/$EXPDIR ; pushd ../../$db/$EXPDIR ; $SCRIPT/makeHtmlSpecies.sh $db $bed $dir $name $SCRIPT $ROOTDIR/retro $WEBSERVER ; popd ; done
        echo "</TR>" >> $ROOTDIR/retro/$dir/../index.html
    else
        mkdir -p $ROOTDIR/retro/$3
        echo "<td>$4</td>" >> $ROOTDIR/retro/$dir/../index.html
        for db in `echo $SPECIES` ; do mkdir -p ../../$db/$EXPDIR ; pushd ../../$db/$EXPDIR ; $SCRIPT/makeHtmlSpeciesOne.sh $db $bed $dir $name $SCRIPT $ROOTDIR/retro $WEBSERVER ; popd ; done
        echo "</TR>" >> $ROOTDIR/retro/$dir/../index.html
    fi
fi
