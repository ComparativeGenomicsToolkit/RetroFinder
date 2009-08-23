#
set -beEu -o pipefail
source $1
bed=$2
dir=$3
name=$4
echo "mkdir -p $ROOTDIR/retro/$3"
mkdir -p $ROOTDIR/retro/$3
echo "<td>$4</td>" >> $ROOTDIR/retro/$dir/../index.html
mkdir -p ../../$DB/$EXPDIR ; pushd ../../$DB/$EXPDIR ; $SCRIPT/makeHtmlSpeciesOne.sh $DB $bed $dir $name $SCRIPT $ROOTDIR/retro ; popd 
echo "</TR>" >> $ROOTDIR/retro/$dir/../index.html
