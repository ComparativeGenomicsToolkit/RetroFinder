#
set -beEu -o pipefail
source $1
mkdir -p $ROOTDIR/retro/$3
echo "<td>$4</td>" >> $ROOTDIR/retro/$3/../index.html
for db in `echo $SPECIES` ; do mkdir -p ../../$db/$EXPDIR ; pushd ../../$db/$EXPDIR ; $SCRIPT/makeHtmlSpecies.sh $db $2 $3 $4 $SCRIPT; popd ; done
echo "</TR>" >> $ROOTDIR/retro/$3/../index.html
