#
set -beEu -o pipefail
source $1
bed=$2
dir=$3
name=$4
echo "mkdir -p $WEBROOT/$3"
mkdir -p $WEBROOT/$3
echo "<td>$4</td>" >> $WEBROOT/$dir/../index.html
mkdir -p ../../$DB/$EXPDIR ; pushd ../../$DB/$EXPDIR ; $SCRIPT/makeHtmlSpeciesOne.sh $DB $bed $dir $name $SCRIPT $WEBROOT ; popd 
echo "</TR>" >> $WEBROOT/$dir/../index.html
