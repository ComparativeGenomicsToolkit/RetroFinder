#
set -beEu -o pipefail
source $1
bed=$2
dir=$3
name=$4
sort -k5,5nr $2 > $2.sort ; mv $2.sort $2
tawk '{print $4}' $2 > tmp.select
grep -F -f tmp.select retroMrnaInfo.txt > $2.txt
bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url $WEBSERVER $1 $2 $WEBROOT/$3
count=`wc -l $2 | awk '{print $1}'`
echo copy $count to $WEBROOT/$3 
echo "<TR><TH>"$4"<TD><a href=../$3>"$count"</a></TR>" >> $WEBROOT/$3/../index.html
