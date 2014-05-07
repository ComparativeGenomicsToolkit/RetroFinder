#!/bin/bash
set -beEu -x -o pipefail
db=$1
bed=$2
dir=$3
name=$4
script=$5
webdir=$6
server=$7
#set -o verbose
echo "Writing $name html into $dir/$db db = $db"
sort -k5,5nr $bed > $bed.sort ; mv $bed.sort $bed
$script/selectById 4 ${bed%%.bed}.recent.bed 1 retroMrnaInfo.txt > $bed.recent.txt
$script/selectById 4 ${bed%%.bed}.ancient.bed 1 retroMrnaInfo.txt > $bed.ancient.txt
mkdir -p $webdir/${dir}Anc/$db
$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $bed.recent.txt -label-tsv retroMrnaInfo.$db.lab -browser-url $server $db ${bed%%.bed}.recent.bed $webdir/${dir}Rec/$db
$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $bed.ancient.txt -label-tsv retroMrnaInfo.$db.lab -browser-url $server $db ${bed%%.bed}.ancient.bed $webdir/${dir}Anc/$db
countR=`wc -l ${bed%%.bed}.recent.bed | awk '{print $1}'`
countA=`wc -l ${bed%%.bed}.ancient.bed | awk '{print $1}'`
echo copy counts to $webdir/$dir/../index.html
echo "<TD><a href=../${dir}Rec/$db>"$countR"</a></Td>" >> $webdir/${dir}Rec/../index.html
echo "<TD><a href=../${dir}Anc/$db>"$countA"</a></Td>" >> $webdir/${dir}Anc/../index.html
#$split $bed ${bed%%.bed}.ancient.bed ${bed%%.bed}.recent.bed
#tawk '{print $4}' $bed > tmp.select
#echo "grep -F -f tmp.select retroMrnaInfo.txt output $bed.txt"
#grep -F -f tmp.select retroMrnaInfo.txt > $bed.txt
