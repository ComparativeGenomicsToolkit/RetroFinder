#!/bin/bash
set -o verbose
db=$1
bed=$2
dir=$3
name=$4
script=$5
webdir=$6
echo "Writing $name html into $dir/$db "
sort -k5,5nr $bed > $bed.sort ; mv $bed.sort $bed
$script/selectById 4 ${bed} 1 retroMrnaInfo.txt > $bed.txt
mkdir -p $webdir/${dir}/$db
$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $bed.txt -label-tsv retroMrnaInfo.$db.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $db ${bed} $webdir/${dir}
countA=`wc -l ${bed}| awk '{print $1}'`
echo "<TD><a href=../${dir}/$db>"$countA"</a></Td>" >> $webdir/${dir}/../index.html
