#!/bin/bash
db=$1
bed=$2
dir=$3
name=$4
script=$5
webdir=$6
echo "Writing $name html into $dir/$db "
sort -k5,5nr $bed > $bed.sort ; mv $bed.sort $bed
#tawk '{print $4}' $bed > tmp.select
#echo "grep -F -f tmp.select retroMrnaInfo.txt output $bed.txt"
#grep -F -f tmp.select retroMrnaInfo.txt > $bed.txt
$script/selectById 4 $bed 1 retroMrnaInfo.txt > $bed.txt
echo "$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $bed.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $db $bed $webdir/$dir/$db"
$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $bed.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $db $bed $webdir/$dir/$db
count=`wc -l $bed | awk '{print $1}'`
echo copy $count to $webdir/$dir/../index.html
echo "<TD><a href=../$dir/$db>"$count"</a></Td>" >> $webdir/$dir/../index.html
