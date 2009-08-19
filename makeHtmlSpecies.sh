#!/bin/bash
db=$1
bed=$2
dir=$3
name=$4
script=$5
webdir=$6
echo "writing html into $dir/$db"
sort -k5,5nr $2 > $2.sort ; mv $2.sort $2
tawk '{print $4}' $2 > tmp.select
echo "makeHtmlSpecies.sh bin $5 file $2 "
pwd
echo "grep -F -f tmp.select retroMrnaInfo.txt output $2.txt"
grep -F -f tmp.select retroMrnaInfo.txt > $2.txt
echo "$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 $webdir/$3/$1"
$script/bedToHtmlDir -page-size 21 -dir-frame-per 60 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 $webdir/$3/$1
count=`wc -l $2 | awk '{print $1}'`
echo copy $count to $webdir/$3/../index.html
echo "<TD><a href=../$3/$1>"$count"</a></Td>" >> $webdir/$3/../index.html
