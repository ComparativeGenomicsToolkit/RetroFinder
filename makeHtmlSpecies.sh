#!/bin/bash
echo "dir = $3/$1"
echo `pwd`
sort -k5,5nr $2 > $2.sort ; mv $2.sort $2
tawk '{print $4}' $2 > tmp.select
echo "makeHtmlSpecies.sh bin $5 file $2 "
pwd
echo "grep -F -f tmp.select retroMrnaInfo.txt output $2.txt"
grep -F -f tmp.select retroMrnaInfo.txt > $2.txt
echo "$5/bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3/$1"
$5/bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3/$1
count=`wc -l $2 | awk '{print $1}'`
echo copy $count to /usr/local/apache/htdocs/retro/$3 
echo "<TD><a href=../$3/$1>"$count"</a></Td>" >> /usr/local/apache/htdocs/retro/$3/../index.html
