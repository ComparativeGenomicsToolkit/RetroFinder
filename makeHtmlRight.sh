#
sort -k5,5nr $2 > $2.sort ; mv $2.sort $2
tawk '{print $4}' $2 > tmp.select
#echo tawk '{print $4}' $2 to tmp.select
#echo grep -F -f tmp.select retroMrnaInfo.txt to  $2.txt
grep -F -f tmp.select retroMrnaInfo.txt > $2.txt
#~markd/bin/bedToHtmlDir -dir-frame-per 75 -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3
./bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3
count=`wc -l $2 | awk '{print $1}'`
echo copy $count to /usr/local/apache/htdocs/retro/$3 
echo "<TR><TH>"$4"<TD><a href=../$3>"$count"</a></TR>" >> /usr/local/apache/htdocs/retro/$3/../index.html
