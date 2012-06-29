#
sort -k4 $2 > $2.sort 
tawk '{print $4}' $2 > tmp.select
echo "grep -F -f tmp.select retroMrnaInfo.txt |sort > $2.txt"
grep -F -f tmp.select retroMrnaInfo.txt |sort > $2.txt
echo "join -1 4 $2.sort $2.txt > $2.join"
join -t '	' -1 4 -o 1.1,1.2,1.3,1.4,1.5,2.4,2.5 $2.sort $2.txt |sort -k6,6 -k4,4> $2.join
echo "sort -k4,4 -k1,1 $2.txt > $2.sort.txt"
sort -k4,4 -k1,1 $2.txt > $2.sort.txt 
echo "mv $2.sort.txt $2.txt"
mv $2.sort.txt $2.txt
wc -l $2.sort $2.txt
#echo "mv $2.join $2"
#v $2.join $2
##echo "cut -f1-57 $2.join > $2.sort"
##cut -f2-61 $2.join > $2.sort
echo "./bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2.join /usr/local/apache/htdocs/retro/$3"
./bedToHtmlDir  -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2.join /usr/local/apache/htdocs/retro/$3
count=`wc -l $2 | awk '{print $1}'`
echo $count copy to /usr/local/apache/htdocs/retro/$3
echo "<TR><TH>"$4"<TD><a href=../$3>"$count"</a></TR>" >> /usr/local/apache/htdocs/retro/$3/../index.html
