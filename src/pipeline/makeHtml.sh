#
tawk '{print $4}' $2 > tmp.select
grep -F -f tmp.select retroMrnaInfo.txt |cut -f1,4 |sort > $2.txt
cut -f4,5,7-20 $2 |sort> $2.append.txt
join $2.txt $2.append.txt |awk '{OFS="\t";$1=$1;print $1,$0}'> $2.join.txt
sort -k4,4nr $2.join.txt > $2.final.txt
echo label $2.final.txt
echo copy to /usr/local/apache/htdocs/retro/$3
echo ../bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.final.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3
../bedToHtmlDir -dir-frame-per 75 -dir-right -context-bases 500 -labels $2.final.txt -label-tsv retroMrnaInfo.lab -browser-url http://hgwdev-baertsch.cse.ucsc.edu $1 $2 /usr/local/apache/htdocs/retro/$3
