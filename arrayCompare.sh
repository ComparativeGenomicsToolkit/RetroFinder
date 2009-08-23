#!/bin/bash
set -o verbose
source DEF
hgsql hgFixed -N -B -e "select h.id+1, m.id+1, m.name from $ARRAY1 h , $ARRAY2 m where h.name like m.name" |sed -e 's/ /_/g' > idmatch.tab
for i in $SPECIES; do F=`awk -F, '{print NF}' $BASE/$i/$INPUT | uniq` ; echo $i $F ; let F=F-1 ; cut -d , -f 2-$F $BASE/$i/$INPUT > retroTissue.$i.tab; done
#cut -d , -f 2-80 retroTissueExtract.hg18.bed> retroTissue.hg18.tab
#cut -d , -f 2-80 retroTissueExtract.mm9.bed> retroTissue.mm9.tab
echo 'source("'$SCRIPT'/overHist.R")'> plotArrayHist.R
echo 'h<-scan("retroTissue.hg18.tab",sep=",")' >> plotArrayHist.R
echo "dim(h)<-c(79,length(h)/79)" >> plotArrayHist.R
echo 'm<-scan("retroTissue.mm9.tab",sep=",")' >> plotArrayHist.R
echo "dim(m)<-c(61,length(m)/61)" >> plotArrayHist.R
tawk '{print "pdf(\"'$WEBDIR'/retro."$3".hg18vsMm9.pdf\")";print "overhist(h["$1",],m["$2",], breaks=100, xlab0=\""$3" (blue=human green=mouse)\")";print "dev.off()"}' idmatch.tab  >> plotArrayHist.R
R < plotArrayHist.R --no-save

cp $SCRIPT/header.html $WEBDIR/index.html
echo "<h2>Welch 2 sided T-test (assuming unequal variance) between human and mouse retrogenes as measured by gnfAtlas2 microarry data</h2>" >> $WEBDIR/index.html
echo "<h3> Histogram of the log ratio of all the retrogenes in that tissue</h3>" >> $WEBDIR/index.html
echo "<h3>  +5 highly expressed -5 repressed relative to the median of all tissues</h3>" >> $WEBDIR/index.html
 
echo 'h<-scan("retroTissue.hg18.tab",sep=",")' > plotArrayTtest.R
echo "dim(h)<-c(79,length(h)/79)" >> plotArrayTtest.R
echo 'm<-scan("retroTissue.mm9.tab",sep=",")' >> plotArrayTtest.R
echo "dim(m)<-c(61,length(m)/61)" >> plotArrayTtest.R
#tawk '{print "hm<-t.test(h["$1",],m["$2",], var.equal=FALSE) ; hm->p.value"}' idmatch.tab  >> plotArrayTtest.R
tawk '{print "ht<-t.test(h["$1",],m["$2",], var.equal=FALSE);df<-data.frame(pvalue=ht$p.value,humanEst=ht$conf.int[1],mouseEst=ht$conf.int[2]);df"}' idmatch.tab >>plotArrayTtest.R
#R< plotArrayTtest.R --slave --no-save 2> ttest.err >> tissueCompare.html
#R CMD BATCH plotArrayTtest.R 
R --slave < plotArrayTtest.R  |grep -v pvalue> plotArrayTtest.Rout

cut -f 3 idmatch.tab | paste plotArrayTtest.Rout - |awk '{print $2"\t"$3,$4,$5}' |sort -g | awk 'BEGIN{print "<tr><th>pvalue<th>human mean<th>mouse mean<th>distribution<th>human<th>mouse</tr>"}{print "<tr><td>"$1"</td><td>"$2"</td><td>"$3"</td><td><a href=\"retro."$4".hg18vsMm9.pdf\">"$4"</a></td><td><a href=\"tissue/hg18/"$4"/index.html\">all</a></td><td><a href=\"tissue/mm9/"$4"/index.html\">all</a></td></tr>"}' >> $WEBDIR/index.html

