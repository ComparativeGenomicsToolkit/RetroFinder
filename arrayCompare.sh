#!/bin/bash
set -o verbose
source DEF
hgsql hgFixed -N -B -e "select h.id+1, m.id+1, m.name from $ARRAY1 h , $ARRAY2 m where h.name like m.name" > idmatch.tab
for i in $SPECIES; do F=`awk -F, '{print NF}' $BASE/$i/$INPUT | uniq` ; echo $i $F ; let F=F-1 ; cut -d , -f 2-$F $BASE/$i/$INPUT > retroTissue.$i.tab; done
#cut -d , -f 2-80 retroTissueExtract.hg18.bed> retroTissue.hg18.tab
#cut -d , -f 2-80 retroTissueExtract.mm9.bed> retroTissue.mm9.tab
echo 'source("'$SCRIPT'/overHist.R")'> plotArrayHist.R
echo 'h<-scan("retroTissue.hg18.tab",sep=",")' >> plotArrayHist.R
echo "dim(h)<-c(79,length(h)/79)" >> plotArrayHist.R
echo 'm<-scan("retroTissue.mm9.tab",sep=",")' >> plotArrayHist.R
echo "dim(m)<-c(61,length(m)/61)" >> plotArrayHist.R
tawk '{print "pdf(\"'$WEBDIR'/retro."$3".hg18vsMm9.pdf\")";print "overhist(h["$1",],m["$2",], breaks=100, xlab0=\""$3"\")";print "dev.off()"}' idmatch.tab  >> plotArrayHist.R
R < plotArrayHist.R --no-save

cp $SCRIPT/header.html tissueCompare.html
#echo 'invisible(options(echo = FALSE))' > plotArrayTtest.R
echo 'h<-scan("retroTissue.hg18.tab",sep=",")' > plotArrayTtest.R
echo "dim(h)<-c(79,length(h)/79)" >> plotArrayTtest.R
echo 'm<-scan("retroTissue.mm9.tab",sep=",")' >> plotArrayTtest.R
echo "dim(m)<-c(61,length(m)/61)" >> plotArrayTtest.R
tawk '{print "t.test(h["$1",],m["$2",], var.equal=FALSE)$p.value"}' idmatch.tab  >> plotArrayTtest.R
R< plotArrayTtest.R --slave --no-save 2> ttest.err >> tissueCompare.html
#R CMD BATCH plotArrayTtest.R 


