hgsql hgFixed -N -B -e "select h.id+1, m.id+1, m.name from gnfHumanAtlas2MedianExps h , gnfMouseAtlas2MedianExps m where h.name like m.name" > idmatch.tab
BASE=/hive/users/baertsch/retro
for i in hg18 mm9 ; do F=`awk -F, '{print NF}' $BASE/$i/exp/retroTissueExtract.bed | uniq` ; echo $i $F ; let F=F-1 ; cut -d , -f 2-$F $BASE/$i/exp/retroTissueExtract.bed > retroTissue.$i.tab; done
cut -d , -f 2-80 retroTissueExtract.hg18.bed> retroTissue.hg18.tab
cut -d , -f 2-80 retroTissueExtract.mm9.bed> retroTissue.mm9.tab
echo 'h<-scan("retroTissue.hg18.tab",sep=",")' > plotArrayHist.R
echo "dim(h)<-c(79,length(h)/79)" >> plotArrayHist.R
echo 'm<-scan("retroTissue.mm9.tab",sep=",")' >> plotArrayHist.R
echo "dim(m)<-c(62,length(m)/62)" >> plotArrayHist.R
tawk '{print "pdf(\"pdf/retro."$3".hg18vsMm9.pdf\")";print "overhist(h["$1",],m["$2",], breaks=100, xlab0=\""$3"\")";print "dev.off()"}' idmatch.tab  >> plotArrayHist.R

tawk '{print "ttest["$3"]<- t.test(h["$1",],m["$2",], var.equal=TRUE)"}' idmatch.tab  > plotArrayTtest.R
