#
#calc divergence histograms 
#$1 = DEF file
source $1
#cd /hive/users/baertsch/retro/$DB
hgsql $DB -e "select (1000-millibad)/10, chrom, type from ${TABLE}" |awk '$2=="chrX"{print $1,$2,$3,$4,$5}$2!="chrX"{print $1,"auto",$2,$3,$4,$5}' | sort -n>diverge.txt
tawk '$1=="chrX"{c="chrX"}$1!="chrX"{c="auto"}{print (1000-$35)/10,c,$15}' exp/age/retroKgCoding.bed > divergeKg.txt
grep pseudogene diverge.txt |awk '{print $1,$2}' > div.pseudo.txt
grep pseudogene diverge.txt |awk '{print $1,$2}' > div.pseudoweak.txt
grep express diverge.txt |grep -v strong |grep -v shuffle |awk '{print $1,$2}' >> div.pseudoweak.txt
grep express diverge.txt |awk '{print $1,$2}' > div.express.txt
grep strong diverge.txt |awk '{print $1,$2}' > div.strong.txt
grep shuffle diverge.txt |awk '{print $1,$2}' >> div.strong.txt
cat divergeKg.txt |awk '{print $1,$2}' > div.strong.kg.txt
grep express diverge.txt|grep -v strong |awk '{print $1,$2}' > div.weak.txt

echo "pdf('div.${DB}.pseudo.pdf')" > div.R
echo "op <- par(mfcol=c(2, 2))" >> div.R
echo "pseudo<-read.table('div.pseudo.txt')">> div.R
echo "express<-read.table('div.express.txt')" >> div.R
echo "strong<-read.table('div.strong.txt')" >> div.R
echo "pseudoweak<-read.table('div.pseudoweak.txt')">> div.R
echo "strongKg<-read.table('div.strong.kg.txt')" >> div.R

echo "pseudogene.${DB}<-subset(pseudo, V1<=$MAXDIVERGENCE)" >> div.R
echo "length(pseudogene.${DB})">>div.R
echo "hp<-hist(pseudogene.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
#echo "hp<-hist(pseudogene.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R

echo "expressedRetro.${DB}<-express[express<=$MAXDIVERGENCE]" >> div.R
echo "length(expressedRetro.${DB})">>div.R
echo "he<-hist(expressedRetro.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}')" >> div.R
#echo "he<-hist(expressedRetro.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R

echo "par(op)">>div.R
echo "dev.off()">> div.R

echo "pdf('div.${DB}.strong.pdf')" >> div.R
echo "op <- par(mfcol=c(1, 3))" >> div.R
echo "pseudoweak.${DB}<-pseudoweak[pseudoweak<=$MAXDIVERGENCE]" >> div.R
echo "length(pseudoweak.${DB})">>div.R
echo "hp<-hist(pseudoweak.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM1),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
#echo "hp<-hist(pseudogene.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
echo "text(hp\$mids, hp\$density, hp\$counts, adj=c(.5, -.5), col='blue3')" >> div.R


echo "expressStrong.${DB}<-strong[strong<=$MAXDIVERGENCE]" >> div.R
echo "expressKg.${DB}<-strongKg[strongKg<=$MAXDIVERGENCE]" >> div.R
echo "length(expressStrong.${DB})">>div.R
echo "hs<-hist(expressStrong.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM2),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
#echo "hs<-hist(expressStrong.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
echo "text(hs\$mids, hs\$density, hs\$counts, adj=c(.5, -.5), col='green3')" >> div.R
echo "hk<-hist(expressKg.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM2),xlab='Range of substitution level (percent)',ylab='Fraction of ${GENOMENAME} known retrogene', main='RetroGenes ${GENOMENAME}',sub='each bar represents ~25MYA')" >> div.R
echo "text(hk\$mids, hk\$density, hk\$counts, adj=c(.5, -.5), col='red3')" >> div.R
echo "par(op)">>div.R
echo "dev.off()">> div.R

echo "q()" >> div.R


R --no-save <div.R
cp -f div.${DB}.strong.pdf ~/.html/pdf
cp -f div.${DB}.pseudo.pdf ~/.html/pdf
#cp -f div.${DB}.weak.pdf ~/.html/pdf
#cp -f div.${DB}.pseudoweak.pdf ~/.html/pdf
