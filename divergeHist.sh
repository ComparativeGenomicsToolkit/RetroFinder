#
#calc divergence histograms 
#$1 = DEF file
source $1
cd /hive/users/baertsch/retro/$DB
hgsql $DB -e "select (1000-millibad)/10, type from ucscRetroInfo" |sort -n>diverge.txt
grep pseudogene diverge.txt |awk '{print $1}' > div.pseudo.txt
grep pseudogene diverge.txt |awk '{print $1}' > div.pseudoweak.txt
grep express diverge.txt |grep -v strong |grep -v shuffle |awk '{print $1}' >> div.pseudoweak.txt
grep express diverge.txt |awk '{print $1}' > div.express.txt
grep strong diverge.txt |awk '{print $1}' > div.strong.txt
grep shuffle diverge.txt |awk '{print $1}' >> div.strong.txt
grep express diverge.txt|grep -v strong |awk '{print $1}' > div.weak.txt

echo "pseudo<-scan('div.pseudo.txt')"> div.R
echo "express<-scan('div.express.txt')" >> div.R
echo "strong<-scan('div.strong.txt')" >> div.R
echo "pseudoweak<-scan('div.pseudoweak.txt')">> div.R

echo "expressedRetro.${DB}<-express[express<=$MAXDIVERGENCE]" >> div.R
echo "length(expressedRetro.${DB})">>div.R
echo "pdf('div.${DB}.express.pdf')" >> div.R
echo "he<-hist(expressedRetro.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}')" >> div.R
#echo "he<-hist(expressedRetro.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}')" >> div.R
echo "dev.off()">> div.R

echo "pseudogene.${DB}<-pseudo[pseudo<=$MAXDIVERGENCE]" >> div.R
echo "length(pseudogene.${DB})">>div.R
echo "pdf('div.${DB}.pseudo.pdf')" >> div.R
echo "hp<-hist(pseudogene.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}')" >> div.R
#echo "hp<-hist(pseudogene.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}')" >> div.R
echo "dev.off()">> div.R

echo "pseudoweak.${DB}<-pseudoweak[pseudoweak<=$MAXDIVERGENCE]" >> div.R
echo "length(pseudoweak.${DB})">>div.R
echo "pdf('div.${DB}.pseudoweak.pdf')" >> div.R
echo "hp<-hist(pseudoweak.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}')" >> div.R
#echo "hp<-hist(pseudogene.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} non-expressed retrocopies', main='Not Expressed ${GENOMENAME}')" >> div.R
echo "dev.off()">> div.R

echo "expressStrong.${DB}<-strong[strong<=$MAXDIVERGENCE]" >> div.R
echo "length(expressStrong.${DB})">>div.R
echo "pdf('div.${DB}.strong.pdf')" >> div.R
echo "hs<-hist(expressStrong.${DB},freq=TRUE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}')" >> div.R
#echo "hs<-hist(expressStrong.${DB},freq=FALSE,breaks=c(${BREAKS}), xlim=range(0,$XLIM),ylim=range(0,$YLIM),xlab='Range of substitution level, representing ~25 MYA each',ylab='Fraction of ${GENOMENAME} expressed retrocopies', main='Expressed ${GENOMENAME}')" >> div.R
echo "dev.off()">> div.R

echo "q()" >> div.R


R --no-save <div.R
cp -f div.${DB}.express.pdf ~/.html/pdf
cp -f div.${DB}.strong.pdf ~/.html/pdf
cp -f div.${DB}.pseudo.pdf ~/.html/pdf
cp -f div.${DB}.pseudoweak.pdf ~/.html/pdf
