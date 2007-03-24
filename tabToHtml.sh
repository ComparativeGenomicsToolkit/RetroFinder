tawk -f tabToHtml2.awk --assign outdir=/cluster/home/baertsch//.html/shuffleJuly/ shufflingAnalysis.txt 
#> hg18.2006-07-22.details.html
#cp hg18.2006-07-22.details.html ~/.html/shuffleJuly/
tawk -f tabToDir.awk < shufflingAnalysis.txt > hg18.2006-07-22.dir.html     
cp hg18.2006-07-22.dir.html ~/.html/shuffleJuly/

