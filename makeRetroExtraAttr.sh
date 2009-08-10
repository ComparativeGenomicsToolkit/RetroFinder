#
set -beEu -o pipefail
source $1
#DB=$1
#TABLE=$2
cd $2
echo "-------- script makeRetroExtraAttr.sh ------------"
echo "extract genbank names"
SQLNAME="select gb.acc, n.name from gbCdnaInfo gb, geneName n, organism o where gb.geneName = n.id and type = 'mRNA' and organism = o.id and o.name = '$GENOMENAME'"
echo "SQL=" $SQLNAME
hgsql $DB -N -B -e "$SQLNAME" > genBankName.txt
echo "drop rbGenBankName"
hgsql $DB -e "drop table rbGenBankName "
echo "create rbGenBankName"
hgsql $DB -e "create table rbGenBankName (acc varchar(255), name varchar(255), PRIMARY KEY(acc));"
echo "load rbGenBankName"
hgsql $DB -e "load data local infile 'genBankName.txt' into table rbGenBankName"

#hgsql $DB -N -B -e "select r.name, r.name, score, concat(f.geneName,',',replace(n.name,'n/a','')), r.type, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad from gbCdnaInfo cd, geneName n,  retroMrnaInfo r left outer join refFlat f on r.refSeq = f.name where cd.geneName = n.id and cd.acc = kgName " > retroMrnaInfo.txt
echo "extract retroMrnaInfo.txt"
echo "hgsql $DB -N -B -e select r.name, r.name, score, n.name,r.type, retroExonCount, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad from $TABLE r left outer join rbGenBankName n on r.refSeq = n.acc " 
hgsql $DB -N -B -e "select r.name, r.name, score, n.name,r.type, retroExonCount, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad from $TABLE  r left outer join rbGenBankName n on r.refSeq = n.acc " > retroMrnaInfo.txt
#hgsql $DB -N -B -e "select r.name, r.name, score, replace(n.name,'n/a',''), r.type, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad from $TABLE  r , rbGenBankName n where r.refSeq = n.acc " > retroMrnaInfo.txt
pwd
hgsql $DB -B -e "select name, name, score, refSeq as parent, type, retroExonCount as Exons, overlapRhesus as Rhesus , overlapMouse as Mus, overlapDog as Dog, overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad from $TABLE  limit 1" > retroMrnaInfo.lab
wc -l retroMrnaInfo.txt retroMrnaInfo.lab

hgsql $DB -B -N -e "select chrom, chromStart, chromEnd, name, score, strand from $TABLE  where overlapRhesus < 20 and overlapDog < 20 and overlapMouse < 20" > retroAncient.bed
wc -l retroAncient.bed
echo "-------- END script makeRetroExtraAttr.sh ------------"
