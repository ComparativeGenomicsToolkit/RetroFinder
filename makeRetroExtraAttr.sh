#!/bin/bash
set -beEu -o pipefail
source $1
#DB=$1
#TABLE=$2
cd $2
echo "-------- script makeRetroExtraAttr.sh ------------"
echo "extract genbank names"
SQLNAME="select gb.acc, n.name , p.name from gbCdnaInfo gb, geneName n, organism o , productName p where gb.geneName = n.id and type = 'mRNA' and organism = o.id and o.name = '$GENOMENAME' and productName = p.id"
echo "SQL=" $SQLNAME
hgsql $DB -N -B -e "$SQLNAME" > genBankName.txt
echo "drop rbGenBankName"
hgsql $DB -e "drop table if exists rbGenBankName "
echo "create rbGenBankName"
hgsql $DB -e "create table rbGenBankName (acc varchar(255), name varchar(255), product varchar(255), PRIMARY KEY(acc));"
echo "load rbGenBankName"
hgsql $DB -e "load data local infile 'genBankName.txt' into table rbGenBankName"

echo "extract retroMrnaInfo.txt"
hgsql $DB -N -B -e "select r.name, r.name, score, n.name,r.type, retroExonCount, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad, n.product from $TABLE  r left outer join rbGenBankName n on r.refSeq = n.acc " > retroMrnaInfo.txt
pwd
hgsql $DB -B -e "select name, name, score, refSeq as parent, type, retroExonCount as Exons, overlapRhesus as Rhesus , overlapMouse as Mus, overlapDog as Dog, overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad, name as product from $TABLE  limit 1" > retroMrnaInfo.$DB.lab
#cp retroMrnaInfo.hg18.lab retroMrnaInfo.$DB.lab
hgsql $DB -B -e "select name, name, score, refSeq as parent, type, retroExonCount as Exons, overlapRhesus as Rat , overlapMouse as Human, overlapDog as Dog, overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad , name as product from $TABLE  limit 1" > retroMrnaInfo.mm9.lab
wc -l retroMrnaInfo.txt retroMrnaInfo.mm9.lab retroMrnaInfo.$DB.lab

echo "-------- END script makeRetroExtraAttr.sh ------------"
