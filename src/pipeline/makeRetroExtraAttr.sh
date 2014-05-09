#!/bin/bash
set -beEu -o pipefail
source $1
#DB=$1
#TABLE=$2
#cd $2
DIR=$2
echo "-------- script makeRetroExtraAttr.sh ------------"
echo "Extract genbank names"
SQLNAME="select gb.acc, n.name , p.name from gbCdnaInfo gb, geneName n, organism o , productName p where gb.geneName = n.id and type = 'mRNA' and organism = o.id and o.name = '$GENOMENAME' and productName = p.id"
echo "SQL=" $SQLNAME
hgsql $DB -N -B -e "$SQLNAME" > $DIR/genBankName.txt
echo "Drop rbGenBankName"
hgsql $DB -e "drop table if exists rbGenBankName "
echo "Create rbGenBankName"
hgsql $DB -e "create table rbGenBankName (acc varchar(255), name varchar(255), product varchar(255), PRIMARY KEY(acc));"
echo "Load rbGenBankName"
hgsql $DB -e "load data local infile '$DIR/genBankName.txt' into table rbGenBankName"

echo "Extract $DIR/retroMrnaInfo.txt"
hgsql $DB -N -B -e "select r.name, r.name, r.score, n.name,r.type, r.retroExonCount, -777, o1.overlap,-1788, r.overName, r.blockCount, r.conservedSpliceSites, r.exonCover, r.coverage, r.milliBad, n.product from $TABLE  r left outer join (rbGenBankName n, $ORTHOTABLE o1) on r.refSeq = n.acc and r.name = o1.name and o1.db = "\"$NET1\""" > $DIR/retroMrnaInfo.txt
#hgsql $DB -N -B -e "select r.name, r.name, score, n.name,r.type, retroExonCount, overlapRhesus, overlapMouse, overlapDog, overName, r.blockCount, conservedSpliceSites, exonCover, coverage, milliBad, n.product from $TABLE  r left outer join rbGenBankName n on r.refSeq = n.acc " > retroMrnaInfo.txt
pwd
hgsql $DB -N -B -e "select r.name as name , r.name as name, score, refSeq as parent, type, retroExonCount as Exons, -777, o1.overlap as "$NET1",-1788, overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad, product from $TABLE  r left outer join (rbGenBankName n, $ORTHOTABLE o1) on r.refSeq = n.acc and r.name = o1.name and o1.db = "\"$NET1\"" limit 1" > $DIR/retroMrnaInfo.$DB.lab
#hgsql $DB -B -e "select name, name, score, refSeq as parent, type, retroExonCount as Exons, overlapRhesus as "$NET3" , overlapMouse as "$NET1", overlapDog as "$NET2", overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad, name as product from $TABLE  limit 1" > retroMrnaInfo.$DB.lab
#cp retroMrnaInfo.hg18.lab retroMrnaInfo.$DB.lab
#hgsql $DB -B -e "select name, name, score, refSeq as parent, type, retroExonCount as Exons, overlapRhesus as Rat , overlapMouse as Human, overlapDog as Dog, overName as exp, blockCount, conservedSpliceSites as consSS, exonCover, coverage, milliBad , name as product from $TABLE  limit 1" > retroMrnaInfo.mm9.lab

echo "-------- END script makeRetroExtraAttr.sh ------------"
