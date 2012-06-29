#
echo reading gp.list, counting genes with one coding exon
for i in `cat gp.list` ; do  echo $i ; wc -l $i |awk '{print $1}' ; genePredFilter -cdsExons=2 $i stdout | wc -l ; done > x
awk 'BEGIN{n=0}n==2{printf("%d\t%d\t%4.2f\t%s\n", tot-$1,tot,(tot-$1)/tot,name);n++}n==1{tot=$1;n++}n==0{name= $1;n++}n==3{n=0}' x
