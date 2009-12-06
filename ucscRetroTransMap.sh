#
# get mrnas from transmap if genbank not setup for this genome or native mRNA are scarse
source DEF
cd $MRNABASE
hgsql $DB -N -B -e "select * from all_mrna" |cut -f2-22 > all_mrna.native.psl
hgsql $DB -N -B -e "select * from transMapAlnMRna" |cut -f2-22 > transMapAlnMRna.psl
hgsql $DB -N -B -e "select * from transMapAlnUcscGenes" |cut -f2-22 > transMapAlnUcscGenes.psl
hgsql $DB -N -B -e "select * from transMapAlnRefSeq" |cut -f2-22 > transMapAlnRefSeq.psl
sort -k10,10 transMapAlnMRna.psl >  transMapAlnMRna.qName.psl 
awk -f ~/baertsch/scripts/stripextversion.awk < transMapAlnMRna.qName.psl > transMapAlnMRna.strip.psl 
pslCDnaFilter -bestOverlap -maxAligns=1 transMapAlnMRna.strip.psl transMapAlnMRna.best.psl

hgsql hgFixed -N -B -e "select id, cds from transMapGeneMRna" > mrnaCds.tab
hgsql hgFixed -N -B -e "select id, cds from transMapGeneUcscGenes" > kgCds.tab
hgsql hgFixed -N -B -e "select id, cds from transMapGeneRefSeq" > refSeqCds.tab
cat *Cds.tab > cds.tab
gzip cds.tab &


mrnaToGene transMapAlnMRna.best.psl transMapAlnMRna.gp -cdsFile=mrnaCds.tab -ignoreUniqSuffix 2> error.mrna
mrnaToGene transMapAlnUcscGenes.psl transMapAlnUcscGenes.gp -cdsFile=kgCds.tab -ignoreUniqSuffix 2> error.ucsc
mrnaToGene transMapAlnRefSeq.psl transMapAlnRefSeq.gp -cdsFile=refSeqCds.tab -ignoreUniqSuffix 2> error.refseq

getRnaPred -cdsUpper panTro2 transMapAlnRefSeq.gp all refseqCds.fa ; tr acgt ACGT < refseqCds.fa > refseq.fa &
getRnaPred -cdsUpper panTro2 transMapAlnUcscGenes.gp all ucscCds.fa ; tr acgt ACGT < ucscCds.fa > ucsc.fa &
mkdir -p $OUTDIR/run.getmrna
mkdir -p $OUTDIR/mrna.psl
rm -f run.getmrna/jobList
for chr in `cut -f1 $GENOME/$DB/chrom.sizes`; do echo "$SCRIPT/genePredToFa.sh ../DEF $chr " >>$OUTDIR/run.getmrna/jobList ;done
ssh -T $CLUSTER "cd $OUTDIR/run.getmrna ; /parasol/bin/para -ram=4g make jobList"
cat run.getmrna/mrnaCds*.fa > mrnaCds.fa
cat run.getmrna/mrna.*.fa > mrna.fa
awk -f $SCRIPT/stripversion.awk < mrna.psl/*.psl > transmap_mrna.psl
cat all_mrna.native.psl transmap_mrna.psl transMapAlnRefSeq.psl > all_mrna.psl
gzip all_mrna.psl

wait
