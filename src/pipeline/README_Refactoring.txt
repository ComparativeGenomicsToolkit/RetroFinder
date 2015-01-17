2015-01-16
Inputs and Outputs:

pslPseudo is the main program that generates the final list of retrogenes and 
associated information. 

pslPseudo is run by the doBuildpk.sh script which is called in ucscRetroStep3.sh 
OUTDIR=/hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/retro/version${VERSION}/${DB}
Inputs to pslPseudo:
-cdsFile=cds.tab.gz - made in ucscRetroStep1.sh, acc, version, name, type from 
refSeqs and mRNAs (or from whatever sequences used as input). name field is the 
CDS in the cds table (one of the genbank tables)
database name - defined as DB in the DEF parameters file 
split/tmp<id>.psl - 
assembly chrom sizes - from browser, usually in /hive/data/genomes/<db>
repeatMasker, rmsk.bed.gz 
net1.txt.gz
net2.txt.gz
simpleRepeat.bed.gz
all_mrna.psl.gz
mrna<id>.psl
pseudo<id>.psl (output?)
<db>.2bit
mrna.2bit
refGene.tab.gz 
gencode.tab.gz (or ensGene.tab.gz)
knownGene.tab.gz
net3.txt.gz
ortho<id>.txt (output?)


Output:
pseudoMrnaLink<id>.txt
pseudo<id>.axt 

