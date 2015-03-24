Instructions for running the RetroFinder pipeline
RetroFinder predicts retrogenes including processed pseudogenes. 

1.     Data Prerequisites
In order to run the retroFinder pipeline, the following data files are 
mandatory: 1) set of mRNA sequences (actual or predicted); 2) repeatMasker 
annotation; 3) simpleRepeat annotation; 4) list of chromosome sizes; 5) 
predicted gene set and alignments; 6) net annotation for three other genomes - 
these should be from fairly complete genomes that are evolutionarily close 
enough to be alignable at the DNA level, the second being the farthest and the 
third the closest; 7) twoBit format file containing the genome sequence. Optionally, we 
can provide two additional gene sets to improve the selection of parent genes.  

The pipeline requires mRNA and/or annotated transcript sequences from the parent genes as the starting point. Including UTR sequence for the input sequences in mRNA.fa 
will improve the retrogene prediction. Most processed pseudogene prediction methods use proteins so they do not include UTR. However, genes under positive 
selection are sometimes missed and could be better detected with protein alignments. 

If the number of mRNAs and gene annotations for the genome are sparse, then use the script 
ucscTransMap.sh to map over predicted mRNAs from another genome that has more transcript data in the UCSC TransMap tracks into the mrna.fa file. Similarly, mRNAs for annotated genes should be stored in a file called refseq.fa. If annotated data is not available, 
simply create an empty file call refseq.fa to avoid script errors. Variables 
GENE1-3 should be set to mySQL tables or flat files containing genePred 
annotations (see UCSC Genome Browser FAQ for format information:
http://www.genome.ucsc.edu/FAQ/FAQformat.html) for three separate gene 
annotations. If only one is available, then set all three variables to 
the same table or there are two then one of them should be repeated.

The RetroFinder pipeline obtains the required data from a combination of flat files or from a mySQL database that uses UCSC data formats 
(http://www.genome.ucsc.edu/FAQ/FAQformat.html). The input parent gene annotation is obtained either from GenBank the user can provide TransMap annotations any other source of mRNA sequences. If mySQL is not available, the scripts can be easily modified to read 
from flat files instead.

2.     Source code and bash scripts
Source code and scripts are in an SVN repository and can be obtained using this command: 
svn co svn+ssh://riverdance.soe.ucsc.edu/projects/compbiousr/svnroot/hausslerlab/retroFinder. The scripts are in branches/version2/src/pipeline directory and 
the C code programs are under retroFinder/branches/version2/src. There 
are dependancies on the kent src tree (See http://genomewiki.ucsc.edu/index.php/CVS_to_Git_Migration and http://genome.ucsc.edu/FAQ/FAQdownloads.html#download27) which should be installed first. The location of the kent/src/ directory should be set in the 
first line of retroFinder/branches/version2/common.mk. Set environment variables USE_BAM=0 
and USE_SSL=0 before compiling the kent src libraries. 
cd retroFinder/branches/version2/src
make clean
make all
The executables are in the retroFinder/branches/version2/bin directory.

3.     Setting up the Parameter file
A default parameter file (DEF) is located in the scripts directory and this 
should be passed as the first parameter to all of the step scripts in the 
pipeline. Example DEF files are in the 
retroFinder/branches/version2/src/pipeline directory and are named as
DEF.<DB> where DB is the UCSC assembly database name e.g. hg38, mm10.
Variables used in the DEF file are documented in the separate README_DEF.txt 
file.

4. Setup before running the pipeline
Set up the working directory:
mkdir -p /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE 
where date is in format 
YYYYMMDD or the same format as DATE variable in the DEF file.
cd /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE
# Create the DEF file in this directory
chmod +x DEF
mkdir -p /hive/data/genomes/$DB/retro
cp DEF /hive/data/genomes/$DB/retro


# Create a tab-separated file of chromosome names and sizes. Remove the mitochondrion genome (chrM) and random chromosomes and any other sequences that are not to be used.
e.g.
cat /hive/data/genomes/$DB/chrom.sizes | grep -v random \
   | grep -v chrUn | grep -v chrM > S1.len
# Also copy to MRNABASE (defined in the DEF file):
mkdir -p /hive/data/genomes/$DB/bed/mrnaBlastz.$VERSION
cp S1.len /hive/data/genomes/$DB/bed/mrnaBlastz.$VERSION


cd /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/mrnaBlastz
# Run scripts from this directory


Pipeline steps: 
NOTE: Only steps ucscRetroStep[1-5].sh, filterEst.sh, filterMrna.sh and analyseExpress.sh need to be run to generate data for the UCSC Genome Browser track and for GENCODE. 
The ucscRetroStep6.sh script is optional and can be run if extra data is required for analysis. 
5.i.   Step 1 - Generate data sets for mRNAs, gene annotation, and expression 
data. Script: ucscRetroStep1.sh

The first step in the pipeline, by default, requires a set of mRNAs and RefSeq annotated transcripts to align to the genome in order to select potential retrogenes and their parents. For genomes with sufficient sequenced mRNAs, these sequences can be extracted from GenBank using the UCSC GenBank pipeline (gbGetSeqs). If that program is not available, 
extract all mRNAs from GenBank manually and put them in a file called 
mrna.fa in the $MRNABASE directory, as defined in the DEF parameter file. Any 
set of mRNAs or curated transcripts can be used as input if they are supplied 
in a file called mrna.fa.

If insufficient mRNAs are available from GenBank, then mRNAs from gene 
predictions can be used or TransMap mRNA can be used instead. The script 
ucscTransMap.sh will take TransMap mRNA data from tables (stored in a MYSQL 
database) and generate an mrna.fa file. 

A cds.tab.gz file must also be created for the input sequences. This file is tab-separated and has four columns: sequence id, version number, cds region (start..end or n/a if there is no CDS defined), sequence type (mRNA). 

If GenBank mRNAs and RefSeqs are to be used as input sequence data then this is the default and the script will fetch this data. If other sequences are to be
used the files must be prepared as described above before running the script. 


If Ensembl transcripts are to be used as input because the RefSeq and/or mRNA data is sparse then in the DEF file, set ALTSEQSOURCE to “ensGene” and set USEALTSEQS to 1. 
ucscRetroStep1.sh then runs a script, getEnsemblTxSeqsWithVersions.pl, to get the Ensembl sequences. 

To run the step 1 script:
retroFinder/branches/version2/src/pipeline/ucscRetroStep1.sh DEF >& step1.log &

Summary of step1: 
The script fetches mRNAs and RefSeq sequences from GenBank using gbGetSeqs. 
The PSL alignments and the CDS regions for these sequences are obtained from 
the appropriate genome assembly database.
A PSL file of only GenBank mRNAs is obtained for a later step. Due to frequent
updates of GenBank sequence alignments by the UCSC Genome Browser team, 
sequence versions and alignments can be out of sync if they are not 
all obtained at the same time. 
PolyA tails are removed from the input sequences and a lift file is created
so that the polyA tails can be added back after the lastz alignments. 
A script, lastz.sh, is called to align the trimmed input sequences to the 
genome using LASTZ. 
The mRNA sequences are in the MRNABASE directory (defined in the DEF file) and they are split up into batches for alignment into the TMPMRNA/split directory (TMPMRNA is defined in the DEF file and it is a working directory for alignments).


LASTZ:
Each sequence in the file mrna.fa against the genome sequence using LASTZ, a tool for aligning two DNA sequences while inferring scoring parameters automatically 
(http://www.bx.psu.edu/miller_lab/). 
Options used:
--ambiguous=iupac - this is used as sequences sometimes contain non-ATGC characters and output format is axt. 
--hspthresh=2000 - set threshold for high scoring pairs (default is 3000), ungapped extensions scoring lower are discarded
--format=axt -  output format is AXT
lastz.sh converts axt format output to PSL using the axtToPsl program. 


This parameter should be tuned depending on 
the evolutionary distance of the pair of species. WHAT DOES THIS MEAN? 
6.     Step 2 - Alignment 
The script, ucscRetroStep2.sh, firstly removes duplicate PSL rows. Resulting alignments are chained together with axtChain  (Kent WJ et. al, 2003). The concatenated chained alignments are sorted then filtered to keep those with a minimum percent identity of 58% and coverage of at least 5% using the utility, pslCDnaFilter in order to select those alignments that will enter the scoring step. The percent identity parameter was set to 58% using a standard ROC curve using HAVANA’s (Wellcome Trust Sanger Institute) manually annotated gene set, Vega (http://vega.sanger.ac.uk/index.html), as the training set.  
Filtered alignments are split into smaller files in the OUTDIR/split directory. pslQueryUniq makes each query name unique by adding a suffix e.g. -1, -2 etc. 
The input mRNA sequences are loaded into the databse: 
* mrna.fa is copied to /gbdb/$DB/blastzRetro${VERSION}
* refseq.fa is copied to /gbdb/$DB/blastzRetro${VERSION}
* hgLoadSeq used to load these files into ucscRetroSeq${VERSION} and ucscRetroExtFile${VERSION} is created to store the file name and its location


To run the step2 script:
retroFinder/branches/version2/src/pipeline/ucscRetroStep2.sh DEF >& step2.log &



7.     Step 3 - Scoring 
The scoring step takes the alignments created in the previous step and runs a cluster job to score the features and create a non-overlapping set of retrogene predictions. 
Input to the scoring program includes a number of annotation that are extracted by this script: RepeatMasker; simple repeats, gene annotations (GENE{1,2,3}); net alignments from three species, two far in evolutionary distance and one near (NET{1-3}). Annotation files are put in the OUTDIR (defined in the DEF file). The mrna.2bit, all_mrna.psl and cds.tab files are all symlinked to the OUTDIR.
The scoring program, pslPseudo, is run by the doBuild.pk script and it has a number of options  that can be set in the parameter file using the RETRO_OPTIONS variable in the DEF file.  


Options for pslPseudo:
 -minAli - the minimum alignment score, this is set so that by default the parent gene must have at least 98% identity with the genome
-nearTop - this allows for parent genes with recent duplications to be included in the pipeline.  It is set to 0.005 so that any mRNA alignment that is within 0.5% percent identity of the parent gene also is included as a parent.


pslPseudo Inputs: 
-cdsFile=cds.tab.gz - made in ucscRetroStep1.sh, acc, version, name, type from
refSeqs and mRNAs (or from whatever sequences used as input). The name field is the
CDS in the cds table (one of the GenBank tables in the UCSC databases)
DB - database name defined in the DEF parameters file
split/tmp<id>.psl - lastz alignments of input mRNAs and RefSeqs (or other input sequences)
chrom.sizes - genome assembly chromosome sizes from the UCSC genome browser, usually in /hive/data/genomes/<db>
rmsk.bed.gz - RepeatMasker repeats annotation for genome assembly, DB
net1.txt.gz - ucscRetroStep3.sh, query to NET1 (defined in DEF file) table in DB
net2.txt.gz - ucscRetroStep3.sh, query to NET2 (defined in DEF file) table in DB
simpleRepeat.bed.gz - ucscRetroStep3, query to simpleRepeat table in DB
all_mrna.psl.gz - BLAT PSL alignments for mRNAs and RefSeqs (or other input mRNA sequences) to the genome assembly, DB, obtained in script, ucscRetroStep1.sh 


pslPseudo Outputs: 
mrna<id>.psl - best alignments with and without introns, in OUT directory (OUTDIR/out)
pseudo<id>.psl - retrogene/pseudogene PSL alignments, in OUT directory (OUTDIR/out)
pseudoMrnaLink<id>.txt - stores link between pseudogene and parent gene, in TMP directory defined in doBuild.pk, moved later by this script
pseudo<id>.axt - axt output of retrogenes/pseudogenes, in TMP directory defined in doBuild.pk, moved later by this script
ortho<id>.txt - break in orthology information for ucscRetroOrtho table, in OUT directory (OUTDIR/out)


doBuildpk.sh script:
The script runs pslPseudo and the output pseudoMrnaLink<id>.txt is filtered on field 32 so only those lines where it is set to -1 or -2 are kept and the output is printed to pseudoMrnaLink<id>.bed. This selects pseudogenes and expressed retrogenes. Field 32 is the label: 1=pseudogene, -1=not pseudogene, -2=expressed retrogene.
OUTDIR=/hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/retro/version${VERSION}/${DB}


Inputs:
id - this is the id used to identify the split input PSL alignment files
DEF - file of variable definitions
pseudoGeneLink<id>.bed - name of output file, stores link between pseudogene and parent 


Outputs:  
pseudoMrnaLink<id>.bed is put in the OUTDIR/result directory 
pseudo<id>.axt is put in the OUTDIR/result/axt directory
pseudo<id>.log is put in the OUTDIR/log directory
Other outputs - see those for pslPseudo


NOTE: The variable PDB in the parameter file defines the name of the UCSC MYSQL protein database. This is not used in this script, only analyseExpress.sh

To run the step 3 script:
retroFinder/branches/version2/src/pipeline/ucscRetroStep3.sh DEF >& step3.log &

8.     Step 4 - Consolidation
This script (ucscRetroStep4.sh) concatenates the results (pseudoGeneLink<id>.bed files from the cluster job and to find the most likely parent of each putative retrogene, it selects those rows where the retrogene has a score > 300 and an axtScore > 10,000 and after removing tandem duplications (removeTandemDups script) and sorting, the output is in pseudoGeneLinkSortFilter.bed in the OUTDIR. The original parent can still be identified for retrogenes that have been tandemly duplicated later. 
removeTandemDups removes retrogenes that have parents on the same chromosome that are within 200k bp of the retrogene to remove false positives due to tandemly duplicated
gene families.
pseudo<id>.psl files are concatenated to pseudo.psl (used in ucscRetroStep5.sh) in OUTDIR
ortho<id>.txt files are concatenated removing “.txt” removed from the name of the net files in these ortho files and the resulting file is ortho.txt in OUTDIR. 
The pseudoGeneLinkSortFilter.bed file is split by chromosome and each file is input to bedOverlap which is run on the cluster to remove overlapping predictions, the highest scoring BED being the one that is kept in each case.


To run the step 4 script:
retroFinder/branches/version2/src/pipeline/ucscRetroStep4.sh DEF >& step4.log &

9.    Step 5 - Analysis and Loading Database Tables
This step (ucscRetroStep5.sh) is used to determine the level of expression 
of retrocopies. Counts of the number of transcripts that overlap the retrocopy 
are noted. To help distinguish protein-coding genes from transcribed 
pseudogenes, an ORF score is calculated (using bestORF) using a window of 250 
bases around the retrocopy.  Optionally, immunoglobin domains and zinc fingers 
are filtered based on the PFAM domain from the parent gene annotation. 

Alignments of mRNA and EST to the genome of interest are used for expression 
analysis and those are normally extracted from the the UCSC mySQL database as BLAT alignments. Alternatively, the user can provide zipped PSL files called all_mrna.psl and estfilter.psl. In order to resolve the best match for ESTs and mRNAs that align to multiple locations in the genome, the pipeline performs a stricter scoring of these alignments using a 
utility called pslCDnaGenomeMatch. It works by counting the number of alignments positions that differ between the two loci. If more than four positions in one loci score better than the other, that loci is picked as the location of the transcript. If there are fewer than four diagnostic positions, then the transcript is discarded.

Each of the predicted retrogenes is classified as one of two types. Type one is exon shuffling or chimeric retrogene, which is defined as a retrogene overlapping at least one exon of a multi-exon gene. Type two is the set of retrogenes that are inserted into locations that previously had no gene. Next, the retros are aged; ancient retrogenes are defined as those that do not have a break in orthology as defined by the UCSC net alignments. Recent retrogenes are ones that were inserted after the speciation event. The main feature table is called ucscRetroInfo (extended BED format) while ucscRetroAli contains the alignment to the 
parent gene. The version number (VERSION) in the DEF parameter file is appended to the table names.

NOTE: Only steps 1-5 are required to generate data for Genome Browser tracks. 
Step6 is an optional extra to generate additional data if it is required 
for analysis.

10.    Step 6 - Web pages for displaying analysis results (optional step)
The final step, ucscRetroStep6.sh, generates web pages for both duplicated 
retrocopies and chimeric retrogenes. It categorizes the retrocopies by age 
using the breaks in orthology feature (see section 2.2.4). Comparing other 
datasets can be also quickly generated with the compareRetroset.sh scripts (see 
section 4 or http://hgwdev.cse.ucsc.edu/~baertsch/retro for example web pages). 
The variable ROOTDIR in the DEF file defines the directory where the HTML 
pages will be created.  The SPECIES parameter in the DEF defines which web 
pages will be generated; entries should contain mysql database names of the 
genomes to be displayed.  Counts are displayed on summary pages with click 
through to each individual example.

11.     Scripts to be after running the  pipeline

filterMrna.sh and filterEst.sh extract both mRNAs and ESTs alignments in PSL 
format from the mySQL database for expression analysis. The parameter file 
contains variable settings to define whether or not the alignment tables are
split by chromosome (recent assemblies are not). 
The prerequisite data can be downloaded from the UCSC Genome 
Browser website (See http://genome.ucsc.edu/FAQ/FAQdownloads.html#download1). 
If there are no MYSQL tables for this data, those steps can be skipped and 
gzipped flat files of the mRNA and EST alignment data in PSL format can be 
utilised by the scripts instead. 

The scripts, filterMrna.sh and filterEst.sh, create sets of filtered BLAT 
alignments of mRNAs and ESTs respectively by picking the best hits. The 
analyseExpress.sh script then finds expressed retrogenes. Those retrogenes 
with at least 5 supporting ESTs and 1 spliced mRNA that do not overlap an 
annotated gene (multi-exonic genes from GENE1 and GENE2, typically UCSC Genes 
(knownGene table) and RefSeqs (refGene table) and have no shuffle events are 
selected.
orfBatch is used to find an ORF and those with a good ORF (as found by 
gene-check) are selected. The resulting genePred is loaded into the 
ucscRetroExpressed$VERSION table.

RETROFINDER DATABASE TABLES:
The following tables are loaded by the RetroFinder pipeline and each are
suffixed with $VERSION as defined in the DEF file:

ucscRetroAli - PSL lastz alignments of the predicted retrogenes

ucscRetroCds - CDS regions for parent mRNAs from GenBank (or other source if
input mRNAs are from another source). NOTE: not all GenBank mRNAs have a CDS
defined and some will be wrong.  

ucscRetroCount - contains the parent gene symbol with >1 retrogene and the 
count of the number of retrogenes for that parent gene, ordered by count in 
ascending order. 

ucscRetroExpressed - genePred format of retrogenes with evidence of expression
 that do not overlap annotated multi-exonic genes and have no evidence of
exon shuffling. 

ucscRetroSeq - sequences tables listing the input mRNA sequence ids, a number
referring to the sequences file in the ExtFile table (extFile column) and 
the offset for the start of the sequence in that file.

ucscRetroExtFile - lists the file ids (extFile field in the Seq table), file 
names, file paths and file sizes for input sequences files.

ucscRetroInfo - retrogene BED format of the PSL lastz alignments and other 
information associated with each retrogene. See below for full description of 
the fields in this table. 

ucscRetroOrtho - for each retrogene id (name column) the nets are listed and
the overlap for that net (this is the break in orthology score displayed in
the "Break in Orthology" table on the Retrogenes track details page for each
retrogene item in the track. For more details on how this score is calculated,
see the netOverlap() function in the
retroFinder/branches/version2/src/pslPseudo/pslPseudo.c source code.

Tables are loaded by the ucscRetroStep5.sh script with the exception of the 
Seq and ExtFile tables which are loaded by the ucscRetroStep2.sh script.
Additionally the input mRNA sequences (e.g. GenBank mRNAs and RefSeq mRNAs), 
aligned by lastz, are symlinked to the /gbdb/<db>/blastzRetro$VERSION
directory. These gbdb directory files are referenced in the ucscRetroExtFile
table. Sequences are in FASTA format and the files reside in the 
/hive/data/genomes/<db>/bed/mrnaBlastz.$VERSION directory.

UCSCRETROINFO TABLE:
The ucscRetroInfo$VERSION table (extended BED12 format) contains both the 
BED format of the retrogene locations on the genome (from the lastz alignments 
of the parent genes to the genome) from the ucscRetroAli$VERSION table and 
additional features of the retrogenes as defined in the UCSC Genome Browser 
kent source tree: kent/src/hg/lib/ucscRetroInfo.as
