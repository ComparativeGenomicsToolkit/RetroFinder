import os, sys, re, subprocess
import time
from commonOps import *
from Files import SeqFiles,GeneralFiles

# Fetches Genbank and RefSeq sequences
genbankSeqProg = "/cluster/data/genbank/bin/x86_64/gbGetSeqs"
gbRoot = "/cluster/data/genbank"
# Fetches Ensembl sequences 
ensemblSeqProg = "/hive/users/hartera/GencodeWG/retroFinder/trunk/src/pipeline/getEnsemblTxSeqsWithVersions.pl"
# Converts genePred to PSL format
gpToPsl = "genePredToPsl"
# Gets cDNA sequence from a genePred table or file
getRnaProg = "getRnaPred"
seqType = "mrna"

class SeqAndAlignData(object):
    def __init__(self, database, alignTable, seqFiles, chromFile=None, ensembl=False, ensDb=None):
        self.database = database
        self.seqFiles = seqFiles
        self.ensembl = ensembl
        self.ensDb = ensDb
        # Get the PSL alignment data for this dataset and for genePred format
        # data, also get the genePred file
        if not self.seqFiles.genePred:
            self.__getPslAlignments(alignTable)
        else:
            self.getPslFromGenePreds(alignTable, self.seqFiles.genePredFile, chromFile)

    def __getDatabasePrefix(self):
        """Gets the letters of the prefix of the database"""
        prefix = ""
        p = re.match(r'^([A-Za-z]+)\d+', self.database)
        if p:
            prefix = p.group(1)
        print "Prefix is: ", prefix
        return prefix

    def getGenbankSeqs(self, source, seqType, outFile):
        print "Source ", source
        print "out file for Genbank seqs", outFile
        gbdb = "-db=" + self.__getDatabasePrefix()
        gbR = "-gbRoot=" + gbRoot
        # Create the output directory for this program 
        makeDir(self.seqFiles.seqsDir)
        subprocess.check_call([genbankSeqProg, "-inclVersion", "-native", gbdb, gbR, source, seqType, outFile])

    def getEnsemblSeqs(self, outFile):
        """Gets the Ensembl sequences"""
        makeDir(self.seqFiles.seqsDir)  
        org = getOrganismName(self.database)
        subprocess.check_call([ensemblSeqProg, "--all", org, "all", outFile])
          
    def getSeqsFromGenePred(self, table, outFile):
        """Gets the sequences from a table or from a file in genePred format"""
        makeDir(self.seqFiles.seqsDir)  
        # table can be a table or a file
        subprocess.check_call([getRnaProg, self.database, table, "all", outFile])

    def __getPslAlignments(self, alignTable):
        """Gets PSL alignments from the database""" 
        makeDir(self.seqFiles.seqsDir)  
        pslSelect = "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from " + alignTable + ";"
        with open(self.seqFiles.alignFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", pslSelect, self.database], stdout=fh)

    def __getGenePredAnnots(self, gpTable, gpFile):
        """Gets genePred annotations from the database"""
        makeDir(self.seqFiles.seqsDir) 
        genePredSelect = "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from " + gpTable + ";"
        with open(gpFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", genePredSelect, self.database], stdout=fh)
        if self.ensembl:
             # This program gets the version numbers for ids and re-writes 
             # the genePred with the ids with version numbers.
             subprocess.check_call(["./getEnsDataWithVersions", self.database, self.ensDb, self.seqFiles.genePredFile, self.seqFiles.seqsDir])
    
    def __convertGenePredToPsl(self, gpFile, chromFile):
        """Converts annotation data from genePred to PSL format."""
        # NOTE: need to know location of chromFile
        makeDir(self.seqFiles.seqsDir) 
        subprocess.check_call([gpToPsl, chromFile, gpFile, self.seqFiles.alignFile])
   
    def getPslFromGenePreds(self, gpTable, gpFile, chromFile):
        self.__getGenePredAnnots(gpTable, gpFile)
        self.__convertGenePredToPsl(gpFile, chromFile)
