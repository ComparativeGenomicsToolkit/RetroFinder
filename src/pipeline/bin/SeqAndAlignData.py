import os, sys, re, subprocess
import time
from commonOps import *
from Files import SeqFiles

# Fetches Genbank and RefSeq sequences
genbankSeqProg = "/cluster/data/genbank/bin/x86_64/gbGetSeqs"
gbRoot = "/cluster/data/genbank"
# Fetches Ensembl sequences 
ensemblSeqProg = "/hive/users/hartera/GencodeWG/retroFinder/trunk/src/pipeline/getEnsemblTxSeqsWithVersions.pl"
# Converts genePred to PSL format
gpToPsl = "genePredToPsl"
# TEMP VARIABLE FOR TESTING
seqType = "mrna"

class SeqAndAlignData(object):
    def __init__(self, database, alignTable, seqFiles):
        self.database = database
        self.seqFiles = seqFiles
        if seqFiles.genePred:
             # Get the genePred format annotations
        # Get the PSL alignment data for this dataset
        self.__getPslAlignments(alignTable)
        
        print "seqFile: %s alignFile: %s genePredFile: %s \n" % (self.seqFiles.seqFile, self.seqFiles.alignFile, self.seqFiles.genePredFile)
        # have an if statement to get Ensembl sequences, do this from 
        # getEnsemblData
       #  self.getEnsemblSeqs(self.ensSeqFile)

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
          
    def __getPslAlignments(self, alignTable):
        """Gets PSL alignments from the database""" 
        makeDir(self.seqFiles.seqsDir)  
        pslSelect = "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from " + alignTable + ";"
        with open(self.seqFiles.alignFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", pslSelect, self.database], stdout=fh)

    def __getGenePredAnnots(self, gpTable):
        """Gets genePred annotations from the database"""
        genePredSelect = "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds"
        with open(self.seqFiles.genePredFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", genePredSelect, self.database], stdout=fh)
     
    def convertGenePredToPsl(self, gpFile, pslFile):
        """Converts annotation data from genePred to PSL format."""
        # NOTE: need to know location of chromFile
        makeDir(self.seqFiles.seqsDir)  
        subprocess.check_call(gpToPsl, self.chromFile, gpFile, pslFile)

     
