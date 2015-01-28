import os, sys, re, subprocess
import time
from commonOps import *

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
    def __init__(self, database, seqType, alignTable, outDir, gp=False):
        # Database name for genome assembly
        self.database = database
        # Sequence type e.g. mrna, refseq, ensembl
        self.seqType = seqType
        # Creates an output directory for sequence and alignment data
        self.seqsDir = makeDir(outDir, self.database, "sequenceData")
        print "seqs dir is", self.seqsDir
        # File for sequence data
        self.seqFile = createFilePath(self.seqsDir, seqType, "fa")
        print "seq file is", self.seqFile
        # File for chromosome names and sizes
        # CONSIDER PUTTING THIS IN A GENERAL DIRECTORY AND THAT MAIN 
        # PIPELINE DRIVER SHOULD GET IT
        self.chromFile = createFilePath(self.seqsDir, "chrom", "sizes")
        # Get chromosome sizes and write to chromFile
        getChromSizes(self.database, self.chromFile)
        # Alignment format is PSL (.psl)
        self.alignFile = createFilePath(self.seqsDir, seqType, "psl")
        # Get the PSL alignment data for this dataset
        self.__getPslAlignments(alignTable)
        # Ensembl sequences file
        self.ensSeqFile = createFilePath(self.seqsDir, "ensembl", "fa")
        # This can be set to None, only required for some datasets
        self.genePredFile = createFilePath(self.seqsDir, seqType, "gp")
        print "seqFile: %s alignFile: %s genePredFile: %s \n" % (self.seqFile, self.alignFile, self.genePredFile)
        self.getEnsemblSeqs(self.ensSeqFile)

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
        subprocess.check_call([genbankSeqProg, "-inclVersion", "-native", gbdb, gbR, source, seqType, outFile])

    def getEnsemblSeqs(self, outFile):
        """Gets the Ensembl sequences"""
        org = getOrganismName(self.database)
        subprocess.check_call([ensemblSeqProg, "--all", org, "all", outFile])
          
    def __getPslAlignments(self, seqTable):
        """Gets PSL alignments from the database"""   
        pslSelect = "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from " + seqTable + ";"
        with open(self.alignFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", pslSelect, self.database], stdout=fh)

    def __getChromSizes(self):
        """Get chromosome sizes"""
        chromStr = "select chrom, size from chromInfo;" 
        with open(self.chromFile, "w") as fh:
            subprocess.check_call(["hgsql", "-Ne", chromStr, self.database], stdout=fh)
 
    def convertGenePredToPsl(self, gpFile, pslFile):
        """Converts annotation data from genePred to PSL format."""
        subprocess.check_call(gpToPsl, self.chromFile, gpFile, pslFile)

     
