import os, sys, re, subprocess
import time
from commonOps import *

class SeqAndAlignData(object):
    def __init__(self, seqType, alignTable, cfgParse, isGenePred=False, ensembl=False, ensDb=None):
        # Config file   
        self.cfg = cfgParse
        self.database = self.cfg.db
        self.seqType = seqType
        self.alignTable = alignTable
        self.cfg = cfgParse
        self.ensembl = ensembl
        self.ensDb = ensDb
        self.tempDir = self.cfg.getTempDir()
        print "temp dir is", self.tempDir
        self.gpFile = createPath(self.tempDir, \
            self.cfg.getGenePredFile(self.seqType))
        print "genePred file is", self.gpFile
        self.cdsFile = createPath(self.tempDir, \
            self.cfg.getCdsFile(self.seqType))
        # Get the PSL alignment data for this dataset 
        self.__getPslAlignmentsAndCds(self.alignTable, isGenePred)
    
    def __getDatabasePrefix(self):
        """Gets the letters of the prefix of the database"""
        prefix = ""
        p = re.match(r'^([A-Za-z]+)\d+', self.database)
        if p:
            prefix = p.group(1)
        print "Prefix is: ", prefix
        return prefix

    def getGenbankSeqs(self, source, outFile):
        """Gets GenBank or RefSeq mRNA sequences"""
        # Source is genbank or refseq
        print "Source ", source
        gbdb = "-db=" + self.__getDatabasePrefix()
        gbR = "-gbRoot=" + self.cfg.getGenVar('gbRoot')

        # Create the output directory for this program 
        makeDir(self.cfg.getTempDir())
        outFile = createPath(self.tempDir, self.cfg.getSeqFile(self.seqType))
        # Program to get GenBank mRNa and RefSeq sequences
        gbSeqProg = self.cfg.getProgVar('genbankSeqProg')
        subprocess.check_call([gbSeqProg, "-inclVersion", "-native", gbdb, gbR, source, "mrna", outFile])

    def getEnsemblSeqs(self, outFile):
        """Gets the Ensembl sequences"""
        # Create the output directory for this program 
        makeDir(self.cfg.getTempDir())  
        org = getOrganismName(self.database)
        dir = self.cfgParse.getGenVar('scriptDir')
        seqProg = self.cfg.getProgVar('ensemblSeqProg')
        getEnsSeqProg = createPath(dir, seqProg)
        # Program to get Ensembl sequences with ids with version numbers
        subprocess.check_call([getEnsSeqProg, "--all", org, "all", outFile])
          
    def getSeqsFromGenePred(self, file, outFile):
        """Gets the sequences from a table or from a file in genePred format"""
        # Create the output directory for this program 
        makeDir(self.cfg.getTempDir()) 
        # table can be a table or a file
        subprocess.check_call(["getRnaPred",self.database,file, "all", outFile])

    def __getPslAlignmentsFromTable(self, alignTable):
        """Gets PSL alignments from the database""" 
        # Create the output directory for this data 
        makeDir(self.cfg.getTempDir())
        pslSelect = "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from " + alignTable + ";"
        outFile = createPath(self.tempDir, self.cfg.getPslFile(self.seqType))
        with open(outFile, "w") as fh:
            queryDb(pslSelect, self.database, fh)
        fh.close()

    def __getGenePredAnnots(self, gpTable):
        """Gets genePred annotations from the database"""
        # Create the output directory for this data 
        makeDir(self.cfg.getTempDir())
        genePredSelect = "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from " + gpTable + ";"
        with open(self.gpFile, "w") as fh:
            queryDb(genePredSelect, self.database, fh)
        if self.ensembl:
             # This program gets the version numbers for ids and re-writes 
             # the genePred with the ids with version numbers.
             getEnsGp = self.cfg.createPath(self.cfg.getGenVar('scriptDir'), \
                 "getEnsDataWithVersions")
             # genePred will initially be written to the temp directory
             subprocess.check_call([getEnsGp, self.database, self.ensDb, \
                 self.gpFile, self.tempDir])
    
    def __convertGenePredToPsl(self):
        """Converts annotation data from genePred to PSL format."""
        makeDir(self.cfg.getTempDir()) 
        subprocess.check_call(["genePredToPsl", self.cfg.chromFile, \
            self.gpFile, self.cfg.getPslFile(self.seqType)])
   
    def __getPslFromGenePreds(self, gpTable):
        """Given a genePred table, creates a file of the data and 
           converts it to PSL format"""
        self.__getGenePredAnnots(gpTable)
        self.__convertGenePredToPsl()

    def __getGenePredCdsRegions(self):
        """Gets the CDS region for a genePred annotation table/file"""
        makeDir(self.cfg.getTempDir())
        scripts = self.cfg.getGenVar('scriptDir')
        # Get CDS regions for sequences with 1-based start coordinates
        getCds = self.cfg.createPath(scripts, "getGenePredCdsRegions")
        subprocess.check_call([getCds,self.database,self.gpFile,self.cdsFile])  
        
    def __getCdsRegions(self):
        """Gets the CDS regions either for GenBank mRNAs and RefSeqs"""
        makeDir(self.cfg.getTempDir()) 
        cdsSelStr = "select acc, version, name, type from " + self.alignTable \
            + " as a, gbCdnaInfo as g, cds as c where qName = acc and \
            cds = c.id"
        with open(self.cdsFile, "w") as fh:
            queryDb(cdsSelStr, self.database, fh)
        fh.close()

    def __getPslAlignmentsAndCds(self, table, isGenePred):
        """Get PSL alignments and CDS regions. Table is either PSL or genePred
           format. GenePreds are converted to PSL format."""
        # Get PSL alignments directly from a table and write to temp dir
        # CDS regions from a database query, write file to temp dir
        if not isGenePred:
            self.__getPslAlignmentsFromTable(table)
            self.__getCdsRegions()
        else:
            # Get genePred file and convert to PSL. Get CDS regions 
            # from genePred file. Write files to temp directory.
            self.__getPslFromGenePreds(table)
            # Get the CDS regsions and write to a tab-separated file
            self.__getGenePredCdsRegions()
