import os, sys, re, subprocess
import time
from commonOps import *
from Files import SeqFiles,GeneralFiles

class SeqAndAlignData(object):
    def __init__(self, seqType, alignTable, cfgParse, chromFile=None, isGenePred=False, ensembl=False, ensDb=None):
        self.database = database
        self.seqType = seqType
        self.alignTable = alignTable
        self.cfg = cfgParse
        self.ensembl = ensembl
        self.ensDb = ensDb
        # Get the PSL alignment data for this dataset and for genePred format
        # data, also get the genePred file
        if not isGenePred:
            self.__getPslAlignments(self.alignTable)
            self.getCdsRegions()
        else:
            self.gpFile = self.cfg.getGenePredFile(self.seqType)
            self.__getPslFromGenePreds(self.alignTable)
            # Get the CDS regsions and write to a tab-separated file
            self.getGenePredCdsRegions()
    
    def __getDatabasePrefix(self):
        """Gets the letters of the prefix of the database"""
        prefix = ""
        p = re.match(r'^([A-Za-z]+)\d+', self.database)
        if p:
            prefix = p.group(1)
        print "Prefix is: ", prefix
        return prefix

    def getGenbankSeqs(self, source, seqType):
        """Gets GenBank or RefSeq mRNA sequences"""
        # Source is genbank or refseq
        print "Source ", source
        print "out file for Genbank seqs", outFile
        gbdb = "-db=" + self.__getDatabasePrefix()
        gbR = "-gbRoot=" + self.cfg.getGenVar('gbRoot')

        # Create the output directory for this program 
        makeDir(self.cfg.getSeqDir())
        outFile = self.cfg.getSeqFile(seqType)
        # Program to get GenBank mRNa and RefSeq sequences
        gbSeqProg = self.cfg.getProgVar('genbankSeqProg')
        subprocess.check_call([gbSeqProg, "-inclVersion", "-native", gbdb, gbR, source, "mrna", outFile])

    def getEnsemblSeqs(self, outFile):
        """Gets the Ensembl sequences"""
        # Create the output directory for this program 
        makeDir(self.cfg.getSeqDir())  
        org = getOrganismName(self.database)
        dir = self.cfgParse.getGenVar('scriptDir')
        seqProg = self.cfg.getProgVar('ensemblSeqProg')
        getEnsSeqProg = self.cfgParse.createPath(dir, seqProg)
        # Program to get Ensembl sequences with ids with version numbers
        subprocess.check_call([getEnsSeqProg, "--all", org, "all", outFile])
          
    def getSeqsFromGenePred(self, table, outFile):
        """Gets the sequences from a table or from a file in genePred format"""
        # Create the output directory for this program 
        makeDir(self.cfg.getSeqDir()) 
        rnaProg = self.cfg.getProgVar('getRnaProg') 
        # table can be a table or a file
        subprocess.check_call([rnaProg, self.database, table, "all", outFile])

    def __getPslAlignments(self, alignTable):
        """Gets PSL alignments from the database""" 
        # Create the output directory for this data 
        makeDir(self.cfg.getSeqDir())
        sqlCmd = self.cfg.getProgVar('hgsqlCmd')  
        pslSelect = "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from " + alignTable + ";"
        getPsl = sqlCmd + [pslSelect, self.database]
        with open(self.cfg.getPslFile(), "w") as fh:
            subprocess.check_call(getPsl, stdout=fh)

    def __getGenePredAnnots(self, gpTable):
        """Gets genePred annotations from the database"""
        # Create the output directory for this data 
        makeDir(self.cfg.getSeqDir())
        sqlCmd = self.cfg.getProgVar('hgsqlCmd')  
        genePredSelect = "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from " + gpTable + ";"
        getGp = sqlCmd + [genePredSelect, self.database]
        with open(gpFile, "w") as fh:
            subprocess.check_call(getGp, stdout=fh)
        if self.ensembl:
             # This program gets the version numbers for ids and re-writes 
             # the genePred with the ids with version numbers.
             getEns = self.cfg.createPath(self.cfg.getGenVar('scriptDir'), \
                 self.cfg.getProgVar('ensemblSeqProg'))
             subprocess.check_call([getEns, self.database, self.ensDb, \
                 self.gpFile, self.cfg.getSeqDir()])
    
    def __convertGenePredToPsl(self):
        """Converts annotation data from genePred to PSL format."""
        makeDir(self.cfg.getSeqDir()) 
        getPsl = self.cfg.getProgVar('gpToPsl')
        subprocess.check_call([getPsl, self.cfg.chromFile, self.gpFile, \
            self.cfg.getPslFile(self.seqType)])
   
    def __getPslFromGenePreds(self, gpTable):
        """Given a genePred table, creates a file of the data and 
           converts it to PSL format"""
        gpFile = self.cfg.getGenePredFile(self.seqType)
        self.__getGenePredAnnots(gpTable)
        self.__convertGenePredToPsl()
 
    def getGenePredCdsRegions(self):
        """Gets the CDS region for a genePred annotation table/file"""
        makeDir(self.cfg.getSeqDir()) 
        subprocess.check_call(["./getGenePredCdsRegions", self.database, \
            self.gpFile, self.cfg.getCdsFile(self.seqType)])  
        
    def getCdsRegions(self):
        """Gets the CDS regions either for GenBank mRNAs and RefSeqs"""
        makeDir(self.cfg.getSeqDir()) 
        cdsSelStr = "select acc, version, name, type from " + alignTable + \
            " as a, gbCdnaInfo as g, cds as c where qName = acc and cds = c.id"
        hgsql = self.cfg.getProgVar('hgsqlCmd')
        query = hgsql + [cdsSelStr, self.database]
        with open(cdsFile, "w") as fh:
            subprocess.check_call(query, stdout=fh)
        fh.close()
