import os, sys, re, subprocess
import time
from commonOps import *

# NOTE: Master script needs to get date once so doesn't change if run scripts on multiple days 
def createRootDirName(outDir, database, date):
    return  outDir + "/" + database + "." + date

def createWorkingDirName(rootDir, dirPath):
    return rootDir + "/" + dirPath

class GeneralFiles(object):
    def __init__(self, rootDir):
        # Chromosome sizes file
        self.chromFile = createFilePath(rootDir, "chrom", "sizes")

class SeqFiles(object):
    def __init__(self, rootDir, database, seqType, isGenePred=False):
        # Database name for genome assembly
        self.database = database
        # Sequence type e.g. mrna, refseq, ensembl
        self.seqType = seqType
        self.seqsDir = createWorkingDirName(rootDir, "sequenceData")
        # Names fo files for sequence-related data
        self.seqFile = createFilePath(self.seqsDir, seqType, "fa")
        print "seq file is", self.seqFile
        # This can be set to None, only required for some datasets
        # It is converted to PSL later
        # Flag to indicate whether this is a genePred annotation
        self.isGenePred = isGenePred
        if (self.isGenePred):
            self.genePredFile = createFilePath(self.seqsDir, seqType, "gp")
        else:
            self.genePredFile = None
        # Alignment format is PSL (.psl)
        self.alignFile = createFilePath(self.seqsDir, seqType, "psl")
        print "seqFile: %s alignFile: %s genePredFile: %s \n" % (self.seqFile, self.alignFile, self.genePredFile)
        # CDS regions file
        cdsStr = "cds." + seqType
        self.cdsFile = createFilePath(self.seqsDir, cdsStr, "tab") 
