import ConfigParser
from commonOps import *

configFile="configFile.cfg"
class ParseConfig:
    def __init__(self, database, configFile):
        # Read in config file
        self.config = ConfigParser.ConfigParser()
        self.config.read(configFile)
        self.db = database
        self.date = getDate()  
        self.version = getGenVar('version')
        # root of the directory for this RetroFinder pipeline run
        self.rootRunDir = self.createRootRunDirName()
        # file for chromosome sizes
        self.chromFile = self.createPath(self.rootRunDir, \
            self.getGenVar('chromFile'))  
       
    def createRootRunDirName(self):
        """Creates path for directory for RetroFinder run"""
        rootWorkDir = self.getGenVar('rootWorkDir')
        workDir = self.db + "/" + self.version + "/" self.date
        return createPath(rootWorkDir, workDir)
    
    def createWorkingDirName(self, dirPath):
        """Creates path and name of directory in root dir"""
        return self.createPath(self.rootRunDir, dirPath)

    def createPath(self, dir, fileOrDir):
        """Adds file or directory to path name to create a full file
           or directory path"""
        return dir + "/" + fileOrDir

    def getVar(self, section, var):
        """Returns variable value from config file"""
        return self.config.get(section, var, 0) 

    def getGenVar(self, var):
        """Returns variable value from General section of config file"""
        return self.getVar('General', var)

    def getProgVar(self, var):    
        """Returns variable value from Programs section of config file"""
        return self.getVar('Programs', var)

    def getSeqVar(self, var):
        """Returns variable value from SequenceData section of config file"""
        return self.getVar('SequenceData', var)

    def getSeqDir(self):
        """Returns the full path of the sequences directory"""
        seqDir = createWorkingDirName(self.rootRunDir, \
            self.getSeqVar('seqDir'))
        return seqDir

    def getFileName(self, seqType, ext):
        """Returns file name for sequence-related data""" 
        return self.getSeqVar(seqType) + "." + ext

    def getSeqFile(self, seqType):
        """Returns full path and file name of FASTA file for sequences"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return createPath(self.getSeqDir(), self.getFileName(seqType, "fa"))
 
    def getPslFile(self, seqType):
        """Returns full path and file name of PSL file for sequences"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return createPath(self.getSeqDir(), self.getFileName(seqType, "psl"))
    
    def getGenePredFile(self, seqType):
        """Returns full path and file name of genePred file for sequences"""
        # seqType can be ensembl or other genePred annotation
        return createPath(self.getSeqDir(), self.getFileName(seqType, "gp"))

    def getCdsFile(self, seqType):
        """Returns full path and file name for CDS regions file"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return createPath(self.getSeqDir(), self.getFileName(seqType, "cds.tab"))

    def getChromFile(self, dirPath)
        return createFilePath(self.rootRunDir, \
           self.getSeqVar('chromFile'))
