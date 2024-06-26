import ConfigParser
from commonOps import *

class ParseConfig:
    def __init__(self, configFile):
        # Read in config file
        self.config = ConfigParser.ConfigParser()
        self.config.read(configFile)
        self.db = self.getGenVar('database')
        self.date = self.getGenVar('date') 
        self.version = self.getGenVar('version')
        # root of the directory for this RetroFinder pipeline run
        self.rootRunDir = self.createRootRunDirName()
        print "root run dir", self.rootRunDir
        # file for chromosome sizes
        self.chromFile = createPath(self.rootRunDir, \
            self.getSeqVar('chromFile'))  
        print "chrom sizes file:", self.chromFile
       
    def createRootRunDirName(self):
        """Creates path for directory for RetroFinder run"""
        rootWorkDir = self.getGenVar('rootWorkDir')
        workDir = self.db + "/version" + self.version + "/" + self.date
        return createPath(rootWorkDir, workDir)
    
    def createWorkingDirName(self, dirPath):
        """Creates path and name of directory in root dir"""
        return createPath(self.rootRunDir, dirPath)

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

    def getAlignVar(self, var):
        """Returns variable value from Alignment section of config file"""
        return self.getVar('Alignment', var)

    def getRetroPredVar(self, var):
        """Returns variable value from RetroPred section of config file"""    
        return self.getVar('RetroPred', var)

    def getTempDir(self):
        """Returns full path of temp directory"""
        return os.environ['TMPDIR'] + "/" + self.getGenVar('tempDirSuffix')

    def getChromFile(self, dirPath):
        """Returns full path for chromosome sizes file"""
        return createPath(self.rootRunDir, self.getSeqVar('chromFile'))
    
    def getSeqDir(self):
        """Returns the full path of the sequences directory"""
        seqDir = self.createWorkingDirName(self.getSeqVar('seqDir'))
        return seqDir

    def getAlignDir(self):
        """Returns the full path of the mRNAs alignment directory"""
        alignDir = self.createWorkingDirName(self.getAlignVar('alignDir'))
        return alignDir

    def getRetroPredAnnotDir(self):
        """Returns the full path of the directory of annotations for
           input to pslPseudo"""
        annotDir = \
            self.createWorkingDirName(self.getRetroPredVar('retroPredAnnot'))
        return annotDir

    def getFileName(self, seqType, ext):
        """Returns file name for sequence-related data""" 
        return self.getSeqVar(seqType) + "." + ext

    def getSeqFile(self, seqType):
        """Returns full path and file name of FASTA file for sequences"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return self.getFileName(seqType, "fa")
 
    def getSeqFileFullPath(self, seqType):
        """Returns file name and path to sequence directory"""
        return createPath(self.getSeqDir(), self.getSeqFile(seqType))
   
    def getSeqVarFullPath(self, var):
        """Returns file name and path to sequence directory. Looks up 
           file name in config file using var."""
        return createPath(self.getSeqDir(), self.getSeqVar(var))
 
    def getPslFile(self, seqType):
        """Returns file name of PSL file for sequences"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return self.getFileName(seqType, "psl")
    
    def getPslFileFullPath(self, seqType):
        """Returns file name and path to sequence directory"""
        return createPath(self.getSeqDir(), self.getPslFile(seqType))
 
    def getGenePredFile(self, seqType):
        """Returns file name of genePred file for sequences"""
        # seqType can be ensembl or other genePred annotation
        return self.getFileName(seqType, "gp")

    def getGenePredFileFullPath(self, seqType):
        """Returns file name and path to sequence directory"""
        return createPath(self.getSeqDir(), self.getGenePredFile(seqType))

    def getCdsFile(self, seqType):
        """Returns file name for CDS regions file"""
        # seqType can be mrna, refSeq, ensembl, all or anther type    
        return self.getFileName(seqType,"cds.tab")

    def getCdsFileFullPath(self, seqType):
        """Returns file name and path to sequence directory"""
        return createPath(self.getSeqDir(), self.getCdsFile(seqType))

    def getAlignFileFullPath(self, fileVar):
        """Returns full path to file or directory in the alignments directory"""
        return createPath(self.getAlignDir(), self.getAlignVar(fileVar))
    
    def getRetroPredFullPath(self, fileVar, suffix):
        """Returns full path to file or directory in the pslPseudo
           input annotations directory"""
        return createFilePath(self.getRetroPredAnnotDir(), \
            self.getRetroPredVar(fileVar), suffix)

    def getRetroPredRepeatFullPath(self, repeat):
        """Returns full path to Repeats file in the pslPseudo input 
           annotations directory"""
        # repeat is the variable name for the repeats table in the config file
        return self.getRetroPredFullPath(repeat, "bed")
     
    def getRetroPredNetFullPath(self, net):
        """Returns full path to Nets file in the pslPseudo input 
           annotations directory"""
        # net is the variable name for the net table in the config file
        return self.getRetroPredFullPath(net, "txt")
    
    def getRetroPredGeneFullPath(self, geneset):
        """Returns full path to Geneset file in the pslPseudo input 
           annotations directory"""
        # geneset is the variable name for the geneset table in the config file
        return self.getRetroPredFullPath(geneset, "gp")
