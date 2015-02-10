import ConfigParser
from commonOps import *

configFile="configFile.cfg"
class ParseConfig:
    def __init__(self, database):
        self.config = ConfigParser.ConfigParser()
        self.config.read(configFile)
        self.db = database
        self.date = getDate()
        self.rootRunDir = self.createRootRunDirName()
        print "Root run directory is", self.rootRunDir

    def writeConfig(self, cfgFile):
        """Writing our configuration file to configFile.cfg"""
        with open(cfgFile, 'a') as cfh:
            self.config.write(cfh)

    def createRootRunDirName(self):
        rootWorkDir = self.config.get('General', 'rootWorkDir', 0) 
        return  rootWorkDir + "/" + self.db + "/" + self.date
   
    def createWorkingDirName(self, rootDir, dirPath):
        return rootDir + "/" + dirPath