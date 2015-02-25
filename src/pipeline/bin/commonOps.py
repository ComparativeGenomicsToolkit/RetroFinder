# Set of common functions used throughout the RetroFinder pipeline
import os,subprocess,re
import time

def queryDb(query, db, fh=None):
    """Make a query to a local database and optionally save output to a file"""
    hgsqlCmd = ['hgsql', '-Ne', query, db]
    print "Cmd is", hgsqlCmd
    if fh != None:
        print "Querying database"
        subprocess.check_call(hgsqlCmd, stdout=fh)
    else:
        subprocess.check_call(hgsqlCmd)

def getDbQueryResult(query, db):
    """Make a query to a local database and return the output"""
    hgsqlCmd = ['hgsql', '-Ne', query, db]
    # Returns result of query and removes newline at the end
    return (subprocess.check_output(hgsqlCmd)[0:-1])

def fileExists(fileName):
    """Checks that the file exists, returns True if it does"""
    return os.path.exists(fileName)
 
def createPath(dir, fileOrDir):
    """Adds file or directory to path name to create a full file
       or directory path"""
    return dir + "/" + fileOrDir

def createFilePath(dir, name, suffix):
    """Creates the file name and path for output files."""
    return dir + "/" + name + "." + suffix

def makeFileList(path, fileExt):
    """Returns a list of files ending in file extension"""
    # return [fn for fn in os.listdir(path) if any([fn.endswith(fileExt)])]
    list = []
    for fn in os.listdir(path): 
        if fn.endswith(fileExt): 
            fn = path + "/" + fn
            list.append(fn)
    return list

def catFiles(outFile, fileList):
    """Cat files in list into outFile"""
    fileStr = ""
    catCmd = ["cat"] + fileList
    with open(outFile, "w") as outFh:
        subprocess.check_call(catCmd, stdout=outFh)
    outFh.close()

def removeFile(file):
    subprocess.check_call(["rm", "-f", file])

def removeDir(dir):
    subprocess.check_call(["rm", "-rf", dir])

def moveFile(oldFile, newFile):
    """Move oldfile to newFile. newFile can be file name or directory"""
    subprocess.check_call(["mv", oldFile, newFile])

def getDate():
    return time.strftime("%Y-%m-%d") 

def getOrganismName(assembly):
    """Give the name of the database return the common name."""
    selectStr = "select organism from dbDb where name =" + "'" + assembly + "';" 
    org = getDbQueryResult(selectStr, "hgcentraltest")
    return org.lower()

def makeDir(dir):
    subprocess.check_call(["mkdir", "-p", dir])
    
def getChromSizes(database, chromFile):
    """Get chromosome sizes"""
    chromStr = "select chrom, size from chromInfo;" 
    with open(chromFile, "w") as fh:
        subprocess.check_call(["hgsql", "-Ne", chromStr, database], stdout=fh)

def removeIdVersion(acc):
    """Removes the .1,.2 etc. suffix from an id"""
    id = ""
    i = re.match(r'^([A-Z0-9]+)\.\d+', acc)
    if i:
        id = i.group(1)
    return id

def getIdVersion(acc):
    """Returns the version number of an id"""
    version = ""
    i = re.match(r'^[A-Z0-9]+\.(\d+)', acc)
    if i: 
        version = i.group(1)
    return version

class TabFileTbl(object):
    """Read from a tabbed file and store in a dictionary keyed by specified
       column."""
    def __init__(self, tabFile, keyCol):
         self.tabFile = tabFile
         # Dictonary of rows from tabbed file
         self.data = None
         # index (0-based) of column to use as hash key
         self.keyCol = keyCol
         self.__buildDictionary()

    def __buildDictionary(self):
        """Create the dictionary with unique keys"""
        self.data = dict()
        for row in TabFileReader(self.tabFile):
            if row[self.keyCol] in self.data:
                raise Exception("This value is already in index: " + row[self.keyCol])
            self.data[row[self.keyCol]] = row
    
class TabFileReader(object):
    """Read from a tabbed file"""
    def __init__(self, tabFile):
        self.fh = open(tabFile, "r");

    def __iter__(self):
        return self
 
    def next(self):
        while True:
            line = self.fh.readline()
            if (line == ""):
                self.fh.close
                raise StopIteration
            if not ((len(line) == 1) or line.startswith('#')):
                line = line[0:-1]  # drop newline
                return line.split("\t")
