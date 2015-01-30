# Set of common functions used throughout the RetroFinder pipeline
import subprocess,re
import time

def catFiles(outFile, fileList):
    """Cat files in list into outFile"""
    fileStr = ""
    for f in fileList:
       if fileList.index(f) == len(fileList) -1:
           fileStr = fileStr + f
       else: 
           fileStr = fileStr + f + '","'
    with open(outFile, "w") as outFh:
        subprocess.check_call(["cat", fileStr], stdout=outFh)

def getDate():
    return time.strftime("%Y-%m-%d") 

def getOrganismName(database):
    """Give the name of the database return the common name."""
    selectStr = "select organism from dbDb where name =" + "'" + database + "';" 
    org = subprocess.check_output(["hgsql", "-Ne", selectStr, "hgcentraltest"])
    org = org.strip()
    return org.lower()

def makeDir(dir):
    subprocess.check_call(["mkdir", "-p", dir])
    
def createFilePath(dir, name, suffix):
    """Creates the file name and path for output files."""
    return dir + "/" + name + "." + suffix

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
