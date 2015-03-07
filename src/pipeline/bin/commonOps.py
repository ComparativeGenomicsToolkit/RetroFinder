# Set of common functions used throughout the RetroFinder pipeline
import os,sys,subprocess,re
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
    exists = os.path.exists(fileName)
    if not exists:
        raise Exception("The file, %s, does not exist.\n" % (fileName))
    return exists
 
def openReadFileHandle(fileName):
    """Tries opening a file, if it does not exist, exits and if it does
       then returns the filehandle for reading."""
    fh = None
    try:
        fh = open(fileName, "r")
    except IOError:
        print "The file, %s, does not exist and can not be opened, \
            exiting...\n" \
            % (fileName)
        sys.exit(0)
    return fh

def getListFromFileColumn(file, col):
    """Returns a list of the items in column, col, in the file."""
    return subprocess.check_output(["cut", "-f" + str(col), file])
    
def createPath(dir, fileOrDir):
    """Adds file or directory to path name to create a full file
       or directory path"""
    return dir + "/" + fileOrDir

def createFilePath(dir, name, suffix):
    """Creates the file name and path for output files."""
    return dir + "/" + name + "." + suffix

def makeFileList(path, fileExt):
    """Returns a list of full paths for files ending in file extension"""
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

def removeFiles(fileList):
    # removes file(s) in the list
    remove = ["rm", "-f"]
    rmCmd = remove + fileList
    subprocess.check_call(rmCmd)

def removeDir(dir):
    subprocess.check_call(["rm", "-rf", dir])

def moveFile(oldFile, newFile):
    """Move oldfile to newFile. newFile can be file name or directory"""
    subprocess.check_call(["mv", oldFile, newFile])

def moveMultFiles(filesPat, dir):
    """Moves multiple files with pattern ending in "*" or just "*" to dir.
       Full paths to files and directory should be provided. """
    mvStr = "mv " + filesPat + " " + dir 
    subprocess.check_call(mvStr, shell=True)

def getSubDirsList(dir):
    """Returns the list of subdirectories contained in dir assuming there 
       is only a single level of subdirectories required"""
    return [name for name in os.listdir(dir)
        if os.path.isdir(os.path.join(dir, name))]

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
