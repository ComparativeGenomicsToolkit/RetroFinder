import os, sys

# Classes for reading and storing information from a tab-separated file

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
