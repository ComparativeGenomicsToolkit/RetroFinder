from __future__ import division
import string

#  Classes for storing information about a genePred and converting chromosome
#  coordinates between those for genome and transcript
class ConvertCoords(object):
    "object that holds information about converting coordinates"
    def __init__(self, genePred, pos, type):
        self.seqId = genePred.name
        self.gp = genePred
        # Type of coordinate e.g. transcript or genomic 
	self.origType = type
        # Orginal coordinate, 0-based
        self.origPos = int(pos)
        self.newType = self.__setNewType()
        # New coordinate, 0-based
        self.newPos = self.__convertCoords()

    def __setNewType(self):
        newType = None
        if self.origType == "transcript":
            newType = "genomic"
        elif self.origType == "genomic":
            newType = "transcript"
        else:
            print "Error: Unrecognised type."
        return newType

    def __convertTxToGenomicPos(self):
        oldPos = 0
        newCoord = -1
        exons = self.gp.exons
        if self.gp.strand == "-":
            # Reverse the array or Exon objects only if not done so already
            if exons[0].iExon == 0:
                exons.reverse()
        #print exons[0].start
        for i in range(0, len(exons)):
            size = exons[i].size()
            newPos = oldPos + size
            if newPos > self.origPos:
                # oldPos is in 1-based coordinates, origPos is 0-based
                diff = self.origPos - oldPos 
                # newCoords should be 0-based
                if self.gp.strand == "+":
                    newCoord = exons[i].start + diff
                else:
                    newCoord = exons[i].end - diff - 1
                return newCoord
            oldPos = newPos
        # Returns -1 if the new positions is not calculated as above
        # i.e. something went wrong.
        return newCoord
    
    def __convertGenomicToTxPos(self):
        """Converts genomic coordinates to transcript coordinates based
           on a genePred format file."""
        oldPos = 0
        diff = 0
        newCoord = -1
        exons = self.gp.exons
        if self.gp.strand == "-":
            # Reverse the array or Exon objects only if not done so already
            if exons[0].iExon == 0:
                exons.reverse()
        for i in range(0, len(exons)):
            size = exons[i].size()
            newPos = oldPos + size
            if self.gp.strand == "+":
                if exons[i].end >= self.origPos:
                    # oldPos is in 1-based coordinates
                    # origPos is 0-based,  
                    diff = self.origPos - exons[i].start
                    newCoord = oldPos + diff
                    return newCoord
            elif self.gp.strand == "-": 
                if exons[i].start <= self.origPos:
                    # oldPos is in 1-based coordinates, origPos is 0-based
                    # newCoords should be 0-based
                    diff = exons[i].end - self.origPos - 1
                    newCoord = oldPos + diff
                    return newCoord
            else:
                print "Error: genePred has incorrect strand type"    
            oldPos = newPos 
        # Returns -1 if the new positions is not calculated as above
        # i.e. something went wrong.
        return newCoord
        
    def __convertCoords(self):
        if self.newType == "genomic":
            return self.__convertTxToGenomicPos()
        elif self.newType == "transcript":
            return self.__convertGenomicToTxPos()

    def __str__(self):
        return "Sequence ID: " + self.seqId + " Orig Type: " + self.origType + " Orig Pos: " + str(self.origPos) + " New Type: " + str(self.newType) + " New Pos: " + str(self.newPos) + "\n"   
