
#! /bin/env python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy

####################################################################################
# path psyco
####################################################################################
#sys.path.append("/home/jsp/prog/utillities/py_modules")
# to speed things up
#import psyco
#psyco.full()
#import bed



######################

def clusterHaplotypes(inFile, outFile):

    # create a list with zeros to store the number of Ns per strain
    for line in inFile:
        ignoreLine=False
        nucs={'A':0,'T':0,'G':0,'C':0,'.':0}
        lineList=line.strip('\n').split(',')
        position=lineList[0]
        for i in range(1, 30):
                if lineList[i] not in nucs.keys():
                    ignoreLine=True
                else: nucs[lineList[i]] +=1
        numNucs=0
        # count number of different alleles to see if the site is polymorphic:
        for n in ['A','T','G','C']:
            if nucs[n] >0:
                numNucs +=1
        if (numNucs >1 or nucs['.']>0) and ignoreLine==False:
            outFile.write(line)



    
###############


def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input.bed> <output.bed> <threshold>
    %prog filters out the lines that don't meet a certain threshold. """

    parser = OptionParser(usage)
   

    return parser



def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    inFN         = args[0]
    outFN        = args[1]

    if inFN == '-':
        inFile = sys.stdin
    else:
        inFile      = open(inFN, 'r')

    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')



    clusterHaplotypes(inFile, outFile)


    

#run main
if __name__ == '__main__':
    main()
