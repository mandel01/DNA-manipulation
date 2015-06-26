#Genome25List.py
#D.R.Mandel/M.J.Mandel
#11/27/2004
#To convert flat file of genome listing into database format of overlapping 25-mers
#in both directions, circularized
#Mandel, M. J. (2005). Nutrient Starvation and Post-Transcriptional Regulation of RpoS in Escherichia coli. PhD Dissertation, Princeton University 1–292.

#10/28/2006--Reminder-input files must be in all lowercase, raw sequence format.

import sys, string, os

def complement(base):
    if base == 'a':
        return 't'
    elif base == 't':
        return 'a'
    elif base == 'g':
        return 'c'
    elif base == 'c':
        return 'g'
    else:
        return 'n'
        print 'found ambiguous base'

def createSingleElements(inp):
    forward = []
    backward = []  # in reverse complement

    fi = open(inp,'r')
    fil = fi.readlines()
    fi.close()

    for x in fil:
        for y in x[:-1]:
            forward.append(y)
            backward.append(complement(y))
    backward.reverse()

    numElems = len(forward)     # length of genome

    # for wraparound append to each to account for circular genome
    for x in forward[:24]:
        forward.append(x)       # genome with 24 bases appended at end to allow
    for x in backward[:24]:     # for circularization
        backward.append(x)
    
    return forward,backward,numElems

def parseElementsForward(lst,fo,numElems):
    y = 0L
    lenLst = len(lst)
    while y < len(lst)-24L:
        b = ''
        for x in range(y,y+25L):
            b = b + lst[x]
            
        indx = y + 13L          # indx starts at 13, which is the middle nucleotide
                                # position of the first 25-mer (1 through 25).
                                # Then indx COUNTS UP.
        if indx > numElems:
            indx = indx - numElems

        fo.write('F' + '\t' + str(indx) + '\t' + b + '\n')
        y = y + 1L

def parseElementsBackward(lst,fo,numElems):
    y = 0L
    lenLst = len(lst)
    while y < len(lst)-24L:
        b = ''
        for x in range(y,y+25L):
            b = b + lst[x]
            
        indx = numElems - y - 12L     # indx starts at (genome size minus 12, which is
                                # the middle nucleotide position of the first 25-mer
                                # of the reverse-complement). Then indx COUNTS DOWN.
        if indx < 1L:
            indx = numElems + indx

        fo.write('R' + '\t' + str(indx) + '\t' + b + '\n')
        y = y + 1L

# start here
try:
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
except:
    print 'format: python genome2db.py {inputfile} {outputfile}'
    sys.exit(0)
    
# creates single element lists
forward, backward, numElems = createSingleElements(inputfile)

fo = open(outputfile,'w')
fo.write('Strand\tPosition\tOligo\n')

# writes lists of 25-mers of forward genome
parseElementsForward(forward,fo,numElems)

# writes lists of 25-mers of reverse genome
parseElementsBackward(backward,fo,numElems)          

fo.close()
