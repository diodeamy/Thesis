#!/usr/bin/python


import os
import sys
import numpy as np
import textFile
from miscellaneous import throwError



help = '''Use this program to combine the results of several files. It copies the specified columns from one file into the output file:
USAGE:  %s  outputFile   file1,noDescLines,noColumns,columnsToCombine   file2,..,..,..  file3,..,..,..  ....
where:
    fileX - is the name of the x-th file
    noDescLines - the number of lines to ignore at the begining of the file
    noColumns -  number of columns in the file
    columnsToCombine - which columns to keep for the output file (e.g. '2:5' '2,3,7,9' ...)
''' % sys.argv[0]


if len(sys.argv)<=2:
    print help
    sys.exit(1)


outputFile = sys.argv[1]
fileList = sys.argv[2:]

outData = None
for f in fileList:
    f1 = f.split(',')
    fileName, noDescLines, noColumns = f1[0], int(f1[1]), int(f1[2])
    f2 = open( fileName, 'r' ); noRows = len(f2.readlines()) - noDescLines; f2.close()
    desc, data = textFile.readTextFile( fileName, noDescLines, noRows, noColumns, np.float64, VERBOSE=True )
    for f3 in f1[3:]:
        if outData is None:
            outData = eval( 'data[:,%s]' % f3 )
        else:
            outData = np.column_stack( (outData, eval( 'data[:,%s]' % f3 )) )
    


descr = "# %s\n\n" % ( " ".join( sys.argv ) )
textFile.writeTextFile( outputFile, descr, outData, VERBOSE=True )
