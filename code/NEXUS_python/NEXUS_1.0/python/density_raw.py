#!/usr/bin/python

import sys
import numpy as np
import numpy.fft as fft
from density import DensityHeader, readDensityData, writeDensityData
from MMF import MMFHeader, readMMFData, writeMMFData, getHeaderType
import MMF
import mathTools
from miscellaneous import throwError



help = """
This program is used to output a density or MMF binary file type into a raw binary file.

Usage: %s  inputFile  outputFile  (file_type)

The additional options 'file_type' are ['--raw' [DEFAULT],'--grid']:
     '--raw' = outputs only the raw data values
     '--grid' = outputs the grid size along each dimension in 'int' format followed by the raw data
     '--box' = outputs: grid size in 'int', box coordinates (xMin,xMax,yMin,yMax,zMin,zMax) in float and the raw data
     '--fortran' = outputs any of the above using the fortran way of reading/writing binary data.\n""" % sys.argv[0].rsplit('/')[-1]
allowedFiles = ['Density', 'MMF']
fileTypeOptions = { '--raw':False, '--grid':False, '--box':False, '--fortran':False }


# read the arguments supplied by the user
if len(sys.argv)<2 or len(sys.argv)>5:
    print help
    sys.exit(1)
elif len(sys.argv)>=4:
    for i in range(3,len(sys.argv)):
        if sys.argv[i] in fileTypeOptions.keys(): fileTypeOptions[ sys.argv[i] ] = True
        else:
            print "Unrecognized value '%s' for the file_type program option." % sys.argv[i]
            print help
            sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
countTrue = fileTypeOptions['--raw']+fileTypeOptions['--grid']+fileTypeOptions['--box']
if 1>countTrue:
    fileTypeOptions['--raw'] = True
elif 2<=countTrue: throwError( "You inserted more than one option for the output file type!" )
programOptions = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])


# read the density file
headerType = getHeaderType(inputFile)
if headerType not in allowedFiles:
    throwError( "Unrecognized input file type in the program '%s'. Allowed file types are '%s', while the input file is of type: '%s'." % (sys.argv[0],str(allowedFiles),headerType) )

header, data = None, None
if headerType is allowedFiles[0]:
    header, data = readDensityData(inputFile)
elif headerType is allowedFiles[1]:
    header, data = readMMFData(inputFile)

#if header.fileType==MMF.MMF_CLEAN_RESPONSE_COMBINED:
#    data = data.astype(np.int32)


# output the data
dataSize = np.array( [data.nbytes], np.int64 )
headerSize = np.array( [0], np.int64 )
if fileTypeOptions['--raw']:
    print "Writing to a raw binary file the grid data ..."
    f = open( outputFile, 'wb' )
    
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    data.tofile( f )
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    f.close()

elif fileTypeOptions['--grid']:
    print "Writing to a raw binary file the grid size and the grid data ..."
    f = open( outputFile, 'wb' )
    
    grid = header.gridSize.astype(np.int32)
    headerSize[0] = grid.nbytes
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    grid.tofile( f )
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    data.tofile( f )
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    f.close()

elif fileTypeOptions['--box']:
    print "Writing to a raw binary file: the grid size, the box coordinates and the grid data ..."
    f = open( outputFile, 'wb' )
    
    grid = header.gridSize.astype(np.int32)
    headerSize[0] = grid.nbytes
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    grid.tofile( f )
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    
    grid = header.box.astype(np.float32)
    headerSize[0] = grid.nbytes
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    grid.tofile( f )
    if fileTypeOptions['--fortran']: headerSize.tofile( f )
    
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    data.tofile( f )
    if fileTypeOptions['--fortran']: dataSize.tofile( f )
    f.close()
