#!/usr/bin/python

import numpy as np
import sys
import math
from gadget import GadgetHeader, readGadgetHeader
from density import DensityHeader, readDensityHeader
from MMF import MMFHeader, readMMFHeader, getHeaderType


fileTypes = { 'Density':readDensityHeader, 'MMF':readMMFHeader, 'Gadget':readGadgetHeader }
BOCDataType = { 'c':'BYTE', 'h':'SHORT', 'i2':'SHORT', 'i':'INT', 'i4':'INT', 'f':'FLOAT', 'f4':'FLOAT', 'd':'DOUBLE', 'f8':'DOUBLE'}



def writeBOCFile(bovFile,time,dataFile,grid,dataType,dataName,brickOrigin=np.zeros(3,'f4'),brickLength=np.ones(3,'f4'),byteOffset=0,dataComp=1):
    f = open( bovFile, 'w' )
    print "Writing the '.bov' file '%s' for file '%s' with grid %s ... " % (bovFile,dataFile,str(grid))
    print >>f, "TIME: %f" % time
    print >>f, "DATA_FILE: %s" % dataFile
    print >>f, "# Data size Nx,Ny,Nz"
    print >>f, "DATA_SIZE: %i %i %i" % (grid[2],grid[1],grid[0])
    print >>f, "# Allowable values for DATA_FORMAT are: BYTE,SHORT,INT,FLOAT,DOUBLE"
    print >>f, "DATA_FORMAT: %s" % BOCDataType[dataType]
    print >>f, "VARIABLE: %s" % dataName
    print >>f, "# Endian representation: Intel is LITTLE, many other processors are BIG."
    print >>f, "DATA_ENDIAN: LITTLE"
    print >>f, "CENTERING: zonal"
    print >>f, "BRICK_ORIGIN: %f %f %f" % (brickOrigin[2],brickOrigin[1],brickOrigin[0])
    print >>f, "BRICK_SIZE: %f %f %f" % (brickLength[2],brickLength[1],brickLength[0])
    print >>f, "BYTE_OFFSET: %i" % byteOffset
    print >>f, "DIVIDE_BRICK: false"
    print >>f, "# DATA_BRICKLETS: 2 2 2"
    print >>f, "DATA_COMPONENTS: %i" % dataComp
    f.close()

def sqrtN(value,power):
    temp = value-1
    temp1 = temp**(1./power)
    return int(math.floor(temp1)) + 1



programName = sys.argv[0].split('/')[-1]
help = '''Use this program to write a '.bov' file that will be sued by VisIt to upload the data for visualization.

Usage:  %s  file_1   (file_2 ... file_N)

Where you can give the name of one or several files. The output is a text file with the same name as the input file and extension '.bov'.
''' % programName

if len(sys.argv)<2:
    print help
    sys.exit(1)


# loop over the input arguments
for inputFile in sys.argv[1:]:
    outputFile = inputFile + '.bov'

    # Read the file header
    fileType = getHeaderType(inputFile)
    header = None
    if fileType in fileTypes.keys():
        header = fileTypes[fileType](inputFile,VERBOSE=False)
        print "Input file is a '%s' file." % fileType
    else:
        print "Unknown file type: '%s'. Can only read '%s' file types." % (fileType, str(fileTypes.keys()) )
        sys.exit(1)


    # Now call the function that will write the boc file
    if fileType=='Gadget':
        variableName = 'gadgetParticles'
        temp = header.npartTotal[1]
        grid = [temp,1,1]
        noComponents = 3
        brickOrigin = [0.,0.,0.]
        brickLength = [1.,1.,1.]
        print "The file '%s' contains the following:" % inputFile
        print "\t variable name        = %s" % variableName
        print "\t number of components = %i" % noComponents
        print "\t data type            =", BOCDataType[ 'f4' ]
        print "\t brick origin         =", brickOrigin
        print "\t brick size           =", brickLength
        writeBOCFile( bovFile=outputFile, time=header.time, dataFile=inputFile, grid=grid, dataType='f4', dataName=variableName, brickOrigin=brickOrigin, brickLength=brickLength, byteOffset=280, dataComp=noComponents )
        
    else:
        variableName = header.DataName()
        noComponents = header.DataComponents()
        brickOrigin = header.box[::2]
        brickLength = header.box[1::2] - header.box[::2]
        print "The file '%s' contains the following:" % inputFile
        print "\t variable name        = %s" % variableName
        print "\t number of components = %i" % noComponents
        print "\t data type            =", BOCDataType[ header.DataType() ]
        print "\t brick origin         =", brickOrigin
        print "\t brick size           =", brickLength
        writeBOCFile( bovFile=outputFile, time=header.time, dataFile=inputFile, grid=header.gridSize, dataType=header.DataType(), dataName=variableName, brickOrigin=brickOrigin, brickLength=brickLength, byteOffset=1048, dataComp=noComponents )
