#!/usr/bin/python

import sys
import numpy as np
import density
from density import DensityHeader, readDensityData, writeDensityData
from miscellaneous import throwError


programName = sys.argv[0].rsplit('/')[-1]
help = """
This program is used to write a binary density file from a raw binary file.

Usage: %s  inputFile  outputFile  gridSize  boxLength  omega0
\n""" % programName



# read the arguments supplied by the user
if len(sys.argv)<6:
    print help
    sys.exit(1)

inputFile  = sys.argv[1]
outputFile = sys.argv[2]
gridSize   = int( sys.argv[3] )
boxLength  = float( sys.argv[4] )
omega0     = float( sys.argv[5] )
programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])


# read the binary raw file
dataSize = gridSize * gridSize * gridSize
data = np.fromfile( inputFile, np.float32, dataSize )
#data = np.abs( data )



#fill in a density header file
header = DensityHeader()
header.gridSize[:] = gridSize
header.totalGrid = np.uint64(dataSize)
header.fileType  = np.int32( density.DENSITY_FILE )
header.box[:] = 0., boxLength, 0., boxLength, 0., boxLength
header.Omega0  = np.float64( omega0 )
header.OmegaLambda  = np.float64( 1.-omega0 )
header.HubbleParam  = np.float64( .7 )
header.observations.put( range(len(programOptionsDesc)), programOptionsDesc)

# write the data
writeDensityData( outputFile, header, data )



