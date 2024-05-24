#!/usr/bin/python

import numpy as np
import sys
import os
import MMF as NEXUS
import miscellaneous as misc
import textFile


# set the details of the input file
noColumns   = 4
posColumns  = [ 0, 1, 2 ]
mpcUnit     = 1000.
noDescriptionLines = 0




# the actual program    

programName = sys.argv[0].rsplit('/')[-1]
help = '''Use this program to read in halo data and split the results according to their NEXUS environment.
USAGE:   %s   input_halo_file   nexus_environment_file   output_file   
''' % programName


if len(sys.argv)<4:
    print help
    sys.exit(1)

programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])
haloFile   = sys.argv[1]
nexusFile  = sys.argv[2]
outputFile = sys.argv[3]



# read the halo data from the input file
if not os.path.isfile( haloFile ):
    misc.throwError( "Cannot find the input halo file '%s'!" % haloFile )

# get the number of rows in the file
lines = 0
for line in open(haloFile):
    lines += 1
N = lines - noDescriptionLines
print "Found %i haloes in the file '%s'. Reading the input data ..." % (N,haloFile)

# read the actual data
halo = np.fromfile( haloFile, np.float64, N*noColumns, sep=' ' )
halo.shape = N, noColumns



# read the nexus data from the input file
header, data = NEXUS.readMMFData( nexusFile )
data[data!=0] -= 1
data.shape = header.gridSize[0], header.gridSize[1], header.gridSize[2]



# get the halo grid indices
offset = header.box[0::2] * mpcUnit    #offset of box coordinates
dx = (header.box[1::2] - header.box[0::2]) / header.gridSize[:] * mpcUnit   #dx for each grid axis
gridPos = [ 0,1,2 ]
for i in range(3):
    gridPos[0] = ( ( halo[:,posColumns[i]] - offset[i] ) / dx[i] ).astype(np.int32)
haloEnv = data[ gridPos[0], gridPos[1], gridPos[2] ]


# loop over the environments and output the haloes in that environment
envs = { 0:'void', 1:'wall', 2:'fila', 3:'node', 4:'all' }
for i in range(5):
    select = haloEnv==i
    outputHaloes = halo[select,:]
    if i==4: outputHaloes = halo
    
    # write the data
    desc = ''
    textFile.writeTextFile( '%s_%s.txt'%(outputFile,envs[i]), desc, outputHaloes )
    
