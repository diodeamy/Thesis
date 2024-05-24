#!/usr/bin/python2.6

import sys
import numpy as np
from MMF import MMFHeader, readMMFData, writeMMFData
import MMF
from miscellaneous import throwError
from scipy.weave import inline, converters


programName = sys.argv[0].rsplit('/')[-1]
help = """ Used to combine the MMF environments: node, filaments and walls into a single file.

Usage:  %s  node_MMF_clean_file  fila_MMF_clean_file  wall_MMF_clean_file  output_file

The input files must be MMF clean response files (1 if a feature is present, 0 otherwise) - in short int format. The combined MMF map is saved in the output file (which contains one short int for each grid cell) using the following convention:
    nodes = 4
    filaments = 3
    walls = 2
    everything else = 0
"""  % programName


if len(sys.argv)!=5:
    print help
    sys.exit(1)

nodeClean = sys.argv[1]
filaClean = sys.argv[2]
wallClean = sys.argv[3]
outputFile = sys.argv[4]
environments = {'node':4, 'fila':3, 'wall':2}
programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])



# read the wall clean response
[header, data] = readMMFData(wallClean)
if header.feature!=MMF.MMF_WALL:   # test that this is a wall file
    throwError( "File '%s' should contain the MMF clean response for walls. But it contain an MMF response for the '%s' feature." % (wallClean,header.Feature()) )
if header.fileType!=MMF.MMF_CLEAN_RESPONSE: # this test if file is a clean response
    throwError( "File '%s' should contain the MMF clean response for walls. But it does not contain an MMF clean response, but it contains '%s'." % (wallClean,header.FileType()) )
ids = data==1
data[ids] = environments['wall']


# read the filament clean response
[tempHeader, tempData] = readMMFData(filaClean)
if tempHeader.feature!=MMF.MMF_FILAMENT:   # test that this is a node file
    throwError( "File '%s' should contain the MMF clean response for filaments. But it contain an MMF response for the '%s' feature." % (filaClean,tempHeader.Feature()) )
if tempHeader.fileType!=MMF.MMF_CLEAN_RESPONSE: # this test if file is a clean response
    throwError( "File '%s' should contain the MMF clean response for filaments. But it does not contain an MMF clean response, but it contains '%s'." % (filaClean,tempHeader.FileType()) )
if not (tempHeader.gridSize==header.gridSize).all():
    throwError( "The filament MMF clean response from file '%s' does not have the same grid dimensions as the wall MMF clean response. Filament grid is '%s' while node grid is '%s'." % (filaClean,str(tempHeader.gridSize),str(header.gridSize)) )
ids = tempData==1
data[ids] = environments['fila']


# read the node clean response
[tempHeader, tempData] = readMMFData(nodeClean)
if tempHeader.feature!=MMF.MMF_NODE:   # test that this is a node file
    throwError( "File '%s' should contain the MMF clean response for nodes. But it contain an MMF response for the '%s' feature." % (nodeClean,tempHeader.Feature()) )
if tempHeader.fileType!=MMF.MMF_CLEAN_RESPONSE: # this test if file is a clean response
    throwError( "File '%s' should contain the MMF clean response for nodess. But it does not contain an MMF clean response, but it contains '%s'." % (nodeClean,tempHeader.FileType()) )
if not (tempHeader.gridSize==header.gridSize).all():
    throwError( "The node MMF clean response from file '%s' does not have the same grid dimensions as the wall MMF clean response. Wall grid is '%s' while node grid is '%s'." % (wallClean,str(tempHeader.gridSize),str(header.gridSize)) )
ids = tempData==1
data[ids] = environments['node']



# Set some values in the output header and write the MMF combined data to the output file
header.feature = np.int32( MMF.MMF_ALL )
header.filter = np.int32(50)
header.fileType = np.int32( MMF.MMF_CLEAN_RESPONSE_COMBINED )
header.AddProgramCommands( programOptionsDesc )

writeMMFData(outputFile,header,data)

