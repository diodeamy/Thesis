#!/usr/bin/python

import sys
import optparse
import numpy as np
import halo
import MMF as NEXUS
from miscellaneous import throwError, throwWarning


# labels to give to the columns that will be added
envLabel       = halo.ENVIRONMENT_LABEL
filaDirLabel   = halo.ENV_FILA_DIRECTION_LABEL
filaThickLabel = halo.ENV_FILA_THICKNESS_LABEL
filaMassLabel  = halo.ENV_FILA_DENSITY_LABEL
wallDirLabel   = halo.ENV_WALL_DIRECTION_LABEL
wallThickLabel = halo.ENV_WALL_THICKNESS_LABEL
wallMassLabel  = halo.ENV_WALL_DENSITY_LABEL



def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputHaloFile  inputEnvironmentFile  outputFile [options] <arg1> ... <argN>', description="Use this program to add the environment information to the halos in the halo file.")
    parser.add_option('-d','--directionFiles', action="store", nargs=2, default=None, type='string', dest="directionFiles", help="Give the 2 files that store the geometrical direction of filaments and walls. It writes this information in the halo file." )
    parser.add_option('-p','--propertiesFiles', action="store", nargs=2, default=None, type='string', dest="propertiesFiles", help="Give the 2 files that store the properties of filaments and wall environments. It writes this information in the halo file." )
    parser.add_option('--newEnv', action="store_true",  default=False, dest="newEnv", help="Rewrite the environment information in the halo file (if already present).")
    parser.add_option('--newDir', action="store_true",  default=False, dest="newDir", help="Rewrite the environment direction information in the halo file (if already present).")
    parser.add_option('--newProp', action="store_true",  default=False, dest="newProp", help="Rewrite the environment properties information in the halo file (if already present).")
    parser.add_option('--envRoot', action="store", default=None, type='string', dest="envRoot", help="Give the directory where the environment files are stored. The program will read all environment files from that directory." )
    
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    if len(args)<3:
        throwError( "The program needs at least 3 arguments giving the input binary halo file, input environment file and the name of the output file. Use '-h' for the help message.\n", EXIT=False)
        parser.print_help()
        sys.exit( 1 )
    options.inputFile = args[0]
    options.environmentFile = args[1]
    options.outputFile = args[2]
    
    if options.envRoot:
        options.environmentFile = options.envRoot + '/' + options.environmentFile
    if options.envRoot and options.directionFiles:
        temp = []
        for i in range( len(options.directionFiles) ):
            temp.append( options.envRoot + '/' + options.directionFiles[i] )
        options.directionFiles = temp
    if options.envRoot and options.propertiesFiles:
        temp = []
        for i in range( len(options.propertiesFiles) ):
            temp.append( options.envRoot + '/' + options.propertiesFiles[i] )
        options.propertiesFiles = temp
    
    return options



# This is the main program part
options = programOptions()
programName = sys.argv[0].rsplit('/')[-1]
options.programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])



# Read the halo data
header, dataIntegers, dataFloats = halo.readHaloData( options.inputFile )
pos = dataFloats[:,header.positionColumns]



# Read the environment data
nexusHeader, nexusData = NEXUS.readMMFData( options.environmentFile )
nexusData.shape = nexusHeader.gridSize

# get the environment at the halo position
gridPos = halo.gridIndices( pos, nexusHeader.box, nexusHeader.gridSize, VERBOSE=False )
envColumn, hasEnv = header.GetColumnIndex( envLabel )[1:]   # returns true if the halo data already has this column
if not hasEnv or options.newEnv:        # write the environment tag into the dataIntegers columns
    print "\nWriting the environment tag in the halo data ..."
    env = nexusData[ gridPos[0],gridPos[1],gridPos[2] ]
    env[ env!=0 ] -= 1      #field=0, wall=1, filament=2 and cluster=3
    if hasEnv:
        dataIntegers[:,envColumn] = env
    else:
        dataIntegers = np.column_stack( (dataIntegers,env) ).copy()
        header.AddColumnIntegers( envLabel )
    



# Read the halo direction data, if given
dirColumn1, hasDir = header.GetColumnIndex( filaDirLabel[0] )[1:]   # returns true if the halo data already has this column
dirColumn2, hasDir = header.GetColumnIndex( wallDirLabel[0] )[1:]
if options.directionFiles and (not hasDir or options.newDir):
    print "\nWriting the filament and wall directions in the halo data ..."
    h1, d1 = NEXUS.readMMFData( options.directionFiles[0] )
    h2, d2 = NEXUS.readMMFData( options.directionFiles[1] )
    filaDirs = NEXUS.NEXUSEnvironmentProperties( nexusData, d1.reshape(-1,h1.gridSize[1]), 3 )
    wallDirs = NEXUS.NEXUSEnvironmentProperties( nexusData, d2.reshape(-1,h2.gridSize[1]), 2 )
    filaDirs.shape = nexusData.shape[0], nexusData.shape[1], nexusData.shape[2], h1.gridSize[1]
    wallDirs.shape = filaDirs.shape
    
    # copy the directions to the halo data array
    if hasDir:
        dataFloats[:,[dirColumn1,dirColumn1+1,dirColumn1+2]] = filaDirs[ gridPos[0],gridPos[1],gridPos[2], : ]
        dataFloats[:,[dirColumn2,dirColumn2+1,dirColumn2+2]] = wallDirs[ gridPos[0],gridPos[1],gridPos[2], : ]
    else:
        dataFloats = np.column_stack( (dataFloats, filaDirs[ gridPos[0],gridPos[1],gridPos[2], : ], wallDirs[ gridPos[0],gridPos[1],gridPos[2], : ] ) ).copy()
        for i in range(3): header.AddColumnFloats( filaDirLabel[i] )
        for i in range(3): header.AddColumnFloats( wallDirLabel[i] )



# Read the halo properties data, if given
dirColumn1, hasDir = header.GetColumnIndex( filaThickLabel )[1:]   # returns true if the halo data already has this column
dirColumn2, hasDir = header.GetColumnIndex( wallThickLabel )[1:]   # returns true if the halo data already has this column
if options.propertiesFiles and (not hasDir or options.newProp):
    print "\nWriting the filament and wall properties in the halo data ..."
    h1, d1 = NEXUS.readMMFData( options.propertiesFiles[0] )
    h2, d2 = NEXUS.readMMFData( options.propertiesFiles[1] )
    filaDirs = NEXUS.NEXUSEnvironmentProperties( nexusData, d1.reshape(-1,h1.gridSize[1]), 3 )
    wallDirs = NEXUS.NEXUSEnvironmentProperties( nexusData, d2.reshape(-1,h2.gridSize[1]), 2 )
    filaDirs.shape = nexusData.shape[0], nexusData.shape[1], nexusData.shape[2], h1.gridSize[1]
    wallDirs.shape = filaDirs.shape
    
    # copy the properties to the halo data array
    if hasDir:
        dataFloats[:,[dirColumn1,dirColumn1+1]] = filaDirs[ gridPos[0],gridPos[1],gridPos[2], : ]
        dataFloats[:,[dirColumn2,dirColumn2+1]] = wallDirs[ gridPos[0],gridPos[1],gridPos[2], : ]
    else:
        dataFloats = np.column_stack( (dataFloats, filaDirs[ gridPos[0],gridPos[1],gridPos[2], :2 ], wallDirs[ gridPos[0],gridPos[1],gridPos[2], :2 ]) ).copy()
        header.AddColumnFloats( filaThickLabel )
        header.AddColumnFloats( filaMassLabel )
        header.AddColumnFloats( wallThickLabel )
        header.AddColumnFloats( wallMassLabel )



# write the output binary halo file
header.AddProgramCommands( options.programOptionsDesc )
halo.writeHaloData( options.outputFile, header, dataIntegers, dataFloats )
#header.PrintValues()

