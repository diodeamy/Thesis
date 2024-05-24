#!/usr/bin/python

import sys
import numpy as np
import MMF
import density
import analysis
import textFile
from miscellaneous import throwError, throwWarning, BoxLength, massInCell


programName = sys.argv[0].rsplit('/')[-1]
help = """ 
This program computes the properties of Cosmic Web filaments and walls.

Usage: %s  inputMMF_cleanFile  inputMMF_directionsFile  input_densityFile  outputFile (options)

where 'options' can be one of the options ['--fila','--wall'] - these are used only when the input file contains all the Cosmic Web environments.
""" % programName




if len(sys.argv)<5:
    print help
    sys.exit(1)

cleanFile = sys.argv[1]
directionsFile = sys.argv[2]
densityFile = sys.argv[3]
outputFile = sys.argv[4]
programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])
options = { '--fila':[4,2], '--wall':[4,3] }
feature = None
if len(sys.argv)>5:
    if sys.argv[5] in options.keys():
        feature = options[sys.argv[5]]
    else:
        throwError( "Unknown option '%s'. Available options are '%s'." % (sys.argv[5],options.keys()) )


# load the clean MMF response file
headerType = MMF.getHeaderType(cleanFile)
if headerType is not 'MMF': throwError( "The first argument supplied to the '%s' program should be the name of a filament clean MMF file." % programName )
cleanHeader, cleanData = MMF.readMMFData(cleanFile)
if cleanHeader.feature not in [MMF.MMF_FILAMENT,MMF.MMF_WALL,MMF.MMF_ALL]:
    throwError( "The first argument supplied to the '%s' program should be the name of a clean MMF file storing filament/wall data. The input file '%s' stores '%s' data." % (programName,cleanFile,cleanHeader.Feature()) )
if cleanHeader.fileType not in [MMF.MMF_CLEAN_RESPONSE,MMF.MMF_CLEAN_RESPONSE_COMBINED]:
    throwError( "The first argument supplied to the '%s' program should be the name of a clean MMF file. The input file '%s' is a '%s' file type." % (programName,cleanFile,cleanHeader.FileType()) )
if cleanHeader.feature==MMF.MMF_ALL:
    if not feature:
        throwError( "When the input file is a file containing all the Cosmic Web environments you must specify an additional option '--fila'/'--wall'." )
    cleanData[ cleanData==feature[0] + cleanData==feature[1] ] = 0
    cleanData[ cleanData!=0 ] = 1
    cleanHeader.feature = (3,2)[int(feature[1]==2)]
cleanData.shape = cleanHeader.gridSize



# load the MMF directions file
headerType = MMF.getHeaderType(directionsFile)
if headerType is not 'MMF': throwError( "The thrid argument supplied to the '%s' program should be the name of a MMF directions file." % programName )
directionsHeader, directionsData = MMF.readMMFData(directionsFile)
if directionsHeader.feature not in [MMF.MMF_FILAMENT,MMF.MMF_WALL]:
    throwError( "The thrid argument supplied to the '%s' program should be the name of a MMF directions file storing filament/wall data. The input file '%s' stores '%s' data." % (programName,directionsFile,directionsHeader.Feature()) )
if directionsHeader.fileType!=MMF.MMF_DIRECTIONS:
    throwError( "The thrid argument supplied to the '%s' program should be the name of a MMF directions file. The input file '%s' is a '%s' file type." % (programName,directionsFile,directionsHeader.FileType()) )
directionsData.shape = (-1,3)



# read the density file
headerType = MMF.getHeaderType(densityFile)
if headerType is not 'Density':
    throwError( "Unrecognized type for the second file name supplied to the program '%s'. Program expects to receive a density file, but it received a: '%s'." % (programName,headerType) )
densityHeader, densityData = density.readDensityData(densityFile)
densityData.shape = densityHeader.gridSize



# get the different objects in the input file
densityData = densityData * massInCell(densityHeader)
sizeData, massData, averageData = MMF.MMFFeatureProperties( cleanData, directionsData, densityData, Feature=cleanHeader.feature, Radius=1., BoxLength=BoxLength(densityHeader.box) )
print np.sum(sizeData[:,1]*sizeData[:,0])/np.sum(sizeData[:,1])



# output the results - average properties
feature = ("filaments","walls")[int(cleanHeader.feature==3)]
desc = """This file contains average properties for Cosmic Web %s.\n""" % (feature)
if feature=="filaments":
    desc = """%sAverage cross-section diameter: %f (Mpc/h)""" % (desc,averageData[0])
    desc = """%sAverage cross-section area: %f (Mpc/h)^2""" % (desc,averageData[0]*averageData[0]*3.14/4.)
    desc = """%sAverage length: %f Mpc/h per (Mpc/h)^3""" % (desc,averageData[1])
if feature=="walls":
    desc = """%sAverage wall height: %f (Mpc/h)""" % (desc,averageData[0])
    desc = """%sAverage wall area: %f (Mpc/h)^2 per (Mpc/h)^3""" % (desc,averageData[1])
output = open(outputFile+'.averageProperties','w')
output.write(desc)
output.close()


# output the results - cross-section diameter for filaments and plane height for walls
desc = """#This file contains the filament diameter / wall height for the Cosmic Web environments.\n#The 3 columns contain:\n#\t1st column = the value of the diameter for filamnets / heigh for walls\n#\t2nd column = the volume fraction in filaments/walls with the given size\n#\t3rd column = the PDF of the given distribution\n#The results were obtained using the program options:  %s\n\n""" % programOptionsDesc
sizePDF = analysis.PDF( sizeData[:,1], (sizeData[0,0],sizeData[-1,0]), "linear" )
results = np.column_stack( (sizeData,sizePDF) )
textFile.writeTextFile( outputFile+'.size', desc, results )


# output the results - linear mass density for filaments and surface mass density for walls
desc = """#This file contains linear mass density for filaments / surface mass density for walls for the Cosmic Web environments.\n#The 3 columns contain:\n#\t1st column = the value of the mass density bin\n#\t2nd column = the volume fraction in filaments/walls with the given mass density\n#\t3rd column = the PDF of the given distribution\n#The results were obtained using the program options:  %s\n\n""" % programOptionsDesc
massPDF = analysis.PDF( massData[:,1], (massData[0,0],massData[-1,0]), "logarithm" )
results = np.column_stack( (massData,massPDF) )
textFile.writeTextFile( outputFile+'.mass', desc, results )
