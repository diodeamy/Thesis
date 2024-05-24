#!/usr/bin/python

import sys
import optparse
import numpy as np
from density import DensityHeader, readDensityData
from MMF import MMFHeader, readMMFData, getHeaderType
from AHF import AHFHalos, readAHFHalos, halosInMassRange, halosInMask
from miscellaneous import throwError, throwWarning



help = """
Use this program the get some overall properties of the MMF environments as: volume fraction, mass fraction, halo content, etc ...

Usage: %s  MMFcleanFile  densityFile  outputFileName  (optional:  AHFhaloFile  )
""" % sys.argv[0]

rhoCritical = 27.7538e+10    # the critical density in units of (M0/h) / (Mpc/h)^3


#read the user options
if len(sys.argv)<4:
    print help
    sys.exit(1)

mmfFile = sys.argv[1]
densityFile = sys.argv[2]
outputFile = sys.argv[3]
ahfFile = None
if len(sys.argv)>4:
    ahfFile = sys.argv[4]
programOptions = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])


# read the MMF file
headerType = getHeaderType(mmfFile)
if headerType is not 'MMF':
    throwError( "Unrecognized type for the first file name supplied to the program '%s'. Program expects to receive a clean MMF response file, but it received a: '%s'." % (sys.argv[0],headerType) )
mmfHeader, mmfData = readMMFData(mmfFile)
if mmfHeader.fileType!=21: # this test if file is a clean response
    throwError( "File '%s' should contain the combined MMF clean response (fileType '21'). But it does not contain a combined MMF clean response, but it contains '%s'." % (mmfFile,header.FileType()) )


# read the density file
headerType = getHeaderType(densityFile)
if headerType is not 'Density':
    throwError( "Unrecognized type for the second file name supplied to the program '%s'. Program expoects to receive a density file, but it received a: '%s'." % (sys.argv[0],headerType) )
densityHeader, densityData = readDensityData(densityFile)
densityData = densityData.astype(np.float64)
cellVolume = (densityHeader.box[1]-densityHeader.box[0]) * (densityHeader.box[3]-densityHeader.box[2]) * (densityHeader.box[5]-densityHeader.box[4]) / densityHeader.totalGrid  # volume of a grid cell
OmegaDM = densityHeader.Omega0    # Omega DM for the given simulation


# read the AHF halo
haloData = None
pos, part, mass, radius = None, None, None, None
if ahfFile:
    haloData = readAHFHalos(ahfFile)
    pos = ".array['pos']"      # 'halosData' array indices giving the position of the halos
    part = ".array['npart']"   # 'halosData' array indices giving the number of particles in each halo
    mass = ".array['Mvir']"    # 'halos' array indices giving the mass of the halos
    radius = ".array['Rvir']"  # 'halos' array indices giving the virial radius of the halos



# get the volume fraction and mass in each environment
results = """# This file contains a series of properties for the MMF environments. \n# Each property is given on two lines: \n#\t 1st line = gives name and short description of property \n#\t 2nd line = gives the property values for: nodes, filaments, walls and fields. \n# The file was obtained using the commands: %s\n\n""" % programOptions
# volume fraction in each environment
density = densityData.astype('f8')
nodeDensity = density[mmfData==4]
filaDensity = density[mmfData==3]
wallDensity = density[mmfData==2]
fieldDensity = density[mmfData==0]
grid = float(mmfHeader.totalGrid)/100.
results += """\n# Volume fraction of environments in %%: \n %f\t %f\t %f\t %f\n""" % (nodeDensity.size/grid, filaDensity.size/grid, wallDensity.size/grid, fieldDensity.size/grid)
#mass fraction in each environment
totalMass = density.sum()/100.
results += """\n# Mass fraction of environments in %%: \n %f\t %f\t %f\t %f\n""" % (nodeDensity.sum()/totalMass, filaDensity.sum()/totalMass, wallDensity.sum()/totalMass, fieldDensity.sum()/totalMass)
# properties of the density distribution in each environment
results += """\n# Properties of the density distribution in the environments:\n"""
results += """# average density: \t%f\t %f\t %f\t %f\n""" % (np.mean(nodeDensity),np.mean(filaDensity), np.mean(wallDensity), np.mean(fieldDensity))
results += """# median density: \t%f\t %f\t %f\t %f\n""" % (np.median(nodeDensity),np.median(filaDensity), np.median(wallDensity), np.median(fieldDensity))
results += """# standard deviation density: \t%f\t %f\t %f\t %f\n""" % (np.std(nodeDensity),np.std(filaDensity), np.std(wallDensity), np.std(fieldDensity))
results += """# standard deviation density = sqrt( sum (x-x_mean)**2 /N )\n"""


# write the number and mass fraction of the halos in the different environments
if haloData:
    # get the halos in each MMF environment
    print "\nFinding the halos in each environment ... "
    mmfData.shape = mmfHeader.gridSize
    validMask = mmfData==4  # find node halos
    nodeHalos, restHalos = halosInMask( haloData, eval('haloData'+pos), mmfHeader.box, validMask) # returns the halos in mask and the ones outside the mask
    validMask = mmfData==3  # find filament halos
    filaHalos, restHalos = halosInMask( haloData, eval('haloData'+pos), mmfHeader.box, validMask)
    validMask = mmfData==2  # find wall halos
    wallHalos, restHalos = halosInMask( haloData, eval('haloData'+pos), mmfHeader.box, validMask)
    validMask = mmfData==0  # find field halos
    fieldHalos, restHalos = halosInMask( haloData, eval('haloData'+pos), mmfHeader.box, validMask)
    mmfData.shape = -1
    
    
    noHalos = float(haloData.noHalos)/100.
    results += """\n# Total number of halos in environments (%i halos in total) + following line %% of that total number of halos:  \n %i\t %i\t %i\t %i\n %f\t %f\t %f\t %f\n""" % (haloData.noHalos, nodeHalos.noHalos, filaHalos.noHalos, wallHalos.noHalos, fieldHalos.noHalos, nodeHalos.noHalos/noHalos, filaHalos.noHalos/noHalos, wallHalos.noHalos/noHalos, fieldHalos.noHalos/noHalos)
    # Get statistics of the mass in halos
    gridCellMass = rhoCritical * OmegaDM * cellVolume  # the mass inside a grid cell
    totalMass = totalMass * 100. * gridCellMass
    totalHaloMass = ( eval('haloData'+mass).astype(np.float64) ).sum()
    nodeHaloMass = ( eval('nodeHalos'+mass).astype(np.float64) ).sum()
    filaHaloMass = ( eval('filaHalos'+mass).astype(np.float64) ).sum()
    wallHaloMass = ( eval('wallHalos'+mass).astype(np.float64) ).sum()
    fieldHaloMass = ( eval('fieldHalos'+mass).astype(np.float64) ).sum()
    results += """\n# The mass statistics accross halos (mass in halos=%.2f %% of total mass in the simulation):\n""" % (totalHaloMass/totalMass*100.)
    results += """# Halo mass distribution with respect to total mass: \t%f\t %f\t %f\t %f\n""" % (nodeHaloMass/totalMass*100.,filaHaloMass/totalMass*100.,wallHaloMass/totalMass*100.,fieldHaloMass/totalMass*100.)
    results += """# Halo mass distribution with respect to total mass in halos: \t%f\t %f\t %f\t %f\n""" % (nodeHaloMass/totalHaloMass*100.,filaHaloMass/totalHaloMass*100.,wallHaloMass/totalHaloMass*100.,fieldHaloMass/totalHaloMass*100.)



# write the results to a file
output = open(outputFile,'w')
output.write(results)
output.close()

