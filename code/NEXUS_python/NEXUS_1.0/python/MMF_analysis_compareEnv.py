#!/usr/bin/python2.6

import sys
import numpy as np
from MMF import MMFHeader, readMMFData, writeMMFData
import MMF
import density
from miscellaneous import throwError
from scipy.weave import inline, converters


programName = sys.argv[0].rsplit('/')[-1]
help = """ Used to compare between two different MMF runs for all the environments: node, filaments and walls.

Usage:  %s  all_env_MMF_clean_file_1  all_env_MMF_clean_file_2  density_file  output_file

The input files must be two MMF clean response files that keep information for all environments.
"""  % sys.argv[0]


if len(sys.argv)!=5:
    print help
    sys.exit(1)

allClean_1 = sys.argv[1]
allClean_2 = sys.argv[2]
densityFile = sys.argv[3]
outputFile = sys.argv[4]
programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])
environments = {4:'node', 3:'fila', 2:'wall'}


# read the first clean response
header1, data1 = readMMFData(allClean_1)
if header1.feature!=MMF.MMF_ALL:   # test that this is a wall file
    throwError( "File '%s' should contain the MMF clean response for all environments. But it contain an MMF response for the '%s' feature." % (allClean_1,header1.Feature()) )
if header1.fileType!=MMF.MMF_CLEAN_RESPONSE_COMBINED: # this test if file is a clean response
    throwError( "File '%s' should contain the MMF combined clean response. But it does not contain an MMF clean response, but it contains '%s'." % (allClean_1,header1.FileType()) )

# read the second clean response
header2, data2 = readMMFData(allClean_2)
if header2.feature!=MMF.MMF_ALL:   # test that this is a wall file
    throwError( "File '%s' should contain the MMF clean response for all environments. But it contain an MMF response for the '%s' feature." % (allClean_2,header2.Feature()) )
if header2.fileType!=MMF.MMF_CLEAN_RESPONSE_COMBINED: # this test if file is a clean response
    throwError( "File '%s' should contain the MMF clean response. But it does not contain an MMF clean response, but it contains '%s'." % (allClean_2,header2.FileType()) )
if (header1.gridSize!=header2.gridSize).any():
    throwError( "You are trying to compare two MMF responses on different grids. The file '%s' has grid size %s while the second file '%s' has grid size %s." % (allClean_1,str(header1.gridSize),allClean_2,str(header2.gridSize)) )
    

# read the density file
headerDensity, density = density.readDensityData(densityFile)
if (header1.gridSize!=headerDensity.gridSize).any():
    throwError( "You are trying to compare two objects on different grids. The MMF clean response file '%s' has grid size %s while the density file '%s' has grid size %s." % (allClean_1,str(header1.gridSize),densityFile,str(headerDensity.gridSize)) )


# Loop over the environments and get statistics of the comparison between the two
results = "# Comparison of the MMF results in 1st file='%s' and the 2nd file='%s'.\n# The program was obtained using the commad: %s\n\n" % (allClean_1,allClean_2,programOptionsDesc)
totalSize = float(data1.size)
totalMass = density.sum()
for i in [4,3,2]:
    env1, env2 = data1==i, data2==i
    envBoth = env1 * env2
    envOnly1, envOnly2 = env1-envBoth, env2-envBoth
    size1, size2, sizeBoth, sizeOnly1, sizeOnly2 = env1.sum(), env2.sum(), envBoth.sum(), envOnly1.sum(), envOnly2.sum()
    den1, denBoth, denOnly1, denOnly2 = density[env1], density[envBoth], density[envOnly1], density[envOnly2]
    mass1, massBoth, massOnly1, massOnly2 = den1.sum(), denBoth.sum(), denOnly1.sum(), denOnly2.sum()
    avgDen1, avgDenBoth, avgDenOnly1, avgDenOnly2 = mass1/size1, massBoth/sizeBoth, massOnly1/sizeOnly1, massOnly2/sizeOnly2
    
    results += '\n# Comparison of %s environments:\n' % environments[i]
    results += '# Comparison of volumes: totalVolume1(%% of box volume)  commonVolume(%% of 1st volume)  volumeOnly_1(%% of 1st volume)  volumeOnly_2(%% of 1st volume)\n'
    results += '\t %f\t %f\t %f\t %f\n' % ( 100.*size1/totalSize, 100.*sizeBoth/float(size1), 100.*sizeOnly1/float(size1), 100.*sizeOnly2/float(size1) )
    
    results += '# Comparison of masses: totalMass1(%% of box mass)  commonMass(%% of 1st value)  massOnly_1(%% of 1st value)  massOnly_2(%% of 1st value)\n'
    results += '\t %f\t %f\t %f\t %f\n' % ( 100.*mass1/totalMass, 100.*massBoth/mass1, 100.*massOnly1/mass1, 100.*massOnly2/mass1 )
    
    results += '# Comparison of avergae densities: avgDen1  commonAvgDen  avgDenOnly_1  avgDenOnly_2 \n'
    results += '\t %f\t %f\t %f\t %f\n' % ( avgDen1, avgDenBoth, avgDenOnly1, avgDenOnly2 )


out = open( outputFile, 'w' )
out.write( results )
out.close()


