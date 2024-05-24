#!/usr/bin/python

import sys
import numpy as np
import density
import textFile
from miscellaneous import throwError, throwWarning




programName = sys.argv[0].rsplit('/')[-1]
programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])

help = '''Use this program to relabel the Watershed Voids according to their volume, with the largest void starting with 1. It outputs the new void indices on a grid as well as a text file with the statistics of the voids.
USAGE:  %s  inputDensityFile  outputFile
''' % programName


if len(sys.argv)!=3:
    print help
    sys.exit(1)
inputFile  = sys.argv[1] 
outputFile = sys.argv[2] 


# read the input file
header, data = density.readDensityData( inputFile )
#data.shape = header.gridSize


# loop over each void index and find the size of each void
print "Computing the size of each void ..."
noVoids = data.max()+1
d2 = data[ data>=0 ]
voidSize = np.bincount( d2 ) # returns the number of points with that given bin value --- it has noVoids entries
totalVolume = np.sum( voidSize.astype(np.float64) )
voidVolumeFraction = totalVolume / data.size * 100.
print "Volume fraction in WVF voids: %.2f%%" % voidVolumeFraction


# Relabel the voids according to their volume
print "Relabeling the voids according to their volume ..."
voidSize[0] = data.size
order = voidSize.argsort()[::-1]
voidSize = voidSize[order]
# need to create the matrix which gives the index where the old void indices point in the volume ordered void set
newIndex = np.empty( (1,order.size), np.int32 )
newIndex[:,order] = np.arange( order.size )
newIndex.shape = -1
select = data>=0
d2 = data.copy()
d2[select] = newIndex[data[select]]
print d2.min(), '  ,  ', d2.max()



# compute the physical size of each void
cellVolume = header.CellVolume()
noValidVoids = np.sum( voidSize[1:]>0 )
print "Found %i valid voids out of %i total voids." % ( noValidVoids, noVoids-1 )
stats = np.zeros( (noValidVoids,4), np.float64 )
stats[:,0] = 1 + np.arange(noValidVoids)    # void label
stats[:,1] = voidSize[1:noValidVoids+1]     # number of cells for each void
stats[:,2] = stats[:,1] * cellVolume        # void volume in (Mpc/h)^3
stats[:,3] = (stats[:,2] * 3./4./np.pi )**(1./3.)       # void radius for spherically symmetric voids


# output the new void indices
header.AddProgramCommands( programOptionsDesc )
density.writeDensityData( outputFile, header, d2 )


# write the void statistics to a text file
desc = '''# Gives the statistics of the Watershed Void Finder (WVF) analysis.
# WVF assigned %.2f of the grid cells to %i distinct voids with labels from 1 to %i.
# The results were obtained using the command: %s 

# (1)void_index   (2)number_grid_cells_in_void   (3)void_volume_in_(Mpc/h)**3    (4)void_radius_in_Mpc/h_for_spherical_voids
''' % ( voidVolumeFraction, noValidVoids, noValidVoids, programOptionsDesc )
textFile.writeTextFile( outputFile+'.stats', desc, stats )
