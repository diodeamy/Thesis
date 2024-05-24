#!/usr/bin/python

import sys
import numpy as np
from density import DensityHeader, readDensityData, writeDensityData
from MMF import MMFHeader, readMMFData, writeMMFData, getHeaderType
from miscellaneous import throwError
from scipy.weave import inline, converters



help = """ Used to downsample a grid quantity by an integer factor.

Usage: %s inputFile outputFile downsampleFactor

    where 'inputFile' can be a density or MMF file
    'downsample' can be a single integer value (downsample by that amount along each direction) or 3 integer values that give the downsample factor accross each axis""" % sys.argv[0]


# Check for the correct number of script parameters
if len(sys.argv)!=4 and len(sys.argv)!=6:
    print help
    sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
downsample = np.zeros( 3, np.uint64 )
if len(sys.argv)==4:
    temp = int(sys.argv[3])
    downsample[0:3] = [temp, temp, temp]
elif len(sys.argv)==6:
    downsample[0:3] = [ int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]) ]




headerType = getHeaderType(inputFile)
if headerType not in ['Density','MMF']:
    throwError( "Unknow file type. This program can be used to downsample grid quantities saved in a density or MMF file type. File type detected '%s'." % headerType )
header = None
data = None
if headerType is 'Density':
    [header, data] = readDensityData(inputFile)
elif headerType is 'MMF':
    [header, data] = readMMFData(inputFile)

# Use C++ code to downsample rapidly the grid data
downsampleCode = """
#line 1000 "downsampleCode"
int n[] = { downsample(0), downsample(1), downsample(2) };
int factor = n[0]*n[1]*n[2];
int newGrid[] = { grid(0)/n[0], grid(1)/n[1], grid(2)/n[2] };
for (int i1=0; i1<newGrid[0]; ++i1)
    for (int i2=0; i2<newGrid[1]; ++i2)
        for (int i3=0; i3<newGrid[2]; ++i3)
        {
            newData(i1,i2,i3) = 0.;
            for (int j1=0; j1<n[0]; ++j1)
                for (int j2=0; j2<n[1]; ++j2)
                    for (int j3=0; j3<n[2]; ++j3)
                        newData(i1,i2,i3) += data( i1*n[0]+j1, i2*n[1]+j2, i3*n[2]+j3 );
            newData(i1,i2,i3) /= factor;
        }
"""
grid = header.gridSize.astype('i4')
newGrid = (grid/downsample).astype('i8')
#newGrid = [ grid[i]/downsample[i] for i in range(3) ]
print "Downsampling the '%s' data from a %s grid to a %s grid ... " % (headerType,str(grid),str(newGrid)),
sys.stdout.flush()  # flush the print statements
data.shape = grid
newData = np.zeros( newGrid, data.dtype )
inline( downsampleCode, ['data','newData','grid','downsample'], type_converters=converters.blitz )
print "Done."

# Change the grid values in the header
header.gridSize = newGrid.astype(np.uint64)
header.totalGrid = np.uint64( newGrid[0]*newGrid[1]*newGrid[2] )

# Write the data to the file
if headerType is 'Density':
    writeDensityData(outputFile,header,newData)
elif headerType is 'MMF':
    writeMMFData(outputFile,header,newData)



