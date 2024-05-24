#!/usr/bin/python

import sys
import numpy as np
import numpy.fft as fft
from density import DensityHeader, readDensityData, writeDensityData
from MMF import MMFHeader, readMMFData, writeMMFData, getHeaderType
import mathTools
from miscellaneous import throwError



help = """
This program is used to compute the Gaussian filter of a 3D density field.

Usage: %s  inputDensityFile  outputFile  filterRadius  (filter_type) (additional options)

Where the output is also a density/MMF file and the filter radius should be given in Mpc.
The option 'filter_type' =['--gaussian'/'-g','--top-hat'/'-t'] can be used to specify the filter to be used 'gaussian' (DEFAULT) or 'top-hat'.
The additional optins are:
\t'--log' : to filter the logarithm of the density\n""" % sys.argv[0].rsplit('/')[-1]
allowedFiles = ['Density', 'MMF']
filterTypeOptions = { '--gaussian':'gaussian', '-g':'gaussian', '--top-hat':'top-hat', '-t':'top-hat' }
options2 = ['--log']


# read the arguments supplied by the user
filterType = None
logTransform = False
if len(sys.argv)<4:
    print help
    sys.exit(1)
elif len(sys.argv)>=5:
    for i in range(4,len(sys.argv)):
        if sys.argv[i] in filterTypeOptions.keys():
            filterType = filterTypeOptions[ sys.argv[i] ]
        elif sys.argv[i] in options2:
            if sys.argv[i] == '--log':
                print "Detected the '--log' option"
                logTransform = True
        else:
            print 'Unrecognized value for the filter_type program option.'
            print help
            sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
filterRadius = float(sys.argv[3])
if filterType is None: filterType = filterTypeOptions[ '--gaussian' ]
programOptions = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])


# read the density file
headerType = getHeaderType(inputFile)
if headerType not in allowedFiles:
    throwError( "Unrecognized input file type in the program '%s'. Allowed file types are '%s', while the input file is of type: '%s'." % (sys.argv[0],str(allowedFiles),headerType) )

header, data = None, None
if headerType is allowedFiles[0]:
    header, data = readDensityData(inputFile)
elif headerType is allowedFiles[1]:
    header, data = readMMFData(inputFile)
data.shape = header.gridSize



# filter the data
if logTransform:
    print "Taking the logarithm of the input data due to the '--log' option..."
    data = np.log10( data )

data = mathTools.Filter( data, (filterRadius,), header.BoxLength(), filterType=filterType )

if logTransform:
    print "Computing back the density after smoothing the logarithm of the density..."
    data = (10.**data).astype(np.float32)



#write the result to file
header.AddProgramCommands( programOptions )
if headerType is allowedFiles[0]:
    writeDensityData(outputFile,header,data)
elif headerType is allowedFiles[1]:
    writeMMFData(outputFile,header,data)
