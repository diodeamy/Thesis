#!/usr/bin/python

import sys
import numpy as np
from density import DensityHeader, readDensityData
from MMF import MMFHeader, readMMFData, getHeaderType
import analysis
from textFile import readTextFile, writeTextFile
from miscellaneous import throwError, throwWarning


environments = {4:'node', 3:'fila', 2:'wall', 0:'field', 10:'all'}
histRange = (1.e-2, 1.e+4)
noHistBins = 60
extension = '.hist'
desc = "# Contains the density histogram for the '%s' MMF environment of the cosmic web. \n# The file contains: \n#\t 1 = the density bin values\n#\t 2 = the count in each bin\n#\t 3 = the PDF in each bin\n# 4 = the cumulative sum of the PDF. \n# The results were obtained using: %s\n\n"


help = """ The program outputs a file with the histogram of the density in the different environments.
%s  MMF_file  density_file  outputFileRoot
""" % sys.argv[0]

if len(sys.argv) is not 4:
    print "Incorrect use of the program."
    print help
    sys.exit(1)


# get the file names
mmfFile = sys.argv[1]
densityFile = sys.argv[2]
outputFileRoot = sys.argv[3]
programOptions = " ".join(sys.argv)


#read the MMF file
headerType = getHeaderType(mmfFile)
if headerType is not 'MMF':
    throwError( "The first argument supplied to the '%s' program should be the name of an MMF file that gives the cosmic web environments." % sys.argv[0] )
mmfHeader, mmfData = readMMFData(mmfFile)
if mmfHeader.fileType not in [20,21]:
    throwError( "The MMF response file supplied must contain an MMF clean response (file type 20) or the combined MMF clean response for all environments (file type 21)." )
if mmfHeader.fileType==20 and mmfHeader.feature not in [2,3,4]:
    throwError( "The MMF response file supplied must contain the MMF clean response for a given feature ('node'=4, 'fila'=3 or 'wall'=2). But the file contains the feature '%s'!" % str(mmfHeader.Feature()) )

#read the density data
headerType = getHeaderType(densityFile)
if headerType is not 'Density':
    throwError( "The second argument supplied to the '%s' program should be the name of the density file that used to compute the MMF response." % sys.argv[0] )
densityHeader, densityData = readDensityData(densityFile)


###### compute the density histogram for the full sample
print "Computing the histogram for all density values ... "
# compute the histogram for the given density values
histData, histBins = analysis.Histogram(densityData, no_bins=noHistBins, Range=histRange, bin_type="logarithm" )
# get the PDF
pdf = analysis.PDF( histData, histRange, bin_type="logarithm" )
# get the cumulative sum for the histogram data
cumSum = analysis.CumulativeSum(histData, order="descending")
norm = np.float64( cumSum[0] )
cumPDF = cumSum[:] / norm

#write the result to a text file
result = np.column_stack( (histBins, histData, pdf, cumPDF) )
allDesc = "# Contains the density histogram for all density values. \n# The file contains: \n#\t 1 = the density bin values\n#\t 2 = the count in each bin\n#\t 3 = the PDF in each bin\n# 4 = the cumulative sum of the PDF. \n# The results were obtained using: %s\n\n" % programOptions
outputFile = outputFileRoot +'_'+ environments[10] + extension
writeTextFile( outputFile, allDesc, result )


#write in a separate array all the density values in valid pixels
if mmfHeader.fileType == 20:    #if this contains the clean response for a single feature
    print "Computing the histogram for a single MMF environment ... "
    # select all the density values corresponding to the feature in question
    densityFeature = np.extract( mmfData==1, densityData )
    # compute the histogram for the given density values
    histData, histBins = analysis.Histogram(densityFeature, no_bins=noHistBins, Range=histRange, bin_type="logarithm" )
    # get the PDF
    pdf = analysis.PDF( histData, histRange, bin_type="logarithm" )
    # get the cumulative sum for the histogram data
    cumSum = analysis.CumulativeSum(histData, order="descending")
    norm = np.float64( cumSum[0] )
    cumPDF = cumSum[:] / norm
    
    #write the result to a text file
    result = np.column_stack( (histBins, histData, pdf, cumPDF) )
    myDesc = desc % (mmfHeader.Feature(),programOptions)
    outputFile = outputFileRoot +'_'+ environments[mmfHeader.feature] + extension
    writeTextFile( outputFile, myDesc, result )
    
elif mmfHeader.fileType == 21:    #if this contains the clean response for multiple environments
    print "Computing the histogram for all MMF environments ... "
    #loop over features
    results = None
    for i in range(4,1,-1):
        # select all the density values corresponding to the feature in question
        densityFeature = np.extract( mmfData==i, densityData )
        # compute the histogram for the given density values
        histData, histBins = analysis.Histogram(densityFeature, no_bins=noHistBins, Range=histRange, bin_type="logarithm" )
        # get the PDF
        pdf = analysis.PDF( histData, histRange, bin_type="logarithm" )
        # get the cumulative sum for the histogram data
        cumSum = analysis.CumulativeSum(histData, order="descending")
        norm = np.float64( cumSum[0] )
        cumPDF = cumSum[:] / norm
        
        #write the result to a text file
        result = np.column_stack( (histBins, histData, pdf, cumPDF) )
        myDesc = desc % (environments[i],programOptions)
        outputFile = outputFileRoot +'_'+ environments[i] + extension
        writeTextFile( outputFile, myDesc, result )
    
    
    # get the histogram for the remaining data - i.e. field
    densityFeature = np.extract( mmfData==0, densityData )
    histData, histBins = analysis.Histogram(densityFeature, no_bins=noHistBins, Range=histRange, bin_type="logarithm" )
    # get the PDF
    pdf = analysis.PDF( histData, histRange, bin_type="logarithm" )
    # get the cumulative sum for the histogram data
    cumSum = analysis.CumulativeSum(histData, order="descending")
    norm = np.float64( cumSum[0] )
    cumPDF = cumSum[:] / norm
    
    #output the results to a text file
    result = np.column_stack( (histBins, histData, pdf, cumPDF) )
    myDesc = desc % (environments[0],programOptions)
    outputFile = outputFileRoot +'_'+ environments[0] + extension
    writeTextFile( outputFile, myDesc, result )