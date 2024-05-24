#!/usr/bin/python

import sys
import optparse
import numpy as np
import numpy.fft as fft
from density import DensityHeader, readDensityData, writeDensityData
from MMF import MMFHeader, readMMFData, writeMMFData, getHeaderType
import MMF
import mathTools
from miscellaneous import throwError, throwWarning





def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputFile  outputFile  [options] <arg1> ... <argN>', description='Use this file to take a slice or output a certain region through a density/NEXUS binary file.')
    parser.add_option('--axis', action="store", type='int', default=None, dest="axis", metavar='label', help="Give the axis along which to take the slice (0=x, 1=y and 2=z). By giving this option you specify that the program will take a slice of the data.")
    parser.add_option('-c','--center', action="store", type='float', default=None, dest="center", metavar='value', help="Give the center of the slice along the axis specified using the '--axis' option.")
    parser.add_option('-t','--thickness', action="store", type='float', default=None, dest="thickness", metavar='value', help="Give the thickness of the slice along the axis coordinate. The slice is takes as 'center-thicknes/2' to 'center+thicknes/2'.")
    parser.add_option('-s','--selection', action="store", nargs=6, type='float', default=None, dest="selection", metavar='coords', help="Select a region in the given data set. You can specify this region using 'xMin xMax yMin yMax zMIn zMax' coordinates. The output is the part of the grid that fits completely inside the above coordinates. The program assumes periodic conditions, so it will automatically periodically translate the grid to the desired coordinates.")
    
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    if len(args)<2:
        throwError( "Imcomplete number of program arguments. The program needs at least 2 arguments giving: inputFile outputFile. Use the '--help' command for additional details." )
    options.inputFile = args[0]
    options.outputFile = args[1]
    
    # check the options for the slice
    countTrue = int(options.axis!=None) + int(options.center!=None) + int(options.thickness!=None)
    if options.axis!=None and 3!=countTrue:
        throwError( "To take a slice of the input data you need to specify all of the '--axis', '--center' and '--thickness' options." )
    options.subregion = np.array( options.selection )
    
    return options





# get the program options
options = programOptions()
options.programOptionsDesc = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])
programName  = sys.argv[0].rsplit('/')[-1]
allowedFiles = [ 'Density', 'MMF' ]



# read the density file
headerType = getHeaderType( options.inputFile )
if headerType not in allowedFiles:
    throwError( "Unrecognized input file type in the program '%s'. Allowed file types are '%s', while the input file is of type: '%s'." % (programName,str(allowedFiles),headerType) )

header, data = None, None
if headerType is allowedFiles[0]:
    header, data = readDensityData(options.inputFile)
elif headerType is allowedFiles[1]:
    header, data = readMMFData(options.inputFile)
comp = header.DataComponents()
data.shape = header.gridSize[0], header.gridSize[1],  header.gridSize[2], comp
result = None



# code for taking a slice of the input data
if options.axis!=None:
    # do some error check
    axis, center, thickness = int(options.axis), float(options.center), float(options.thickness)
    if axis not in range(0,3):
        throwError( "Unrecognized value for the '--axis' option. This option can take values in the list '%s', but you inserted the value %i." % (str(range(0.3)),axis) )
    if center<header.box[2*axis] or center>header.box[2*axis+1]:
        throwError( "The value for the '--center' option is outside the boundaries of the data box. You inserted the value %f for the center, while the box on axis %i extends from %f to %f." % (center,axis,header.box[2*axis],header.box[2*axis+1]) )
    if thickness>(header.box[2*axis+1]-header.box[2*axis]):
        throwError( "The value for the '--thickness' option is larger than the length of the data box. You inserted the value %f for the thickness, while the box has a length of %f on axis %i." % (thickness,header.box[2*axis+1]-header.box[2*axis],axis) )


    # compute the slice
    dx = (header.box[2*axis+1]-header.box[2*axis]) / header.gridSize[axis]
    minPos, maxPos = center-thickness/2, center+thickness/2
    minIndex, maxIndex = int( round( (minPos-header.box[2*axis]) / dx) ), int( round( (maxPos-header.box[2*axis]) / dx) )
    if minIndex<0:
        throwWarning( "The lower limit of the slice extends beyond the boundaries of the box. The slice will be take only until the limits of the box." )
        minIndex = 0
        minPos = header.box[2*axis]
    if maxIndex>= header.gridSize[axis]:
        throwWarning( "The upper limit of the slice extends beyond the boundaries of the box. The slice will be take only until the limits of the box." )
        maxIndex = header.gridSize[axis]
        maxPos = header.box[2*axis+1]
    print "Taking a slice along the %s-axis between positions %s to %s ..." % (['x','y','z'][axis],minPos,maxPos)
    command = '[minIndex:maxIndex,:,:,:]'
    if axis==1: command = '[:,minIndex:maxIndex,:,:]'
    elif axis==2: command = '[:,:,minIndex:maxIndex,:]'
    result = eval( 'data' + command ).sum(axis=axis)
    if headerType is allowedFiles[0]:
        result /= maxIndex - minIndex   #normalize the density
    elif headerType is allowedFiles[1]:
        if header.fileType==MMF.MMF_CLEAN_RESPONSE or header.fileType==MMF.MMF_CLEAN_RESPONSE_COMBINED:
            result[ result>0 ] = 1
            result = result.astype(np.int16)
        else:
            result /= maxIndex - minIndex
    
    # make the changes in the output header
    header.gridSize[axis] = 1
    header.totalGrid = np.int64( header.gridSize[0]*header.gridSize[1]*header.gridSize[2] )
    header.box[2*axis:2*(axis+1)] = minPos, maxPos



# code for outputing a region of the input data
if options.selection!=None:
    #do some error checks
    for i in range(3):
        if options.subregion[2*i]>=options.subregion[2*i+1]:
            throwError( "The coordinates of the subbox along coordinate %s are invalid. The lower diemsnion of the box=%s is higher than the upper diemnsion of the box=%s." % (['x','y','z'][axis],options.subregion[2*i],options.subregion[2*i+1]) )
        p1, p2 = -2.*header.box[2*i+1]+header.box[2*i], 2.*header.box[2*i+1]-header.box[2*i]
        if options.subregion[2*i]>=p2 or options.subregion[2*i+1]<p1 :
            throwError( "The coordinates of the subbox along coordinate %s are invalid. There is no intersection between the data box and the region of interest. Data box=(%s,%s), while region of interest=(%s,%s)" % (['x','y','z'][axis],header.box[2*i],header.box[2*i+1],options.subregion[2*i],options.subregion[2*i+1]) )
    
    # select the grid inside the subregion
    dx = (header.box[1::2]-header.box[::2]) / header.gridSize[:]
    iMin = np.round( (options.subregion[::2]-header.box[::2])/dx[:] ).astype(np.int32)
    iMax = np.round( (options.subregion[1::2]-header.box[::2])/dx[:] ).astype(np.int32 )
    for i in range(3):
        if iMin[i]<-1*header.gridSize[i]:
            iMin[i] = -1*header.gridSize[i]
        if iMax[i]>2*header.gridSize[i]:
            iMax[i] = header.gridSize[i]
    xMin = header.box[0::2] + iMin*dx
    xMax = header.box[0::2] + iMax*dx
    newGrid = (iMax - iMin).astype(np.int64)
    print "The program will output the region [ (%s, %s),  (%s, %s),  (%s, %s) ] using a grid of %i x %i x %i cells ..." % (xMin[0],xMax[0],xMin[1],xMax[1],xMin[2],xMax[2],newGrid[0],newGrid[1],newGrid[2])
    
    # select the region along the x-axis
    result = np.empty( (newGrid[0], header.gridSize[1], header.gridSize[2], comp), data.dtype )
    index = 0
    i1, i2 = iMin[index], iMax[index]
    j1, j2 = int(i1>=0), int(i2<header.gridSize[index])
    k1 = [ 0, i1 ][j1]
    k2 = [ header.gridSize[index], i2 ][j2]
    result[    0:k1-i1,:,:,:] = data[header.gridSize[index]+i1:header.gridSize[index]-k1,:,:,:]
    result[k1-i1:k2-i1,:,:,:] = data[k1:k2,:,:,:]
    result[k2-i1:i2-i1,:,:,:] = data[k2-header.gridSize[index]:i2-header.gridSize[index],:,:,:]
    
    # select the region along the y-axis
    tempRes = np.empty( (newGrid[0], newGrid[1], header.gridSize[2], comp), data.dtype )
    index = 1
    i1, i2 = iMin[index], iMax[index]
    j1, j2 = int(i1>=0), int(i2<header.gridSize[index])
    k1 = [ 0, i1 ][j1]
    k2 = [ header.gridSize[index], i2 ][j2]
    tempRes[:,    0:k1-i1,:,:] = result[:,header.gridSize[index]+i1:header.gridSize[index]-k1,:,:]
    tempRes[:,k1-i1:k2-i1,:,:] = result[:,k1:k2,:,:]
    tempRes[:,k2-i1:i2-i1,:,:] = result[:,k2-header.gridSize[index]:i2-header.gridSize[index],:,:]

    # select the region along the z-axis
    result = np.empty( (newGrid[0], newGrid[1], newGrid[2], comp), data.dtype )
    index = 2
    i1, i2 = iMin[index], iMax[index]
    j1, j2 = int(i1>=0), int(i2<header.gridSize[index])
    k1 = [ 0, i1 ][j1]
    k2 = [ header.gridSize[index], i2 ][j2]
    result[:,:,    0:k1-i1,:] = tempRes[:,:,header.gridSize[index]+i1:header.gridSize[index]-k1,:]
    result[:,:,k1-i1:k2-i1,:] = tempRes[:,:,k1:k2,:]
    result[:,:,k2-i1:i2-i1,:] = tempRes[:,:,k2-header.gridSize[index]:i2-header.gridSize[index],:]
    
    #set output header data
    header.gridSize[:] = newGrid[:]
    header.totalGrid = np.int64( header.gridSize[0]*header.gridSize[1]*header.gridSize[2] )
    header.box[::2] = xMin
    header.box[1::2] = xMax



# output the data
header.AddProgramCommands( options.programOptionsDesc )
if headerType is allowedFiles[0]:
    writeDensityData(options.outputFile,header,result)
elif headerType is allowedFiles[1]:
    writeMMFData(options.outputFile,header,result)
