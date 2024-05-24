#!/usr/bin/python

import sys
import optparse
import numpy as np
import density
import MMF as NEXUS
from miscellaneous import throwError





def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputFile outputFile [options] <arg1> ... <argN>', description='Use this file to rotate a scalar on a grid to a new grid in a different coordinate system. First the z-rotation is done, than the y-rotation and at the end the x-rotation.')
    parser.add_option('-x','--xRotation', action="store", type='float', default=None, dest="xRotation", metavar='angle', help="Give the angle in degrees for a rotation along the x-axis.")
    parser.add_option('-y','--yRotation', action="store", type='float', default=None, dest="yRotation", metavar='angle', help="Give the angle in degrees for a rotation along the y-axis.")
    parser.add_option('-z','--zRotation', action="store", type='float', default=None, dest="zRotation", metavar='angle', help="Give the angle in degrees for a rotation along the z-axis.")
    parser.add_option('-v','--vector', action="store", nargs=3, type='float', default=None, dest="vector", metavar='value', help="Give a vector as its 3 components: vx, vy and vz. The program will compute the rotation angles and rotate the grid such that this vector, after rotation, will lie along the z-axis.")
    parser.add_option('-c','--center', action="store", nargs=3, type='float', default=None, dest="center", metavar='value', help="Give the center of rotation as its 3 components: x, y and z. The program will rotate the grid around this center using the rotation coordinates specified.")
    parser.add_option('-n','--noSamples', action="store", type='int', default=3, dest="noSamples", metavar='value', help="Give the number of samples used along each axis to compute the rotated field. The total number of used samples is 'noSamples'^3.")
    
    
    options, args = parser.parse_args()
    if options.vector and (options.xRotation or options.yRotation or options.zRotation):
        throwError( "You cannot specify an additional x,y or z-roattion when using the '--vector' option." )
    if (options.vector or (options.xRotation or options.yRotation or options.zRotation)) and options.center==None:
        throwError( "You need to specify a rotation center using the '--center' option" )
    if len(args)<2:
        throwError( "You must supply at least two arguments to the program giving the input and output files." )
    
    options.inputFilename = args[0]
    options.outputFilename = args[1]
    
    return options
    


# function used to rotate around a given axis
conversion = np.pi / 180.
def xAxisRotation(pos, angle):
    """ Roates the positions around the x-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp     = cosA * pos[:,1] - sinA * pos[:,2]
    pos[:,2] = sinA * pos[:,1] + cosA * pos[:,2]
    pos[:,1] = temp

def yAxisRotation(pos, angle):
    """ Roates the positions around the y-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp     =  cosA * pos[:,0] + sinA * pos[:,2]
    pos[:,2] = -sinA * pos[:,0] + cosA * pos[:,2]
    pos[:,0] = temp
    
def zAxisRotation(pos, angle):
    """ Roates the positions around the y-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp     = cosA * pos[:,0] - sinA * pos[:,1]
    pos[:,1] = sinA * pos[:,0] + cosA * pos[:,1]
    pos[:,0] = temp


def computeRotationAngle(vector):
    """ Computes the rotation angles along the x- and y-directions such that the input vector, after rotation, lies along the z-direction. """
    tanY = -vector[0] / vector[2]
    angleY = np.arctan( tanY )
    tanX = vector[1] / ( -vector[0] * np.sin(angleY) + vector[2] * np.cos(angleY) )
    angleX = np.arctan( tanX )
    return ( angleX/conversion, angleY/conversion )






# read program options
options = programOptions()
options.programOptionsDesc = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])
programName  = sys.argv[0].rsplit('/')[-1]
allowedFiles = ['Density', 'MMF']



# read the input file
headerType = NEXUS.getHeaderType( options.inputFilename )
if headerType not in allowedFiles:
    throwError( "Unrecognized input file type in the program '%s'. Allowed file types are '%s', while the input file is of type: '%s'." % (programName,str(allowedFiles),headerType) )

header, data = None, None
if headerType is allowedFiles[0]:
    header, data = density.readDensityData( options.inputFilename )
elif headerType is allowedFiles[1]:
    header, data = NEXUS.readMMFData( options.inputFilename )




# if the user gave the '--vector' option, than compute the rotation angles along the x and y directions
if options.vector:
    angleX, angleY = computeRotationAngle( options.vector )
    options.xRotation = angleX
    options.yRotation = angleY
    options.zRotation = None


# do the requested rotations
doRotation = options.xRotation!=None or options.yRotation!=None or options.zRotation!=None
result = None
if not doRotation:
    result = data
else:
    totalSize = header.totalGrid
    result = np.zeros( totalSize, np.float32 )
    data.shape = header.gridSize[0], header.gridSize[1], header.gridSize[2]
    
    # get the coordinates of the input and output grid centers
    g = np.mgrid[ 0:header.gridSize[0], 0:header.gridSize[1], 0:header.gridSize[2] ]
    offset = header.box[::2].copy()
    boxLength = header.box[1::2] - header.box[::2]
    dx = boxLength / header.gridSize
    offset -= options.center
    pos = np.empty( (header.totalGrid,3), np.float32 )
    for i in range(3):
        pos[:,i] = offset[i] + g[i].flatten() * dx[i]     # coordinates with respect to the rotation center
    del g
    
    # get the coordinates for the sampling points
    g = np.mgrid[ 0:options.noSamples, 0:options.noSamples, 0:options.noSamples ]
    noSamples = options.noSamples * options.noSamples * options.noSamples
    temp = np.empty( (noSamples,3), np.float32 )
    for i in range(3):
        temp[:,i] = g[i].flatten() * dx[i]/options.noSamples + dx[i]/2.
    data /= noSamples
    
    
    # loop over the samples
    if options.xRotation!=None: print "Rotating %.2f degrees around the x-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.xRotation,options.center[0],options.center[1],options.center[2])
    if options.yRotation!=None: print "Rotating %.2f degrees around the y-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.yRotation,options.center[0],options.center[1],options.center[2])
    if options.zRotation!=None: print "Rotating %.2f degrees around the z-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.zRotation,options.center[0],options.center[1],options.center[2])
    print "Computing the rotation by using %i samples inside each grid cell:" % noSamples
    
    for i in range(noSamples):
        print "\t computing rotation <%i> ..." % (i+1)
        pos2 = pos + temp[i,:]
    
        # do the rotations
        if options.xRotation!=None: xAxisRotation(pos2, -options.xRotation)
        if options.yRotation!=None: yAxisRotation(pos2, -options.yRotation)
        if options.zRotation!=None: zAxisRotation(pos2, -options.zRotation)
        
        # get the contributions
        grid = ((pos2 - offset)/ dx).astype(np.int32)
        selection = (grid[:,0]>=0) * (grid[:,0]<header.gridSize[0]) * (grid[:,1]>=0) * (grid[:,1]<header.gridSize[1]) * (grid[:,2]>=0) * (grid[:,2]<header.gridSize[2])
        grid = grid[selection,:]
        result[selection] += data[ grid[:,0], grid[:,1], grid[:,2] ]
    


# output the results
if headerType is allowedFiles[0]:
    density.writeDensityData( options.outputFilename, header, result )
elif headerType is allowedFiles[1]:
    NEXUS.writeMMFData( options.outputFilename, header, result )



