#!/usr/bin/python

import numpy as np
import optparse
import sys
from gadget import GadgetHeader, readGadgetData, writeGadgetData, selectParticlesInBox
from miscellaneous import throwError
import MMF
import textFile


MPC_UNIT = 1000.


def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputFile outputFile [options] <arg1> ... <argN>', description='Use this file to output the gadget particle data after rotating the particles along a given direction/directions. First the z-rotation is done, than the y-rotation and at the end the x-rotation.')
    parser.add_option('--subBox', action="store", nargs=6, default=None, type="float", dest="subBox", help="Select only the particles in the given subBox (subBox position as ratio of boxLength). Must supply 6 arguments giving xmin, xmax, ymin, ymax,...")
    parser.add_option('--subBoxMpc', action="store", nargs=6, default=None, type="float", dest="subBoxMpc", help="Select only the particles in the given subBox (subBox position in Mpc). Must supply 6 arguments giving xmin, xmax, ymin, ymax,... in Mpc.")
    parser.add_option('-x','--xRotation', action="store", type='float', default=None, dest="xRotation", metavar='angle', help="Give the angle in degrees for a rotation along the x-axis.")
    parser.add_option('-y','--yRotation', action="store", type='float', default=None, dest="yRotation", metavar='angle', help="Give the angle in degrees for a rotation along the y-axis.")
    parser.add_option('-z','--zRotation', action="store", type='float', default=None, dest="zRotation", metavar='angle', help="Give the angle in degrees for a rotation along the z-axis.")
    parser.add_option('-v','--vector', action="store", nargs=3, type='float', default=None, dest="vector", metavar='value', help="Give a vector as its 3 components: vx, vy and vz. The program will compute the rotation angles and rotate the particles such that this vector, after rotation, will lie along the z-axis.")
    parser.add_option('-c','--center', action="store", nargs=3, type='float', default=None, dest="center", metavar='value', help="Give the center of rotation as its 3 components: x, y and z. The program will rotate the positions around this center using the rotation coordinates specified.")
    
    
    options, args = parser.parse_args()
    if options.subBox and options.subBoxMpc:
        throwError( "You can specify only option '--subBox' or '--subBoxMpc' at a time, but not both at the same time." )
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
def xAxisRotation(particles, angle):
    """ Roates the positions around the x-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp               = cosA * particles.pos[:,1] - sinA * particles.pos[:,2]
    particles.pos[:,2] = sinA * particles.pos[:,1] + cosA * particles.pos[:,2]
    particles.pos[:,1] = temp
    # rotate the velocities
    temp               = cosA * particles.vel[:,1] - sinA * particles.vel[:,2]
    particles.vel[:,2] = sinA * particles.vel[:,1] + cosA * particles.vel[:,2]
    particles.vel[:,1] = temp

def yAxisRotation(particles, angle):
    """ Roates the positions around the y-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp               =  cosA * particles.pos[:,0] + sinA * particles.pos[:,2]
    particles.pos[:,2] = -sinA * particles.pos[:,0] + cosA * particles.pos[:,2]
    particles.pos[:,0] = temp
    # rotate the velocities
    temp               =  cosA * particles.vel[:,0] + sinA * particles.vel[:,2]
    particles.vel[:,2] = -sinA * particles.vel[:,0] + cosA * particles.vel[:,2]
    particles.vel[:,0] = temp
    
def zAxisRotation(particles, angle):
    """ Roates the positions around the y-axis using the given rotation angle. """
    a = angle * conversion
    cosA, sinA = np.cos( a ), np.sin( a )
    # rotate the positions
    temp               = cosA * particles.pos[:,0] - sinA * particles.pos[:,1]
    particles.pos[:,1] = sinA * particles.pos[:,0] + cosA * particles.pos[:,1]
    particles.pos[:,0] = temp
    # rotate the velocities
    temp               = cosA * particles.vel[:,0] - sinA * particles.vel[:,1]
    particles.vel[:,1] = sinA * particles.vel[:,0] + cosA * particles.vel[:,1]
    particles.vel[:,0] = temp


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



# read the gadget data
particles = readGadgetData( options.inputFilename, ID=False )
particles.AddIds( np.arange(particles.noParticles).astype(np.int32) )


if options.subBox:
    boxLength = particles.Header().BoxSize
    options.subBoxLength = np.array( options.subBox ) * boxLength
elif options.subBoxMpc:
    options.subBoxLength = np.array( options.subBoxMpc ) * MPC_UNIT




# if the user gave the '--vector' option, tahn compute the rotation angles along the x and y directions
if options.vector:
    angleX, angleY = computeRotationAngle( options.vector )
    options.xRotation = angleX
    options.yRotation = angleY
    options.zRotation = None


# do the requested rotations
doRotation = options.xRotation!=None or options.yRotation!=None or options.zRotation!=None
if doRotation:  #take out the center coordinate
    center = np.array( options.center[:] ) * MPC_UNIT
    particles.pos[:,:] -= center

if options.zRotation!=None:
    print "Rotating %.2f degrees around the z-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.zRotation,options.center[0],options.center[1],options.center[2])
    zAxisRotation(particles, options.zRotation)

if options.yRotation!=None:
    print "Rotating %.2f degrees around the y-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.yRotation,options.center[0],options.center[1],options.center[2])
    yAxisRotation(particles, options.yRotation)
    
if options.xRotation!=None:
    print "Rotating %.2f degrees around the x-axis through the rotation center (%.2f, %.2f, %.2f) ..." % (options.xRotation,options.center[0],options.center[1],options.center[2])
    xAxisRotation(particles, options.xRotation)

if doRotation:  #add back the center coordinate
    center = np.array( options.center[:] ) * MPC_UNIT
    particles.pos[:,:] += center



# if select particles in subBox, do so
if options.subBox or options.subBoxMpc:
    box = options.subBoxLength
    periodicLength = particles.header.BoxSize
    particles = selectParticlesInBox( particles, box, periodicLength)




#write the data to a custom or a gadget file
particles.Header().npart[1] = particles.noParticles
particles.Header().npartTotal[1] = particles.noParticles
writeGadgetData(options.outputFilename,particles)
