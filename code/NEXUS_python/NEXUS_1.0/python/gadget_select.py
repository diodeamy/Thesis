#!/usr/bin/python

import numpy as np
import optparse
import sys
from gadget import GadgetHeader, readGadgetData, writeGadgetData
from miscellaneous import throwError
import MMF
import textFile


MPC_UNIT = 1000.


def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputFile outputFile [options] <arg1> ... <argN>')
    parser.add_option('-t','--type', type='choice', default='1', choices=['1','2'], dest="snapshotType", metavar='TYPE', help="Choose the file type for the output Gadget snapshot (available values: 1 or 2).")
    parser.add_option('--subBox', action="store", nargs=6, default=None, type="float", dest="subBox", help="Select only the particles in the given subBox (subBox position as ratio of boxLength). Must supply 6 arguments giving xmin, xmax, ymin, ymax,...")
    parser.add_option('--subBoxMpc', action="store", nargs=6, default=None, type="float", dest="subBoxMpc", help="Select only the particles in the given subBox (subBox position in Mpc). Must supply 6 arguments giving xmin, xmax, ymin, ymax,... in Mpc.")
    parser.add_option('-d','--downsample', action="store", default=None, type="float", dest="downsample", metavar='number', help="Randomly choose only this fraction of particles for output.")
    parser.add_option('--environment', action="store", default=None, type="string", dest="envFile", metavar='file', help="Give the name of the input clean response file containing the Cosmic Web environments. It uses this file to locate the particle in its given environment and to output the environment tag of each particle.")
    parser.add_option('-o','--output', type='choice', default='1', choices=['1','2','3'], dest="output", metavar='TYPE', help="Select a custom output file format: 1=gadget file, 2=text file with particle positions, 3=binary file with particle positions.")
    
    
    options, args = parser.parse_args()
    if options.subBox and options.subBoxMpc:
        throwError( "You can specify only option '--subBox' or '--subBoxMpc' at a time, but not both at the same time." )
    if (not options.downsample is None) and (options.downsample>1. or options.downsample<0.):
        throwError( "The argument suplied to the '--downsample' option must be in the range [0.,1.]" )
    if len(args)<2:
        throwError( "You must supply at least two arguments to the program giving the input and output files." )
    
    options.inputFilename = args[0]
    options.outputFilename = args[1]
    
    return options
    


# function used to write a different output than a gadget file
def writeOutputFile(options,particles):
    if options.output == '2':  # write the results to a text file
        results = np.column_stack( (particles.pos[:,0],particles.pos[:,1],particles.pos[:,2]) )
        results /= MPC_UNIT
        if options.envFile:
            results = np.column_stack( (results,particles.environment) )
        textFile.writeTextFile(options.outputFilename,"",results)
    
    elif options.output == '3':    #write the positions to a binary file
        results = None
        if options.envFile:
            results = np.empty( (particles.noParticles,4), np.float32 )
            results[:,0:3] = particles.pos / MPC_UNIT
            results[:,3] = particles.environment
        else:
            results = np.empty( (particles.noParticles,3), np.float32 )
            results[:,0:3] = particles.pos / MPC_UNIT
        results.tofile( options.outputFilename )





# read program options
options = programOptions()

# read the gadget data
particles = readGadgetData( options.inputFilename )

if options.subBox:
    boxLength = particles.Header().BoxSize
    options.subBoxLength = np.array( options.subBox ) * boxLength
elif options.subBoxMpc:
    options.subBoxLength = np.array( options.subBoxMpc ) * MPC_UNIT



# if select the nevironment option, find to which environment each particle corresponds to
if options.envFile:
    print "\nFinding the environment of each particle using the Cosmic Web environments given by the file %s ..." % options.envFile
    # read the environment file
    header, data = MMF.readMMFData( options.envFile )
    if header.feature!=MMF.MMF_ALL:   # test that this is a wall file
        throwError( "File '%s' should contain the MMF clean response for all environments. But it contain an MMF response for the '%s' feature." % (allClean_1,header1.Feature()) )
    if header.fileType!=MMF.MMF_CLEAN_RESPONSE_COMBINED: # this test if file is a clean response
        throwError( "File '%s' should contain the MMF combined clean response. But it does not contain an MMF clean response, but it contains '%s'." % (allClean_1,header1.FileType()) )
    box = header.box * MPC_UNIT
    dx = (box[1::2] - box[::2]) / header.gridSize[:]
    pos = particles.pos
    gridIndex = np.empty( (particles.noParticles,3), np.int32 )
    gridIndex[:,0] = (pos[:,0]-box[0])/dx[0]
    gridIndex[:,1] = (pos[:,1]-box[2])/dx[1]
    gridIndex[:,2] = (pos[:,2]-box[4])/dx[2]
    select = (gridIndex[:,0]>=0) * (gridIndex[:,0]<header.gridSize[0]) * (gridIndex[:,1]>=0) * (gridIndex[:,1]<header.gridSize[1]) * (gridIndex[:,2]>=0) * (gridIndex[:,2]<header.gridSize[2])
    NYZ, NZ = header.gridSize[1]*header.gridSize[2], header.gridSize[2]
    gridI = gridIndex[select,0]*NYZ + gridIndex[select,1]*NZ + gridIndex[select,2]
    particles.environment = np.zeros( particles.noParticles, np.int32 )
    particles.environment[select] = data[gridI]


# if select particles in subBox, do so
if options.subBox or options.subBoxMpc:
    box = options.subBoxLength
    print "\nFinding the particles in the subBox %s ..." % str(box)
    pos = particles.pos
    particleInside = (pos[:,0]>=box[0]) * (pos[:,0]<=box[1]) * (pos[:,1]>=box[2]) * (pos[:,1]<=box[3]) * (pos[:,2]>=box[4]) * (pos[:,2]<=box[5])
    particles.Update(particleInside)
    if options.envFile: particles.environment = particles.environment[particleInside]
    print "\tfound %i particles in the subbox." % particles.noParticles


# if select 'n' random particles from the box, than do so
if options.downsample is not None:
    fraction, initialParticles = options.downsample, particles.noParticles
    print "\nDownsampling the number of particle via the factor '%f' ... " % fraction
    select = np.random.random_sample( particles.noParticles ) <= options.downsample
    particles.Update(select)
    if options.envFile: particles.environment = particles.environment[select]
    print "\tselect %i particles during the downsampling -- downsample factor = %f." % (particles.noParticles, (1.*particles.noParticles)/initialParticles )


#write the data to a custom or a gadget file
if options.output in ['2','3']:
    writeOutputFile(options,particles)
else:
    particles.Header().npart[1] = particles.noParticles
    particles.Header().npartTotal[1] = particles.noParticles
    writeGadgetData(options.outputFilename,particles)
