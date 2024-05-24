#!/usr/bin/python

import sys
import optparse
import math
import numpy as np
import halo
import analysis
import textFile
from miscellaneous import throwError, throwWarning


# labels to give the column names in the halo data
numPartLabel = halo.NUMBER_PARTICLES_LABEL
envLabel     = halo.ENVIRONMENT_LABEL
hostLabel    = halo.HOST_HALO_LABEL
radiusLabel  = halo.VIRIAL_RADIUS_LABEL




def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputHaloFile  outputFile [options] <arg1> ... <argN>', description="Use this program to do operations on the halo information stored in a binary halo file.")
    parser.add_option('-m','--massRange', action="store", nargs=2, default=None, type='float', dest="massRange", help="Choose only the halos in the given mass range. Needs two arguments to specify the mimimum and maximum mass range. Use '-1' to specify that the program should take the lowest or highest halo mass.")
    parser.add_option('-n','--numPart', action="store", nargs=2, default=None, type='int', dest="numberParticlesRange", help="Choose only the halos with the given number of particle range. Needs two arguments to specify the mimimum and maximum number of particles. Use '-1' to specify that the program should take the lowest or highest number of particles.")
    parser.add_option('--subBox', action="store", nargs=6, default=None, type="float", dest="subBox", help="Select only the halo in the given subBox (subBox position as ratio of boxLength). Must supply 6 arguments giving xmin, xmax, ymin, ymax,...")
    parser.add_option('--subBoxMpc', action="store", nargs=6, default=None, type="float", dest="subBoxMpc", help="Select only the halo in the given subBox (subBox position in Mpc). Must supply 6 arguments giving xmin, xmax, ymin, ymax,... in Mpc.")
    
    parser.add_option('--environment', type="choice", default=None, choices=['field', 'node', 'fila', 'wall'], dest="envSelection", metavar='env', help="Choose to output only the halos in a given environment. Use this option to give the name of the environment for which to output the halo: 'node', 'fila', 'wall' and 'field'.")
    
    parser.add_option('--halo', action="store_true",  default=False, dest="halo", help="Output only the halos.")
    parser.add_option('--subhalo', action="store_true",  default=False, dest="subhalo", help="Output only the subhalos.")
    
    parser.add_option('-o','--output', type='choice', default='1', choices=['1','2','3'], dest="output", metavar='TYPE', help="Select a custom output file format: 1=halo binary file (DEFAULT), 2=AHF halo file, 3=text file with positions, masses and radia.")
    parser.add_option('--select', action="store", type='string', nargs=2, default=None, dest="select", metavar='list', help="Select this option to output to a text file the columns that you supplied in this selection. Needs two lists as arguments, with the first being the elements from the integer column and the second the columns in the floating point data.")
    parser.add_option('--condition', action="store", type='string', default=None, dest="condition", metavar='exp.', help="Give an expression to be used as a condition for selecting only certain halos. Example: 'dataFloat[:,10]==1'.")
    
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    if len(args)<2:
        throwError( "The program needs at least 2 arguments giving the input binary halo file and the name of the output file. Use '-h' for the help message.\n", EXIT=False)
        parser.print_help()
        sys.exit( 1 )
    options.inputFile = args[0]
    options.outputFile = args[1]
    
    if options.massRange:
        if options.massRange[0]==-1: options.massRange = None, options.massRange[1]
        if options.massRange[1]==-1: options.massRange = options.massRange[0], None
    if options.numberParticlesRange:
        if options.numberParticlesRange[0]==-1: options.numberParticlesRange = None, options.numberParticlesRange[1]
        if options.numberParticlesRange[1]==-1: options.numberParticlesRange = options.numberParticlesRange[0], None
    if options.massRange and options.numberParticlesRange:
        throwError( "You can choose only one of the options '--massRange' or '--numPart' since both options determine a mass cut selection." )
    
    return options



def outputAHFFile(filename,header,dataIntegers,dataFloats):
    """Writes the halo data to a AHF text file. """
    print "Writing the output to an AHF file."
    columnNames = '#'
    for i in range( header.columnNamesIntegers.shape[0] ):
        columnNames = '%s  (%i)%s' % ( columnNames, i+1, header.columnNamesIntegers[i].tostring().rstrip('\x00') )
    for i in range( header.columnNamesFloats.shape[0] ):
        columnNames = '%s  (%i)%s' % ( columnNames, i+1+header.noColumnsIntegers, header.columnNamesFloats[i].tostring().rstrip('\x00') )
    columnNames = '%s \n' % columnNames
    data = np.empty( (dataFloats.shape[0],header.noColumns), dataFloats.dtype )
    data[:,:header.noColumnsIntegers] = dataIntegers
    data[:,header.noColumnsIntegers:] = dataFloats
    textFile.writeTextFile( filename, columnNames, data )
    

def outputTextFile(filename, options, header, dataFloats):
    """Writes some of the halo data to a text file."""
    print "Writing the output to a text file."
    descr = '''# The file contains the halo information. It was obtained using: %s\n\n# (1)Xc   (2)Yc   (3)Zc   (4)Mvir%s\n''' % (options.programOptionsDesc,'%s')
    columns = list(header.positionColumns) + [header.massColumn]
    radiusColumn, hasRadius = header.GetColumnIndex( radiusLabel )[1:]
    if not hasRadius:
        print "No radius information found in the halo file. No radius column will be written to the output file."
    else:
        columns = columns  + [radiusColumn]
        descr = descr % ('   (5)Rvir%s')
    res = dataFloats[:,columns]
    textFile.writeTextFile( filename, descr%'', res )

def outputData(filename,options,header,dataIntegers,dataFloats):
    if options.output == '1':
        halo.writeHaloData( filename, header, dataIntegers, dataFloats )
    elif options.output == '2':
        outputAHFFile( filename, header, dataIntegers, dataFloats )
    elif options.output == '3':
        outputTextFile(filename, options, header, dataFloats)
    else: throwError( "Unknow output option. Please select a valid output type - see help." )





# This is the main program part
options = programOptions()
programName = sys.argv[0].rsplit('/')[-1]
options.programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])



# Read the halo data
header, dataIntegers, dataFloats = halo.readHaloData( options.inputFile )



# do some error checking - check if the halo data has hostHalo and environment columns
hostColumn, hasHost = header.GetColumnIndex( hostLabel )[1:]   # checks if halo has a host column
if not hasHost and (options.halo or options.subhalo):
    throwError( "Cannot split the halo information into halos and subhalos since the halo data has no host halo column with the name '%s'." % options.hostHaloLabel )
envColumn, hasEnv = header.GetColumnIndex( envLabel )[1:]   # checks if halo has an environment column
if not hasEnv and options.envSelection:
    throwError( "Cannot split the halos according to the environment they reside in since the halo data has no environment column with the name '%s'." % envLabel )
numPartColumn, hasNumPart = header.GetColumnIndex( numPartLabel )[1:]   # checks if halo has a number of particles column
if not hasNumPart and options.numberParticlesRange:
    throwError( "Cannot split the halos according to the number of particles they have since there is no number of particle column with the name '%s'." % numPartLabel )





# select only the halos in the mass range
if options.massRange:
    mass = dataFloats[:,header.massColumn]
    minVal, maxVal = options.massRange
    if minVal is None: minVal = mass.min()
    if maxVal is None: maxVal = mass.max()
    print "Selecting the halos in the mass range [%.2e,%.2e] ..." % (minVal,maxVal)
    select = (mass>=minVal) * (mass<=maxVal)
    dataIntegers = dataIntegers[select,:]
    dataFloats   = dataFloats  [select,:]

if options.numberParticlesRange:
    mass = dataIntegers[:,numPartColumn]
    minVal, maxVal = options.numberParticlesRange
    if minVal is None: minVal = mass.min()
    if maxVal is None: maxVal = mass.max()
    print "Selecting the halos with the number of particles in the range [%.2e,%.2e] ..." % (minVal,maxVal)
    select = (mass>=minVal) * (mass<=maxVal)
    dataIntegers = dataIntegers[select,:]
    dataFloats   = dataFloats  [select,:]



# Select halo/subhalo populations
if options.halo or options.subhalo:
    label = ['subhalos', 'halos'] [int(options.halo)]
    print "Selecting only %s objects ..." % label
    
    select = dataIntegers[:,hostColumn]==-1
    if options.subhalo:
        select = - select
    
    # keep only the halo/subhalos
    dataIntegers = dataIntegers[select,:].copy()
    dataFloats   = dataFloats  [select,:].copy()



# select only the halos in a given subregion of the box
options.subBoxLength = None
if options.subBox:
    boxLength = header.BoxLength()
    options.subBoxLength = header.box.copy()
    options.subBoxLength[::2] = header.box[::2] +  boxLength*options.subBox[::2]
    options.subBoxLength[1::2] = header.box[::2] +  boxLength*options.subBox[1::2]
elif options.subBoxMpc:
    options.subBoxLength = options.subBoxMpc * header.mpcUnit
if options.subBox or options.subBoxMpc:
    box = options.subBoxLength
    print "\nFinding the halos in the subBox %s Mpc/h ..." % str(box/header.mpcUnit)
    pos = dataFloats[:,header.positionColumns]
    select = (pos[:,0]>=box[0]) * (pos[:,0]<=box[1]) * (pos[:,1]>=box[2]) * (pos[:,1]<=box[3]) * (pos[:,2]>=box[4]) * (pos[:,2]<=box[5])
    dataIntegers = dataIntegers[select,:].copy()
    dataFloats   = dataFloats  [select,:].copy()
    print "\tfound %i halos in the subbox." % np.sum(select)



# implement the condition selection
if options.condition:
    print "Computing the results of condition '%s' ..." % options.condition
    select = eval( options.condition )
    dataIntegers = dataIntegers[select,:].copy()
    dataFloats   = dataFloats  [select,:].copy()



# write the output data
header.AddProgramCommands( options.programOptionsDesc )
header.noHalos = dataIntegers.shape[0]
if options.select:
    res1, res2 = None, None
    has1, has2 = False, False
    columns = options.select[0].strip("[()]")
    if columns is not '':
        res1 = dataIntegers[:,eval(columns)]
        has1 = True
    columns = options.select[1].strip("[()]")
    if columns is not '':
        res2 = dataFloats[:,eval(columns)]
        has2 = True
    if not has2: res2 = res1.astype( dataFloats.dtype )
    elif has1: res2 = np.column_stack( (res1,res2) )
    textFile.writeTextFile( options.outputFile, '# Results obtained using: %s\n\n' % options.programOptionsDesc, res2 )
else:
    outputData( options.outputFile, options, header, dataIntegers, dataFloats )


