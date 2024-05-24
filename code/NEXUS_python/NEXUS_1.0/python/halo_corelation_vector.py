#!/usr/bin/python

import sys
import optparse
import numpy as np
import halo
import analysis
import textFile
from miscellaneous import throwError, throwWarning


# labels to give the column names in the halo data
envLabel     = halo.ENVIRONMENT_LABEL
filaDirLabel = halo.ENV_FILA_DIRECTION_LABEL
wallDirLabel = halo.ENV_WALL_DIRECTION_LABEL
hostLabel    = halo.HOST_HALO_LABEL

# some default settings
defaultMassBins = [6.88e9, 1.e10, 3.e10, 1.e11, 3.e11, 1.e12, 3.e12, 1.e13, 3.e13, 1.e14, 3.e14]

tempText = '[%.2e' % defaultMassBins[0]
for m in defaultMassBins[1:]: tempText = '%s, %.2e' % (tempText,m)
tempText = '%s%s' % (tempText,']')



def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputHaloFile  vectorColumn  outputFile  [options] <arg1> ... <argN>', description="Use this program to compute the corelation alignment of the halo direction with the direction of the environment in which the halo resides.")
    
    parser.add_option('-c','--condition', action="store", type='string', default=None, dest="condition", metavar='exp.', help="Give an expression to be used as a condition for selecting only certain halos. Example: 'dataFloat[:,10]==1'.")
    parser.add_option('--halo', action="store_true",  default=False, dest="halo", help="Computations only for halos.")
    parser.add_option('--subhalo', action="store_true",  default=False, dest="subhalo", help="Computations only for subhalos.")
    
    parser.add_option('--mB','--massBins', action="store", type='string', default=None, dest="massBins", metavar='list', help="Give the mass bins used to split the halos in different mass bins. Example: '--massBins 1.e10,1.e11,1.e12' to have 2 mass bins. DEFAULT value: %s"%tempText )
    
    parser.add_option('-d','--dependence', action="store",  default=None, dest="dependence", metavar='exp.', help="Give the variable used to invertigate its corelation with halo alignment. Example: 'dataFloat[:,10]/np.mean(dataFloat[:,10])' - a variable already expressed in terms of its mean.")
    parser.add_option('--dB','--dependenceBins', action="store", type='string', default=None, dest="dependenceBins", metavar='list', help="Give the bin values used for the dependence correlation check. Use ',' separated values.")
    parser.add_option('-q','--quantity', action="store", type='string', default=None, dest="quantity", metavar='exp.', help="Give the quantity whose average to compute when using the '--dependence' option. Allowed values: 'filaCosinus' and 'wallCosinus'.")
    
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    if len(args)<3:
        throwError( "The program needs at least 3 arguments giving the input binary halo file, the first column with the vector direction and the name of the output file. Use '-h' for the help message.\n", EXIT=False)
        parser.print_help()
        sys.exit( 1 )
    options.inputFile = args[0]
    options.vectorColumn = int(args[1])
    options.outputFile = args[2]
    
    
    # read some of the other variables
    if options.massBins:
        temp = options.massBins.strip( '([]),' ).split(',')
        if len(temp)<=1: throwError( "The '--massBins' option needs at least 2 elements in the list to specify the mimum and maximum of the mass bin. Use comma ',' to separate the value." )
        print "Found %i mass bins given by %s." % (len(temp),temp)
        options.massBins = []
        for t in temp:
            options.massBins.append( float(t) )
    else: options.massBins = defaultMassBins
    if options.dependenceBins:
        temp = options.dependenceBins.strip( '([]),' ).split(',')
        if len(temp)<=1: throwError( "The '--dependenceBins' option needs at least 2 elements in the list to specify the mimum and maximum of the bin. Use comma ',' to separate the values." )
        print "Found %i dependence bins given by %s." % (len(temp),temp)
        options.dependenceBins = []
        for t in temp:
            options.dependenceBins.append( float(t) )
    elif options.dependence:
        throwError( "The options '--dependence' needs also to be given bin boundaries using the '--dependenceBins' option." )
    if options.dependence and not options.quantity:
        throwError( "The options '--dependence' needs also to be given a quantity using the '--quantity' option." )
    
    return options







# This is the main program part
options = programOptions()
programName = sys.argv[0].rsplit('/')[-1]
options.programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])



# Read the halo data
header, dataIntegers, dataFloats = halo.readHaloData( options.inputFile )



# do some error checking - check if the halo data has hostHalo and environment columns
hostColumn, hasHost = header.GetColumnIndex( hostLabel )[1:]   # checks if halo has a host column
if not hasHost and (options.halo or options.subhalo):
    throwError( "Cannot split the halo information into halos and subhalos since the halo data has no host halo column with the name '%s'." % hostLabel )
envColumn, hasEnv = header.GetColumnIndex( envLabel )[1:]   # checks if halo has an environment column
if not hasEnv:
    throwError( "Cannot split the halos according to the environment they reside in since the halo data has no environment column with the name '%s'." % envLabel )

filaDirColumn, hasFilaDir = header.GetColumnIndex( filaDirLabel[0] )[1:]   # checks if halo has a filament environment column
if not hasFilaDir:
    throwError( "Cannot compute the alignment with environment since no environment direction was found. The halo data has no filament direction column with the name '%s'." % filaDirLabel )
wallDirColumn, hasWallDir = header.GetColumnIndex( wallDirLabel[0] )[1:]   # checks if halo has a wall environment column
if not hasWallDir:
    throwError( "Cannot compute the alignment with environment since no environment direction was found. The halo data has no wall direction column with the name '%s'." % wallDirLabel )




# select only the halos in filaments and walls
select = (dataIntegers[:,envColumn]==2) + (dataIntegers[:,envColumn]==1)
dataIntegers = dataIntegers[select,:]
dataFloats   = dataFloats  [select,:]




# select only the halos in the mass range
massRange = options.massBins[0], options.massBins[-1]
if massRange:
    mass = dataFloats[:,header.massColumn]
    minVal, maxVal = massRange
    print "Selecting the halos in the mass range [%.2e,%.2e] ..." % (minVal,maxVal)
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




# implement the condition selection
if options.condition:
    print "Computing the results of condition '%s' ..." % options.condition
    select = eval( options.condition )
    dataIntegers = dataIntegers[select,:].copy()
    dataFloats   = dataFloats  [select,:].copy()



# get the filament and wall directions for each halo and compute the angle of the halo with its environment
c1     = options.vectorColumn
dirs   = dataFloats[:,[c1,c1+1,c1+2]]
vectorName = header.columnNamesFloats[c1,:].tostring().rstrip('\x00')

c1     = filaDirColumn
envDir = dataFloats[:,[c1,c1+1,c1+2]]
filaCosinus = np.abs( analysis.scalarProduct( dirs, envDir ) )

c1     = wallDirColumn
envDir = dataFloats[:,[c1,c1+1,c1+2]]
wallCosinus = np.abs( analysis.scalarProduct( dirs, envDir ) )
#del dirs; del envDir






# functions used to output the data
def outputGeneral(outFile,data,massBins,varName,envName,returnData=False):
    massBins = np.array( massBins )
    massBinValues = np.sqrt( massBins[:-1] * massBins[1:] )
    
    # compute the statistics
    results = halo.computeStatistics( data, massBins, massBinValues, noBootstrapSamples=100 )
    
    # write a description and output the data
    desc = '''# Contains the general results of the alignment of halo vector '%s' with '%s' direction.\n# The results were obtained using:   %s\n\n# (1)Halo_mass   (2)halo_count   (3)median   (4)median_error   (5)16th_percentile   (6)84th_percentile   (7)average   (8)average_error\n  ''' % ( varName, envName, options.programOptionsDesc )
    if returnData:
        return results
    else:
        textFile.writeTextFile( outFile, desc, results )






# output general statistics if no 'dependence' variable is present
if not options.dependence:
    # first compute the results for filaments
    env = 'fila'
    select = dataIntegers[:,envColumn]==2
    data = np.column_stack( (dataFloats[select,header.massColumn],filaCosinus[select]) )
    outputGeneral( options.outputFile+'.'+env, data, options.massBins, vectorName, env )
    
    # now compute the results for filaments but using wall environment directions
    env = 'fila-wall'
    data = np.column_stack( (dataFloats[select,header.massColumn],wallCosinus[select]) )
    outputGeneral( options.outputFile+'.'+env, data, options.massBins, vectorName, env )
    
    # now compute the results for walls
    env = 'wall'
    select = -select
    data = np.column_stack( (dataFloats[select,header.massColumn],wallCosinus[select]) )
    outputGeneral( options.outputFile+'.'+env, data, options.massBins, vectorName, env )
    




# if an additional 'dependence' was given
if options.dependence:
    var = eval( options.dependence )
    data = eval( options.quantity )
    mass = dataFloats[:,header.massColumn]
    
    bins = np.array( options.dependenceBins )
    massBins = np.array( options.massBins )
    
    print "Computing the alignment angle as a function of the '%s' variable using the bins '%s'. The quatity used is '%s' ..." % ( options.dependence, options.dependenceBins, options.quantity )
    
    items = bins.size - 1
    noMassBins = massBins.size - 1
    res = np.zeros( (noMassBins,1+7*items), np.float32 )
    # get full data
    res[:,0] = np.sqrt( massBins[:-1] * massBins[1:] )
    
    # loop over all the bins and compute their value as a function of mass
    count = 1
    noEnvHalos = {0:data.shape[0]}
    for i in range(items):
        select = (var>=bins[i]) * (var<bins[i+1])
        noEnvHalos[i+1] = np.sum( select )    #number of halos in this dependence bin
        tempData = np.column_stack( (mass[select], data[select]) )
        res[:,count:count+7] = outputGeneral( None, tempData, massBins, vectorName, 0, True )[:,1:]
        count = count + 7
    
    
    # output the data
    noObjects = 'all:  %i' % (noEnvHalos[0])
    for i in range(items): noObjects = '%s    bin_%i:  %i - %.3f%%' % (noObjects,i+1,noEnvHalos[i+1],(100.*noEnvHalos[i+1])/noEnvHalos[0])
    label = '(1)massBin'
    for i in range(items): label = '%s    (%i)bin_%i_halo_count   (%i)bin_%i_median   (%i)bin_%i_median_error   (%i)bin_%i_16th_percentile   (%i)bin_%i_84th_percentile   (%i)bin_%i_average   (%i)bin_%i_average_error' % ( label, 1+i*7+0,i+1, 1+i*7+1,i+1, 1+i*7+2,i+1, 1+i*7+3,i+1, 1+i*7+4,i+1, 1+i*7+5,i+1, 1+i*7+6,i+1 )
    
    #'''# Contains the general results of the alignment of halo vector '%s' with '%s' direction.\n# The results were obtained using:   %s\n\n# (1)Halo_mass   (2)halo_count   (3)median   (4)median_error   (5)16th_percentile   (6)84th_percentile   (7)average   (8)average_error\n  ''' % ( varName, envName, options.programOptionsDesc )
    
    desc = '''# Contains the alignment results of halo vector '%s'. The valid halos were split according to the '%s' bins in the variable '%s'.\n#The number of halos in each bin: %s \n# The results were obtained using the program options: %s\n\n#%s\n''' % ( vectorName, options.dependenceBins,options.dependence, noObjects, options.programOptionsDesc, label )
    textFile.writeTextFile( options.outputFile, desc, res )

