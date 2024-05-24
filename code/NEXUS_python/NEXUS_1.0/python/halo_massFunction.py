#!/usr/bin/python

import sys
import optparse
import math
import numpy as np
import halo
import analysis
import textFile
from miscellaneous import throwError, throwWarning


# labels to give to the columns that will be added
envLabel   = halo.ENVIRONMENT_LABEL
hostLabel  = halo.HOST_HALO_LABEL


def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  inputHaloFile  outputFile [options] <arg1> ... <argN>', description="Use this program to compute the mass function for data in a binary halo file.")
    parser.add_option('--diff', action="store_true",  default=False, dest="differentialMassFunc", help="Choose this option if to compute the differential mass function. DEFAULT: computes cumulative mass function.")
    parser.add_option('--env', action="store_true",  default=False, dest="environment", help="Choose this option if to compute the mass function as a function of environment. DEFAULT: disabled.")
    parser.add_option('--halo', action="store_true",  default=False, dest="halo", help="Choose this option if to compute the mass function only for halos. DEFAULT: combined halos and subhalos.")
    parser.add_option('--subhalo', action="store_true",  default=False, dest="subhalo", help="Choose this option if to compute the mass function only for subhalos. DEFAULT: combined halos and subhalos.")
    
    parser.add_option('--range', action="store", nargs=2, default=None, type='float', dest="range", help="Specify the lower and upper range for the halo mass function [DEFAULT: 1e10 to 1e15 Msolar/h]. Needs two arguments to specify the mimimum and maximum mass range. Use '-1' to specify that the program should take the lowest or highest halo mass.")
    parser.add_option('--bins', action="store", default=None, type='int', dest="no_bins", help="Set the number of bins for the halo mass function (DEFAULT: 10 bins per mass decade).")
    parser.add_option('--hostHalo', action="store", default=hostLabel, type='string', dest="hostHaloLabel", help="Give a different column name to be used as the host halo column. DEFAULT value: %s" %hostLabel)
    
    
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    if len(args)<2:
        throwError( "The program needs at least 2 arguments giving the input binary halo file and the name of the output file. Use '-h' for the help message.\n", EXIT=False)
        parser.print_help()
        sys.exit( 1 )
    options.inputFile = args[0]
    options.outputFile = args[1]
    
    options.minMass, options.maxMass = None, None
    options.maxMass = None
    if options.range is not None:
        options.minMass, options.maxMass = options.range
        if options.minMass==-1: options.minMass = None
        if options.maxMass==-1: options.maxMass = None
    else:
        options.minMass, options.maxMass = 1.e10, 1.e15
    
    return options


def cumulativeMassFunction(massData,massRange,noBins,boxVolume):
    """Computes the cumulative mass function for the given data."""
    haloHist, histBins = analysis.Histogram( massData, no_bins=noBins, Range=massRange, bin_type="logarithm" )
    # compute the cumulative sum
    cumHist = analysis.CumulativeSum( haloHist, order="descending") #this gives N(>M)
    haloMassNumber = cumHist
    haloMassFunct = haloMassNumber / boxVolume
    return np.column_stack( (histBins,haloHist,haloMassNumber,haloMassFunct) )

def differentialMassFunction(massData,massRange,noBins,boxVolume):
    """Computes the cumulative mass function for the given data."""
    haloHist, histBins = analysis.Histogram( massData, no_bins=noBins, Range=massRange, bin_type="logarithm" )
    # divide by the bin size to get the mass function "dN / dlog(M)"
    dx = math.log10( massRange[1]/massRange[0] ) / noBins
    haloMassFunct = haloHist / dx / boxVolume
    return np.column_stack( (histBins,haloHist,haloMassFunct) )
    


# This is the main program part
options = programOptions()
programName = sys.argv[0].rsplit('/')[-1]
options.programOptionsDesc = programName + ' ' + ' '.join(sys.argv[1:])



# Read the halo data
header, dataIntegers, dataFloats = halo.readHaloData( options.inputFile )
mass = dataFloats[:,header.massColumn]
boxVolume = header.BoxLength().prod() * header.mpcUnit**3



# do some error checking - check if the halo data has hostHalo and environment columns
hostColumn, hasHost = header.GetColumnIndex( options.hostHaloLabel )[1:]   # checks if halo has a host column
if not hasHost and (options.halo or options.subhalo):
    throwError( "Cannot split the halo information into halos and subhalos since the halo data has no host halo column with the name '%s'." % options.hostHaloLabel )
envColumn, hasEnv = header.GetColumnIndex( envLabel )[1:]   # checks if halo has an environment column
if not hasEnv and options.environment:
    throwError( "Cannot split the halo information according to the environment they reside in since the halo data has no environment column with the name '%s'." % envLabel )




# select only the halos in the mass range
if options.minMass is None: options.minMass = mass.min()
if options.maxMass is None: options.maxMass = mass.max()
if options.no_bins is None: options.no_bins = math.ceil( math.log10(options.maxMass/options.minMass) * 10 )
noBins = options.no_bins
select = (mass>=options.minMass) * (mass<=options.maxMass)
dataIntegers = dataIntegers[select,:]
dataFloats   = dataFloats  [select,:]
mass         = mass[select]
    



# Select halo/subhalo populations
if options.halo or options.subhalo:
    label = ['subhalos', 'halos'] [int(options.halo)]
    print "Selecting only the %s for mass function computation ..." % label
    
    select = dataIntegers[:,hostColumn]==-1
    if options.subhalo:
        select = - select
    
    # keep only the halo/subhalos
    dataIntegers = dataIntegers[select,:].copy()
    dataFloats   = dataFloats  [select,:].copy()
    mass         = mass        [select].copy()



# compute the cumulative mass function
envs = {3:'node', 2:'filament', 1:'wall', 0:'field' }
if not options.differentialMassFunc:
    print "Computing the cumulative mass function in the mass range (%.2e,%.2e) using %i bins for %i halos in a %.2e (Mpc/h)^3 volume ..." % (options.minMass,options.maxMass,noBins,mass.shape[0],boxVolume)
    desc = '''# Gives the cumulative mass function.\n# These results were obatined using: %s\n\n# (1)halo_mass   (2)num_halos_in_mass_bin  (3)cumulative_halo_number  (4)cum_mass_funct   %s\n''' % (options.programOptionsDesc,'%s')
    
    res = cumulativeMassFunction( mass, (options.minMass,options.maxMass), noBins, boxVolume )  # mass function for all the halos
    if options.environment: 
        envCol = dataIntegers[:,envColumn]
        for i in [3,2,1,0]: # loop over the environments and output the mass function in each one
            select = envCol==i
            temp = cumulativeMassFunction( mass[select], (options.minMass,options.maxMass), noBins, boxVolume )
            res = np.column_stack( ( res, temp[:,1:] ) )
            desc = desc % ( 'columns_(2)_to_(4)_for_%s   %s' % (envs[i],'%s') )
    
    textFile.writeTextFile( options.outputFile, desc%'', res )

else:       #if differential mass function
    print "Computing the differential mass function in the mass range (%.2e,%.2e) using %i bins for %i halos in a %.2e (Mpc/h)^3 volume ..." % (options.minMass,options.maxMass,noBins,mass.shape[0],boxVolume)
    desc = '''# Gives the differential mass function.\n# These results were obatined using: %s\n\n# (1)halo_mass   (2)num_halos_in_mass_bin  (4)diff_mass_funct   %s\n''' % (options.programOptionsDesc,'%s')
    
    res = differentialMassFunction( mass, (options.minMass,options.maxMass), noBins, boxVolume )  # mass function for all the halos
    if options.environment: 
        envCol = dataIntegers[:,envColumn]
        for i in [3,2,1,0]: # loop over the environments and output the mass function in each one
            select = envCol==i
            temp = differentialMassFunction( mass[select], (options.minMass,options.maxMass), noBins, boxVolume )
            res = np.column_stack( ( res, temp[:,1:] ) )
            desc = desc % ( 'columns_(2)_to_(4)_for_%s   %s' % (envs[i],'%s') )
    
    textFile.writeTextFile( options.outputFile, desc%'', res )







