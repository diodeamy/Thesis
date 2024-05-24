#!/usr/bin/python2.6

import sys
import optparse
import numpy as np
import density
import MMF
import analysis
from AHF import AHFHalos, readAHFHalos
from textFile import readTextFile, writeTextFile, writeTextFile_gnuplot3D
from miscellaneous import throwError, throwWarning

rhoCritical = 27.7538e+10    # the critical density in units of (M0/h) / (Mpc/h)^3



def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  cleanMMF_file outputFile [options] <arg1> ... <argN>', description="Use this program the analyze the MMF node output between different MMf runs and also to compare with AHF halos.")
    parser.add_option('--compareHalos', action="store", default=None, dest="compareHalos", metavar='(density_file,AHF_file,MMF_clean_file)', help="Insert a tuple giving: '(density_file, AHF file with the corresponding halos, additional MMF_file)'. If you want to omit any of the last two, just use 'None' - you must give the density file.")
    
    # Check that the options satisfy some conditions
    options, args = parser.parse_args()
    options.args = args
    options.ahfFile = None
    options.mmfFile = None
    options.outputFile = None
    options.densityFile = None
    options.mmfFile2 = None
    
    if len(args)<2:
        throwError( "Imcomplete number of program arguments. The program needs at least 2 arguments giving: cleanMMF_file outputFileName. Use the '--help' command for additional details.", EXIT=False )
        parser.print_help()
        sys.exit(1)
    options.mmfFile = args[0]
    options.outputFile = args[1]
    
    if options.compareHalos is not None:
        temp = options.compareHalos.strip('(), ').split(',')
        if len(temp)!=3:
            throwError( "The option '--compareHalos' needs a 3 value tuple giving: '(density_file,AHF_file,2nd_MMF_file)'. You inserted '%s' which was interpreted by the program as '%s'." % (options.compareHalos,temp) )
        for i in range(len(temp)):
            if temp[i]=="None": temp[i] = None
        options.densityFile = temp[0]
        options.ahfFile = temp[1]
        options.mmfFile2 = temp[2]
        if options.densityFile is None: throwError( "The option '--compareHalos' needs a density file for input. You cannot use None." )
    
    return options


# This is the main program part
options = programOptions()
options.programOptionsDesc = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])


# open the MMF data files
mmfHeader, mmfData = None, None
boxLength = None
if options.mmfFile is not None:
    headerType = MMF.getHeaderType(options.mmfFile)
    if headerType is not 'MMF': throwError( "The first argument supplied to the '%s' program should be the name of an MMF file that gives the cosmic web environments." % sys.argv[0] )
    mmfHeader, mmfData = MMF.readMMFData(options.mmfFile)
    mmfData.shape = mmfHeader.gridSize
    boxLength = np.zeros(3,np.float32)
    boxLength[:] = mmfHeader.box[1::2] - mmfHeader.box[0::2]

mmfHeader2, mmfData2 = None, None
if options.mmfFile2 is not None:
    headerType = MMF.getHeaderType(options.mmfFile2)
    if headerType is not 'MMF': throwError( "The second MMF file '%s' is not an MMF response file." % options.mmfFile2 )
    mmfHeader2, mmfData2 = MMF.readMMFData(options.mmfFile2)
    mmfData2.shape = mmfHeader2.gridSize

# Chech that the files contain the clean response for nodes
if options.mmfFile is not None:
    if mmfHeader.fileType not in [20,21]: throwError( "The MMF response file supplied must contain an MMF clean response (file type 20) or the combined MMF clean response for all environments (file type 21)." )
    if mmfHeader.feature not in [MMF.MMF_NODE,MMF.MMF_ALL]: throwError( "The MMF response file supplied must contain the MMF response for nodes." )
    if mmfHeader.feature==MMF.MMF_ALL:
        mmfData[mmfData==MMF.MMF_FILAMENT] = 0
        mmfData[mmfData==MMF.MMF_WALL] = 0
        mmfData[mmfData==MMF.MMF_NODE] = 1
if options.mmfFile2 is not None:
    if mmfHeader2.fileType not in [20,21]: throwError( "The second MMF response file supplied must contain an MMF clean response (file type 20) or the combined MMF clean response for all environments (file type 21)." )
    if mmfHeader2.feature not in [MMF.MMF_NODE,MMF.MMF_ALL]: throwError( "The second MMF response file supplied must contain the MMF response for nodes." )
    if mmfHeader2.feature==MMF.MMF_ALL:
        mmfData2[mmfData2==MMF.MMF_FILAMENT] = 0
        mmfData2[mmfData2==MMF.MMF_WALL] = 0
        mmfData2[mmfData2==MMF.MMF_NODE] = 1


# read the density information
densityData, densityHeader = None, None
if options.densityFile is not None:
    headerType = MMF.getHeaderType(options.densityFile)
    if headerType is not 'Density': throwError( "Expected that the file '%s' to be a density file. Found that this file is a '%s' file." % (options.densityFile,headerType) )
    densityHeader, densityData = density.readDensityData(options.densityFile)
    densityData.shape = densityHeader.gridSize


# read the halo information, if any
ahfData = None
if options.ahfFile is not None:
    ahfData = readAHFHalos(options.ahfFile)
id = ".array['id']"        # 'halosData' array indices giving the id of the halos
pos = ".array['pos']"      # 'halosData' array indices giving the position of the halos
part = ".array['npart']"   # 'halosData' array indices giving the number of particles in each halo
mass = ".array['Mvir']"    # 'halos' array indices giving the mass of the halos
vel = ".array['vel']"      # 'halos' CM velocities
radius = ".array['Rvir']"  # 'halos' array indices giving the virial radius of the halos



# do the actual computations
if options.compareHalos is not None:
    
    # get the masses and centers of the objects
    mmfObjects1, mmfObjects2 = None, None
    mmfObjectsVolume1, mmfObjectsVolume2 = None,None
    mmfObjectsMass1, mmfObjectsMass2 = None,None
    mmfObjectsCenter1, mmfObjectsCenter2 = None,None
    if mmfData is not None:
        mmfObjects1, mmfObjectsVolume1 = MMF.identifyMMFObjects(mmfData)
        mmfObjectsMass1 = MMF.MMFObjectsMass(mmfObjects1,densityData,boxLength) * rhoCritical * mmfHeader.Omega0
        mmfObjectsCenter1 = MMF.MMFObjectsCenter(mmfObjects1,densityData,boxLength)
    if mmfData2 is not None:
        mmfObjects2, mmfObjectsVolume2 = MMF.identifyMMFObjects(mmfData2)
        mmfObjectsMass2 = MMF.MMFObjectsMass(mmfObjects2,densityData,boxLength) * rhoCritical * mmfHeader2.Omega0
        mmfObjectsCenter2 = MMF.MMFObjectsCenter(mmfObjects2,densityData,boxLength)
    
    if ahfData is None: throwError( "The option '--compareHalos' needs an input AHF halo file to be able to work properly." )
    ahfCenters = eval('ahfData'+pos).copy()
    ahfMass = eval('ahfData'+mass).copy()
    ahfVolume = eval('ahfData'+radius)**3 * np.float32(4*3.14/3. / 1000.**3)
    
    # find the corresponding objects in the different sets
    description = "# This file contains a match of the corresponding objects in different analysis runs.\n# \tFirst 3 columns contains the AHF halo information: halo id, halo mass and halo volume (only for halos with matched objects).\n"
    results = np.column_stack( (eval('ahfData'+id),ahfMass,ahfVolume) )
    merger1, merger2 =  np.zeros( ahfCenters.shape[0], np.int32 ), np.zeros( ahfCenters.shape[0], np.int32 )
    mergerCenterDistance1, mergerCenterDistance2 = np.zeros( ahfCenters.shape[0], np.int32 ), np.zeros( ahfCenters.shape[0], np.int32 )
    if mmfObjectsCenter1 is not None:
        merger1, mergerCenterDistance1 = MMF.matchObjectsAccordingToCenter( ahfCenters, mmfObjectsCenter1, tolerance=2.*boxLength[0]/mmfData.shape[0], boxLength=boxLength )
        print "The code matched %i (out of %i total objects) for the first analysis" % ( (merger1!=-1).sum(), mmfObjects1.max()+1 )
    if mmfObjectsCenter2 is not None:
        merger2, mergerCenterDistance2 = MMF.matchObjectsAccordingToCenter( ahfCenters, mmfObjectsCenter2, tolerance=2.*boxLength[0]/mmfData2.shape[0], boxLength=boxLength )
        print "The code matched %i (out of %i total objects) for the second analysis" % ( (merger2!=-1).sum(), mmfObjects2.max()+1 )
    validIndices = (merger1!=-1)*(merger2!=-1)
    results = results[ validIndices ].copy()
    mergerCenterDistance1 = mergerCenterDistance1[validIndices].copy()
    mergerCenterDistance2 = mergerCenterDistance2[validIndices].copy()
    match1 = merger1[ validIndices ].copy()
    match2 = merger2[ validIndices ].copy()
    print "<<<<  Common objects in both matches: %i\n" % match1.size
    
    histogramResult = None
    if mmfObjectsCenter1 is not None:
        objectId = np.arange( mmfObjects1.max()+1 ).astype(np.int32)
        temp = np.column_stack( (objectId[match1],mmfObjectsMass1[match1],mmfObjectsVolume1[match1],mergerCenterDistance1) )
        description += "# \tThe next 4 columns contain the MMF node information for the first MMF input file: node label (according size), node mass, node volume and CM distance from the center of the halo it ahs been identified with.\n"
        results = np.column_stack( (results,temp) )
        tempData, tempBins = analysis.Histogram(mergerCenterDistance1, no_bins=20, bin_type="linear" )
        histogramResult = np.column_stack( (tempBins,tempData) )
    if mmfObjectsCenter2 is not None:
        objectId = np.arange( mmfObjects2.max()+1 ).astype(np.int32)
        temp = np.column_stack( (objectId[match2],mmfObjectsMass2[match2],mmfObjectsVolume2[match2],mergerCenterDistance2) )
        description += "# \tThe next 4 columns contain the MMF node information for the second MMF input file: node label (according size), node mass, node volume and CM distance from the center of the halo it ahs been identified with.\n"
        results = np.column_stack( (results,temp) )
        tempData, tempBins = analysis.Histogram(mergerCenterDistance2, no_bins=20, bin_type="linear" )
        if histogramResult is None: histogramResult = np.column_stack( (tempBins,tempData) )
        else: histogramResult = np.column_stack( (histogramResult,tempData) )
    
    # write the results to file
    description += "\n# The results were obtained using the command: %s \n\n" % options.programOptionsDesc
    writeTextFile(options.outputFile,description,results)
    
    # write the histogram of the center distance
    description = "# This file contains the histogram of the distance between the centers of MMF nodes and the corresponding AHF halo.\n# The first column gives the histogram bins while the 2nd and/or 3rd columns gives the histogram data.\n# The results were obtained using the command: %s \n\n" % options.programOptionsDesc
    writeTextFile(options.outputFile+'.distanceHist',description,histogramResult)
    
    
    # Compute the point distribution for the mergerCenterDistance and massDifference and also VolumeDifference
    print "Computing the point distribution for center distance as a function of mass difference '(MMF mass - AHF mass) / AHF mass' ... "
    data1 = np.column_stack( ( np.abs(results[:,1]-results[:,4])/results[:,1], results[:,6]/(results[:,2]/(4*3.14/3))**.33  ) ).copy()  # first column stores the mass difference while the second stores the distance difference in virial radius fractions
    data1Distribution = analysis.PointDistribution(data1,noBins=(22,22),xRange=(-0.1,1),yRange=(-0.1,1))
    data1Density = data1Distribution[:,:,2].copy()
    markers = analysis.FindContourValue(data1Density,enclosed_fraction=[0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99])
    
    description = "# This file gives the point distribution (count in cell) of the DeltaM/M to distance/R_vir for the match of MMF nodes and AHF halos.\n# The first columns gives the DeltaM/M bin value, the second the distance/R_vir bin value while the 3rd gives the number of points in that cell.\n# The following are contour line values enclosing a given fraction of the total number of points - [desired fraction, count threshold, obtained fraction] -: %s \n# The results were obtained using the command: %s \n\n" % (str(markers).replace('\n',','),options.programOptionsDesc)
    writeTextFile_gnuplot3D(options.outputFile+'.massVariationDistribution',description,data1Distribution)
    
    print "Computing the point distribution for center distance as a function of volume difference '(MMF volume - AHF volume) / AHF volume' ... "
    data1 = np.column_stack( ( np.abs(results[:,2]-results[:,5])/results[:,2], results[:,6]/(results[:,2]/(4*3.14/3))**.33  ) ).copy()  # first column stores the mass difference while the second stores the distance difference in virial radius fractions
    data1Distribution = analysis.PointDistribution(data1,noBins=(22,22),xRange=(-0.1,1),yRange=(-0.1,1))
    data1Density = data1Distribution[:,:,2].copy()
    markers = analysis.FindContourValue(data1Density,enclosed_fraction=[0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99])
    
    description = "# This file gives the point distribution (count in cell) of the DeltaV/V to distance/R_vir for the match of MMF nodes and AHF halos.\n# The first columns gives the DeltaV/V bin value, the second the distance/R_vir bin value while the 3rd gives the number of points in that cell.\n# The following are contour line values enclosing a given fraction of the total number of points - [desired fraction, count threshold, obtained fraction] -: %s \n# The results were obtained using the command: %s \n\n" % (str(markers).replace('\n',','),options.programOptionsDesc)
    writeTextFile_gnuplot3D(options.outputFile+'.volumeVariationDistribution',description,data1Distribution)
