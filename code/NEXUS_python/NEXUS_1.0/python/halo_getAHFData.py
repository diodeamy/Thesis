#!/usr/bin/python

import sys
import numpy as np
import AHF
import halo
from miscellaneous import throwError, throwWarning


box = [ 0., 100., 0., 100., 0., 100. ]
mpcUnit = 1000.
positionColumns = [ 1, 2, 3 ]
massUnit = 1.
massColumn = 0
radiusColumn = 7
hostHaloColumn = 1 #None

# options for AHFhalos files
numberZeros = 4
maxFileCount = 2#608#None
hasID = True  #None


def readAHFHalosFile(file,VERBOSE=True):
    # find the root name for multiple text files
    largeNumber = 10**numberZeros
    zeroStr = str(largeNumber)[1:] #now this is '0..0' numberZeros
    zeroStr += '.'
    splitFilename = file.split(zeroStr)
    if len(splitFilename)==1:
        print "Single AHF file is: '%s'" % splitFilename[0]
    elif len(splitFilename)!=2:
        throwError( "Split the AHF '_halos' filename using the '%s' string in %i parts. There was an error somewhere since the program expected only 2 parts." % (zeroStr,len(splitFilename)) )
    
    # count how many lines there are in the AHF files, this tells the program the number of halos
    noValidFiles, noTotalHalos, noHalosInFile = AHF.numberOfAHFHalos( splitFilename, largeNumber, maxFileCount )
    if VERBOSE: print "Found %i valid AHF '_halos' file(s) with %i halos." % (noValidFiles,noTotalHalos)
    
    # read the halos properties in each file
    noHalos1, noHalos2 = 0, 0
    noColumns = AHF.noColumnsAHFTextFile
    columnsDescription = None
    data = None
    for i in range(noValidFiles):
        if noHalosInFile[i]<=0: continue
        noHalos2 += noHalosInFile[i]
        
        tempName = AHF.AHFMultipleFiles( splitFilename, largeNumber+i, i )
        if VERBOSE: print "Reading file %i of %i files with name '%s' - it has %i halos ... " % ( i+1,noValidFiles, tempName, noHalosInFile[i] ),
        
        if i is 0:
            f = open( tempName, 'r' )
            columnsDescription = f.readline()    #skip first line which contains the names of the columns (only for the first file)
            noColumns = len( columnsDescription.split() )
            print "\n\t found that the data has %i columns ..." % noColumns,
            data = np.empty( (noTotalHalos,noColumns), np.float32 )
            data[noHalos1:noHalos2,:] = np.fromfile( f, np.float32, noColumns*noHalosInFile[i], sep=" " ).reshape(-1,noColumns)
            f.close()
            if hostHaloColumn:
                select = data[noHalos1:noHalos2,hostHaloColumn]>=noHalos2       # set a different host halo index for halos who's hosts are in a different file
                data[noHalos1:noHalos2,hostHaloColumn][select] = noTotalHalos
        else:
            data[noHalos1:noHalos2,:] = np.fromfile( tempName, np.float32, noColumns*noHalosInFile[i], sep=" " ).reshape(-1,noColumns)
            if hasID:   # relabel the last read halos
                data[noHalos1:noHalos2,0] += noHalos1
                if hostHaloColumn:
                    select = data[noHalos1:noHalos2,hostHaloColumn]!=-1
                    data[noHalos1:noHalos2,hostHaloColumn][select] += noHalos1
                    select = data[noHalos1:noHalos2,hostHaloColumn]>=noHalos2       # set a different host halo index for halos who's hosts are in a different file
                    data[noHalos1:noHalos2,hostHaloColumn][select] = noTotalHalos
        noHalos1 = noHalos2
        if VERBOSE: print "Done."
    
    data.shape = -1, noColumns
    columnsDescription = columnsDescription[1:]
    if not hasID:
        data = np.column_stack( (np.arange(data.shape[0]),data) )
        columnsDescription = "ID(1)  %s" % columnsDescription
    
    #if True:
        #fileID = np.empty( data.shape[0] )
        #count = 0
        #for i in range(noValidFiles):
            #fileID[count:count+noHalosInFile[i]] = i
            #count += noHalosInFile[i]
        #columnsDescription = "%s fileID(0)" % columnsDescription
        #data = np.column_stack( (data,fileID) )
    
    return data, columnsDescription


def getColumnNames(description):
    s1 = description.split()
    result = np.zeros( (len(s1),halo.HaloHeader.columnNameLength), dtype='c' )
    for i in range(len(s1)):
        temp = s1[i].split( '(' )[0]
        choose = int( len(temp) > halo.HaloHeader.columnNameLength )
        l = [ len(temp), halo.HaloHeader.columnNameLength ] [ choose ]
        result[i,:l] = temp[:l]
    return result



AHFDescription = { 
                   41:'''ID(0)  npart(1)       nvpart(2)       Xc(3)   Yc(4)   Zc(5)   VXc(6)  VYc(7)  VZc(8)  Mvir(9) Rvir(10)        Vmax(11)        Rmax(12)        sigV(13)        lambda(14)      Lx(15)  Ly(16)  Lz(17)  a(18)   Eax(19) Eay(20) Eaz(21) b(22)   Ebx(23) Eby(24) Ebz(25) c(26)   Ecx(27) Ecy(28) Ecz(29) ovdens(30)      Redge(31)       nbins(32)       Ekin(33)        Epot(34)        mbp_offset(35)  com_offset(36)  r2(37)  lambdaE(38)     v_esc(39)       Phi0(40)''', 
                   43:'''ID(1)  hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)    lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)    SurfP(41)       Phi0(42)        cNFW(43)''' 
                }

IntegerColumns = { 
                    41:[ 0, 1, 31, 32 ],
                    43:[ 0, 1, 2, 4, 36 ],
                    44:[ 0, 1, 2, 4, 36, 43 ]
                }
FloatColumns = { 41:range(41), 43:range(43), 44:range(44) }
for i in IntegerColumns.keys():
    for e in IntegerColumns[i]:
        FloatColumns[i].remove(e)





help = '''Use this program to read AHF halo data and rewrite it into the halo binary format.
USAGE:  %s   input_AHF_file   output_halo_file
''' % sys.argv[0]

if len(sys.argv)!=3:
    print help
    sys.exit(1)

inFile = sys.argv[1]
outFile = sys.argv[2]
programOptionsDesc = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])



# read the data
data, description = None, None
if '.AHF_halos' in inFile:
    data, description = readAHFHalosFile(inFile)
else:
    data = AHF.readAHFHalos_binary(inFile).data
    description = AHFDescription[ data.shape[1] ]



# prepare the header of the output file
header = halo.HaloHeader()
header.noHalos = data.shape[0]
header.noColumnsIntegers = len( IntegerColumns[data.shape[1]] )
header.noColumnsFloats = len( FloatColumns[data.shape[1]] )
header.noColumns = data.shape[1]

header.mpcUnit = 1.
header.box[:] = box[:]
header.positionColumns[:] = positionColumns[:]

header.massUnit = massUnit
header.massColumn = massColumn

names = getColumnNames(description)
header.columnNamesIntegers = names[ IntegerColumns[data.shape[1]] ]
header.columnNamesFloats   = names[ FloatColumns[data.shape[1]] ]


# prepare the data for writing
dataIntegers = data[ :, IntegerColumns[data.shape[1]] ].astype(np.int32).copy()
dataFloats   = data[ :, FloatColumns[data.shape[1]] ].astype(np.float32).copy()

posCol = header.positionColumns
dataFloats[:,posCol] /= mpcUnit

massCol = header.massColumn
#dataIntegers[:,hostHaloColumn] = halo.matchSubstructures( haloPos=dataFloats[:,posCol], hostID=dataIntegers[:,hostHaloColumn], haloRadia=dataFloats[:,radiusColumn]/1000.*header.mpcUnit, haloMass=dataFloats[:,massCol] )
#print np.min(dataIntegers[:,hostHaloColumn]), np.max(dataIntegers[:,hostHaloColumn])
dataIntegers, dataFloats = halo.sortHalos( dataIntegers, dataFloats, dataFloats[:,massCol], order='descending', relabelHalos=True, parentHaloColumn=None, VERBOSE=True )


# write the data
header.AddProgramCommands(programOptionsDesc)
halo.writeHaloData( outFile, header, dataIntegers, dataFloats )
    


