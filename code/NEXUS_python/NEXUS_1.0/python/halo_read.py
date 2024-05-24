#!/usr/bin/python

import sys
import numpy as np
import halo
from miscellaneous import throwError, throwWarning


programName = sys.argv[0].rsplit('/')[-1]
help = '''Use this program to read the data in a halo file.
USAGE:   %s   input_halo_file  
''' % programName


if len(sys.argv) not in [2,3,4]:
    print help
    sys.exit(1)

inFile = sys.argv[1]
column = None
position = ':10'
if len(sys.argv)>=3:
    column = int(sys.argv[2])
if len(sys.argv)>=4:
    position = sys.argv[3]



header, dataIntegers, dataFloats = halo.readHaloData( inFile )

if column is not None:
    print "The %s elements for column %i: " % (position,column)
    print eval( 'dataFloats['+position+',column]' )

#columns = [4, 5, 6, 15, 20, 21, 38, 39, 40]
#texts   = ['Vcx', 'Vcy', 'Vcz', 'lambda ', 'b_shape', 'c_shape', 'env_dir_x', 'env_dir_y', 'env_dir_z']
#for i in range( len(columns)-3 ):
    #length = len( texts[i] )
    #header.columnNamesFloats[columns[i],:length] = texts[i]


#header.PrintValues()
#halo.writeHaloData( inFile, header, dataIntegers, dataFloats )

#print 'no substructures:', dataIntegers[:13,2]
#print 'Mass:', dataFloats[:10,header.massColumn]
#print 'Rvir:', dataFloats[:10,7]


#print np.sum( dataIntegers[:,1]!=-1 ), np.sum( dataIntegers[:,2] )
#print np.sum( dataIntegers[:,5]!=-1 ), np.sum( dataIntegers[:,6] )
#print dataFloats[:10,:]

