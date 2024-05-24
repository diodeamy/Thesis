#!/usr/bin/python


import array
import os
import sys


fillSize = 256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8


if len(sys.argv)!=4:
    print "The program needs two options:\n\
    \t1: The name of the Gadget snapshot file that has to be converted from one endian style to another.\n\
    \t2: The name of the output file where to write the result.\n\
    \t3: The number of files\n."
    sys.exit(1)
inputFile = sys.argv[1]
outputFile = sys.argv[2]
noFiles = int(sys.argv[3])

#llop over number of files
for i in range(noFiles) :
    inFile = inputFile + "%i" % i
    print "\nReading file no %i : '%s' ..." % (i,inFile)
    f = open(inFile, 'rb')
    
    #reading header
    buffer11 = array.array('i')
    buffer11.fromfile(f, 1)
    buffer11.byteswap()
    
    noParticle = array.array('i')
    noParticle.fromfile(f,6)
    noParticle.byteswap()
    
    mass = array.array('d')
    mass.fromfile(f,8) #6 mass + time + redshift
    mass.byteswap()
    
    noTotParticles = array.array('i')
    noTotParticles.fromfile(f,10) # 2 ints + 6 totParticles + 2 ints
    noTotParticles.byteswap()
    
    boxLength = array.array('d')
    boxLength.fromfile(f,4) # 4 doubles
    boxLength.byteswap()
    
    fill = array.array('c')
    fill.fromfile(f,fillSize)
    fill.byteswap()
    
    buffer12 = array.array('i')
    buffer12.fromfile(f, 1)
    buffer12.byteswap()
    
    #print "\n\nHeader buffer 1 = %i,  buffer 2 = %i" % (buffer11[0], buffer12[0])
    print "No particles this file: %i, %i, %i, %i, %i, %i" % (noParticle[0],noParticle[1],noParticle[2],noParticle[3],noParticle[4],noParticle[5])
    #print "Scaling factor: %d" % mass[6]
    
    
    
    #reading positions
    buffer21 = array.array('i')
    buffer21.fromfile(f, 1)
    buffer21.byteswap()
    
    positions = array.array('f')
    positions.fromfile(f,3*noParticle[1])
    positions.byteswap()
    
    buffer22 = array.array('i')
    buffer22.fromfile(f, 1)
    buffer22.byteswap()
    
    #print "\nPositions particles 1 = %i,  particles 2 = %i" % (buffer21[0]/12, buffer22[0]/12)
    #print "Expected particles = %i " % (noParticle[1])
    
    
    
    #reading velocities
    buffer31 = array.array('i')
    buffer31.fromfile(f, 1)
    buffer31.byteswap()
    
    velocities = array.array('f')
    velocities.fromfile(f,3*noParticle[1])
    velocities.byteswap()
    
    buffer32 = array.array('i')
    buffer32.fromfile(f, 1)
    buffer32.byteswap()
    
    #print "\nVelocities particles 1 = %i,  particles 2 = %i" % (buffer31[0]/12, buffer32[0]/12)
    #print "Expected particles = %i " % (noParticle[1])
    
    
    
    #reading identities
    buffer41 = array.array('i')
    buffer41.fromfile(f, 1)
    buffer41.byteswap()
    
    identities = array.array('l')
    identities.fromfile(f,noParticle[1])
    identities.byteswap()
    
    buffer42 = array.array('i')
    buffer42.fromfile(f, 1)
    buffer42.byteswap()
    
    #print "\nIdentities particles 1 = %i,  particles 2 = %i" % (buffer41[0]/8, buffer42[0]/8)
    #print "Expected particles = %i " % (noParticle[1])
    
    
    
    #reading masses
    if ( mass[1]<=0. ) :
        buffer51 = array.array('i')
        buffer51.fromfile(f, 1)
        buffer51.byteswap()
        
        masses = array.array('f')
        masses.fromfile(f,3*noParticle[1])
        masses.byteswap()
        
        buffer52 = array.array('i')
        buffer52.fromfile(f, 1)
        buffer52.byteswap()
        
        #print "\nMasses particles 1 = %i,  particles 2 = %i" % (buffer51[0]/12, buffer52[0]/12)
        #print "Expected particles = %i " % (noParticle[1])
    
    f.close()
    
    
    #Writting the Endian swap results to an output file
    outFile = outputFile + "%i" % i
    print "Writting the file '%s' ..." % outFile
    f = open(outFile, 'wb')
    
    buffer11.tofile(f)
    noParticle.tofile(f)
    mass.tofile(f)
    noTotParticles.tofile(f)
    boxLength.tofile(f)
    fill.tofile(f)
    buffer12.tofile(f)
    
    buffer21.tofile(f)
    positions.tofile(f)
    buffer22.tofile(f)
    
    buffer31.tofile(f)
    velocities.tofile(f)
    buffer32.tofile(f)
    
    buffer41.tofile(f)
    identities.tofile(f)
    buffer42.tofile(f)
    
    if ( mass[1]<=0. ) :
        buffer51.tofile(f)
        masses.tofile(f)
        buffer52.tofile(f)
    
    f.close()