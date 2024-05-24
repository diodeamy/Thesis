#!/usr/bin/python

import numpy as np
import os
import sys

typeList = ['c', 'h', 'H', 'i', 'I', 'l', 'L', 'f', 'd']
typeDict = {'c':"char", 'h':"short", 'H':"unsigned short", 'i':"integer", 'I':"unsigned integer", 'l':"long int", 'L':"unsigned long int", 'f':"float", 'd':"double"}
endianess = sys.byteorder
native_code = endianess and '<' or '>'
swapped_code = endianess and '>' or '<'



def throwError(str1, showList=False):
    print "~~~~~ ERROR ~~~~~ %s" % str1
    if showList:
        print "Available format types: ",
        print typeList
    sys.exit(1)

def changeType(Type):
    if Type == 'c': Type = 'a1'
    if Type == 'h': Type = 'i2'
    if Type == 'H': Type = 'u2'
    if Type == 'i': Type = 'i4'
    if Type == 'I': Type = 'u4'
    if Type == 'l': Type = 'i8'
    if Type == 'L': Type = 'u8'
    if Type == 'f': Type = 'f4'
    if Type == 'd': Type = 'f8'
    return Type


def readElements(f,noElements,Type,printElements):
    oldType = Type
    Type = changeType( Type )
    temp = np.fromfile( f, np.dtype(Type), noElements )
    if printElements:
        print "Printing %i '%s' elements:" % (noElements,typeDict[oldType])
        print "\tactual read : ",
        print temp
        temp2 = temp.newbyteorder()
        print "\tswaped bytes: ",
        print temp2
    else: print "Skipping %i '%s' elements:" % (noElements,typeDict[oldType])


def readType(str,j):
    totalLength = len(str)
    noElements = 1
    temp = ''           #read the digits in front of the format type
    while ( j<totalLength and str[j].isdigit() ):
        temp = temp + str[j]
        j += 1
    if len(temp)>0: noElements = int( temp )
    if j>=totalLength: throwError("The format string cannot end with a number!")
    Type = str[j]       # read the format and check if valid
    if not Type in typeList: throwError("The format '%s' is not a recognized format!" % Type, True)
    j += 1
    printElements = True
    if j<totalLength:   # check if to print or not the read elements
        temp = str[j]
        if temp in ['-']:
            printElements = False
            j += 1 
    return [ j, noElements, Type, printElements ]


if len(sys.argv)<=2:
	print "The program needs at least two options:\n\
    \t1: The name of the binary file from where to read.\n\
    \t2: The type of the elements to be read ('c' char, 'h' short, 'i' integer, 'l' for long int, 'f' and 'd' float and double - capital letters for unsigned types). Can also use 'nc' to read n char characters (similar for 'h', 'i', 'l', 'f' and 'd') Use '-' after the elements that should not be printed on the screen. For example: i4f8i-2d will print an int, 4 floats, will skip over the 8 ints and print 2 more doubles.\n"
	sys.exit(1)

inputFile = sys.argv[1]
f = open(inputFile, 'rb')
print "Reading the variables from the binary file. "

for i in range(len(sys.argv)-2) :
    formatString = sys.argv[i+2]
    j = 0
    while (j<len(formatString)):
        [ j, noElements, Type, printElements ] = readType( formatString, j )
        readElements( f, noElements, Type, printElements )
        

f.close()
