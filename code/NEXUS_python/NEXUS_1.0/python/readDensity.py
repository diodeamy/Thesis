#!/usr/bin/python

from density import DensityHeader, readDensityHeader, readDensityData, writeDensityData
import sys

file = sys.argv[1]
[header, denData ] = readDensityData(file)
header.PrintValues()
print denData

#writeDensityData("test.dat",header,denData)