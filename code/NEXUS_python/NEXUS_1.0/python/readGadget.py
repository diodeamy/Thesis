#!/usr/bin/python

from gadget import GadgetParticles, GadgetHeader, readGadgetHeader, readGadgetData, writeGadgetData
import sys

file = sys.argv[1]
particles = readGadgetData(file)
particles.Header().PrintValues()
print particles.Pos()
print particles.Vel()
print particles.Ids()
print particles.Mass()

#writeGadgetData("test.dat",particles)