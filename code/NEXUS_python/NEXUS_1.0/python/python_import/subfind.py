
import numpy as np
import os
import scipy.optimize
import analysis
from miscellaneous import throwError, throwWarning


def ReadData(file, dtype, count, byteswap):
    if byteswap:
        return np.fromfile(file, dtype, count).byteswap()
    else:
        return np.fromfile(file, dtype, count)

def SkipInt32(file,byteswap):
    if byteswap:
        return (np.fromfile(file, np.int32, 1).byteswap())[0]
    else:
        return np.fromfile(file, np.int32, 1)[0]


class MergerTree:
    """Class that reads the Subfind merger trees from an input file. See for additional details:  http://astro.dur.ac.uk/~jch/password_pages/merger_trees.html"""
    byteswap = False
    noReadTrees = 0
    className = 'MergerTree'
    
    def __init__(self,inputFile=None,particleMass=1):
        """Open the input merger tree file and reads the number of merger trees and redshifts in the file."""
        self.partMass = particleMass
        if inputFile is None:
            return
        
        self.fileName = inputFile
        self.file = open(inputFile,"rb")
        self.noReadTrees = 0
        # read the buffer and determine the endian-ness
        buffer1 = np.fromfile(self.file, np.int32, 1)
        if buffer1[0]<0 or buffer1[0]>65535:
            self.byteswap=True
            buffer1 = buffer1.byteswap()
        
        #get the number of merger trees and the redshifts from the file
        self.noTrees, self.firstSnapID, self.lastSnapID = ReadData( self.file, np.int32, 3, self.byteswap )[:]
        self.noSnap = self.lastSnapID-self.firstSnapID+1
        buffer2 = ReadData( self.file, np.int32, 2, self.byteswap )
        if buffer1[0]!=buffer2[0] or buffer1[0]!=12:
            throwError( "While reading the '%s' merger tree file header. Found buffer1=%i and buffer2=%i when both should be 12." % (self.fileName,buffer1[0],buffer2[0]) )
        self.redshift = ReadData( self.file, np.float32, self.noSnap, self.byteswap )
        buffer1 = SkipInt32(self.file,self.byteswap)
        self.aFactor = 1. / (1+self.redshift)
    
    def CheckCorrectProgress(self,buffer1,buffer2,expectedValue,quantity):
        """Checks that the input file respects the expected file format."""
        if buffer1!=buffer2 or buffer1!=expectedValue:
            throwError( "In class '%s' while reading the '%s' data block for merger tree %i. The buffer preceding and following the data block do not match. Buffer1=%i while buffer2=%i, while the expected value is %i." % (self.className, quantity, self.noReadTrees, buffer1, buffer2, expectedValue) )
    
    def Next(self):
        """Reads the next merger tree from file."""
        if self.noReadTrees==self.noTrees:
            throwError( "You already read all the merger trees in the current file." )
        
        # read the number of halos in the tree
        buffer1 = SkipInt32(self.file,self.byteswap)
        self.treeFirstSnap, self.treeLastSnap, self.noHalo, self.noSubhalo = ReadData(self.file, np.int32, 4, self.byteswap)[:]
        buffer2 = SkipInt32(self.file,self.byteswap)
        self.CheckCorrectProgress(buffer1, buffer2, 4*4, "merger tree info")
        self.treeNoSnap = self.treeLastSnap-self.treeFirstSnap+1
        
        # read how many halos there are at each time step
        buffer1 = SkipInt32(self.file,self.byteswap)
        self.noHaloStep = ReadData(self.file, np.int32, self.treeNoSnap, self.byteswap)
        buffer2 = SkipInt32(self.file,self.byteswap)
        self.CheckCorrectProgress(buffer1, buffer2, 4*self.treeNoSnap, "merger tree snapshot ids")
        
        # read halo data
        buffer1 = SkipInt32(self.file,self.byteswap)
        self.haloId         = ReadData(self.file, np.int32, self.noHalo, self.byteswap)
        self.haloNoPart     = ReadData(self.file, np.int32, self.noHalo, self.byteswap)
        self.haloDesc       = ReadData(self.file, np.int32, self.noHalo, self.byteswap)
        self.haloStepDesc   = ReadData(self.file, np.int32, self.noHalo, self.byteswap)
        self.haloNoSubhalos = ReadData(self.file, np.int32, self.noHalo, self.byteswap)
        self.haloFirstSubahlo = ReadData(self.file, np.int32, self.noHalo, self.byteswap) - 1
        buffer2 = SkipInt32(self.file,self.byteswap)
        self.CheckCorrectProgress(buffer1, buffer2, 4*6*self.noHalo, "merger tree halos")
        
        # read the subhalo data
        buffer1 = SkipInt32(self.file,self.byteswap)
        self.subhaloId          = ReadData(self.file, np.int32, self.noSubhalo, self.byteswap)
        self.subhaloNoPart      = ReadData(self.file, np.int32, self.noSubhalo, self.byteswap)
        self.subhaloDesc        = ReadData(self.file, np.int32, self.noSubhalo, self.byteswap)
        self.subhaloStepDesc    = ReadData(self.file, np.int32, self.noSubhalo, self.byteswap)
        self.subhaloPos         = ReadData(self.file, np.float32, 3*self.noSubhalo, self.byteswap).reshape(-1,3)
        self.subhaloVel         = ReadData(self.file, np.float32, 3*self.noSubhalo, self.byteswap).reshape(-1,3)
        self.subhaloVelDisp     = ReadData(self.file, np.float32, self.noSubhalo, self.byteswap)
        self.subhaloVMax        = ReadData(self.file, np.float32, self.noSubhalo, self.byteswap)
        self.subhaloSpin        = ReadData(self.file, np.float32, 3*self.noSubhalo, self.byteswap).reshape(-1,3)
        self.subhaloMostBoundId = ReadData(self.file, np.int64,   self.noSubhalo, self.byteswap)
        self.subhaloRHalf       = ReadData(self.file, np.float32, self.noSubhalo, self.byteswap)
        if buffer1!=4*18*self.noSubhalo:
            self.file.seek(buffer1-4*18*self.noSubhalo, 1)
        buffer2 = SkipInt32(self.file,self.byteswap)
        self.CheckCorrectProgress(buffer1, buffer2, buffer1, "merger tree subhalos")
        
        # increase the tree count
        self.noReadTrees += 1
        self.MMPComputed = False
    
    def PresentAtFinalTime(self):
        """Checks is the last halo in the merger tree is recovered at the last snapshot (i.e. at present)."""
        if self.treeLastSnap==self.lastSnapID:
            return True
        else:
            return False
        
    def FinalMass(self):
        """Returns the mass of the last halo (most recent) in the merger tree."""
        return self.partMass * self.haloNoPart[-1]
    
    def FinalPosition(self):
        """Returns the position of the most recent halo (the top of the merger tree)."""
        tempSubhalos = np.arange(self.haloFirstSubahlo[-1], self.haloFirstSubahlo[-1]+self.haloNoSubhalos[-1])  # all the subhalos in the halo
        tempNoPart = self.subhaloNoPart[tempSubhalos].sum()                                 # gets the number of particles in all subhalos
        #print self.subhaloPos[tempSubhalos], self.subhaloNoPart[tempSubhalos].reshape(1,-1)
        tempPos = (self.subhaloPos[tempSubhalos]*self.subhaloNoPart[tempSubhalos].reshape(1,-1).transpose()).sum(0) / tempNoPart    # get the center of mass of the halo
        return tempPos
    
    def MostMassiveSubhalo(self,haloOffset):
        """Returns the offset of the most massive subhalo for the given halo ID."""
        tempSubhalos = np.arange(self.haloFirstSubahlo[haloOffset], self.haloFirstSubahlo[haloOffset]+self.haloNoSubhalos[haloOffset])
        tempMasses = self.subhaloNoPart[tempSubhalos]
        maxMass = tempMasses.max()
        return ( tempSubhalos[tempMasses==maxMass] )[0]
    
    def MostMassiveProgenitor(self):
        """Computes the most massive progenitor subhalo of the halo. It computes the subhalos ids for all available snapshots, the mass as well as the other properties for each snapshot."""
        if self.MMPComputed:
            return
        
        # find the most massive progenitor line
        indices = np.zeros( self.treeNoSnap+1, np.int32 )
        indices[1:] = self.noHaloStep.cumsum()  # computes the cummulative sum - gives the indices for the haloes at each snapshot
        tempIDs = np.zeros( self.treeNoSnap, np.int32 ) # will keep track of the position of the main progenitor at each snapshot
        tempIDs[-1] = self.noHalo-1
        tempHaloID = self.haloId[-1]            # will keep track of the descendant halo id - start with the most recent halo
        for i in range(self.treeNoSnap-2,-1,-1):
            if self.noHaloStep[i]==1:
                tempIDs[i] = indices[i]             # keep track of the halo position
                tempHaloID = self.haloId[indices[i]]# keep track of the halo id for finding the most massive parents
                continue
            # if more than one halo at a given time step, find the main progenitor
            tempIndices = np.arange(indices[i],indices[i+1])[ self.haloDesc[indices[i]:indices[i+1]]==tempHaloID ]
            if tempIndices.size==1:         #if one progenitor
                tempIDs[i] = tempIndices[0]
                tempHaloID = self.haloId[tempIndices[0]]
            elif tempIndices.size==0:       #if no progenitors
                tempIDs[0:i+1] = -1
                break
            else:                           #if multiple progenitors
                tempI   = tempIndices[ self.haloNoPart[tempIndices]==self.haloNoPart[tempIndices].max() ][0]
                tempIDs[i] = tempI
                tempHaloID = self.haloId[tempI]
        
        #get the halo properties for the main progenitor branch
        tempValid = tempIDs!=-1         # keep track of the snapshots where no valid progenitor was found
        validIDs = tempIDs[tempValid]
        self.MMP_offset = self.noSnap - tempValid.sum()  # keep track of the number of snapshots where no MMP was found (save with respect to the full number of available snapshots, not only the ones in the tree)
        offset = self.MMP_offset
        offset2 = self.noSnap - self.treeNoSnap
        
        # keep track of the halo ids for the MMP branch
        self.MMP_haloOffset = (-1 * np.ones(self.noSnap) ).astype(np.int32)
        self.MMP_haloID     = (-1 * np.ones(self.noSnap) ).astype(np.int32)
        self.MMP_haloOffset[offset:] = validIDs
        self.MMP_haloID[offset:]     = self.haloId[validIDs]
        
        # compute the number of particles in the MMP tree
        self.MMP_noPart = np.zeros( self.noSnap, np.int32 )
        self.MMP_noPart[offset:] = self.haloNoPart[validIDs]
        
        # compute the number of halo mergers
        self.MMP_noMajorMergers = np.zeros( self.noSnap, np.int32 )         # ratio 1:3 to 1:1
        self.MMP_noIntermediateMergers = np.zeros( self.noSnap, np.int32 )  # ratio 1:10 to 1:3
        self.MMP_noMinorMergers = np.zeros( self.noSnap, np.int32 )         # ratio 1:100 to 1:10
        tempHaloID = self.MMP_haloID[-1]
        for i in range(self.treeNoSnap-2,-1,-1):
            tempIndices = np.arange(indices[i],indices[i+1])[ self.haloDesc[indices[i]:indices[i+1]]==tempHaloID ]  #the progenitors of the halo at snapshot i
            if tempIndices.size==0:
                break
            elif tempIndices.size==1:
                continue
            self.MMP_noMajorMergers[offset2+i+1] = ( (self.haloNoPart[tempIndices]>=self.MMP_noPart[offset2+i]/3) * (self.haloNoPart[tempIndices]<self.MMP_noPart[offset2+i]) ).sum()
            self.MMP_noIntermediateMergers[offset2+i+1] = ( (self.haloNoPart[tempIndices]>=self.MMP_noPart[offset2+i]/10) * (self.haloNoPart[tempIndices]<self.MMP_noPart[offset2+i]/3) ).sum()
            self.MMP_noMinorMergers[offset2+i+1] = ( (self.haloNoPart[tempIndices]>=self.MMP_noPart[offset2+i]/100) *  (self.haloNoPart[tempIndices]<self.MMP_noPart[offset2+i]/10) ).sum()
        
        self.MMPComputed = True
    
    def FormationRedshift(self,minNoSnapshots=5,checkMajorMargers=False):
        """Computes the redshift of formation using different methods. Output:
            ( True/False, zHalf, zStar, zStarHalf, reducedChi )
        where:
            True/False - true if the mass accretion history is well fitted with an exponential
            zHalf - the redshift where the halo got half of its mass for the first time
            zStar - from exponential fit to the mass accretion rate (MAR) (z value where a slope of MAR is 2)
            zStarHalf - from exponential fit - when had half current mass
            reducedChi - chi divided by the number of degrees of freedom
        It uses the 'minNoSnapshots' parameter to decide betwen valid fits to the MAR history. It also checks the chi paramater to check the quality of the fit."""
        if not self.MMPComputed:
            self.MostMassiveProgenitor()
        
        # Compute the simple expression for formation time
        indexHalfMass = ( np.arange(self.noSnap)[ self.MMP_noPart>self.MMP_noPart[-1]/2 ] )[0]
        a1, a2 = self.aFactor[indexHalfMass-1:indexHalfMass+1]      # expansion before and after halo had half of the mass
        N1, N2 = self.MMP_noPart[indexHalfMass-1:indexHalfMass+1]   # mass (number particles) before and after halo had half of the mass
        slope = (N2-N1) / (a2-a1)
        const = N1 - slope*a1
        aHalf = (self.MMP_noPart[-1]/2 - const) / slope # linear interpolated expansion factor when halo got half of its current mass
        zHalf = 1./aHalf - 1.
        
        # check how many major mergers did the halo have
        valid = True
        if checkMajorMargers:
            noMergers = (self.MMP_noMajorMergers[self.MMP_noPart>self.MMP_noPart[-1]/2.]).sum()
            if noMergers>0:
                valid = False
        
        # use exponential fit to MAR
        zStar, zStarHalf, reducedChi = 0., 0., 0.
        fitFunc = lambda p, x: np.exp( p[0]*(1.-1./x) ) # the function to fit
        errFunc = lambda p, x, y: fitFunc(p,x) - y      # the error function
        p0 = [np.log(2)/zHalf]      # take redshift zHalf as initial guess
        mask   = self.MMP_noPart>0.
        aValue = self.aFactor[mask]
        MAR    = self.MMP_noPart[mask]/(self.MMP_noPart[-1]*1.)
        if MAR.size>=minNoSnapshots and valid:
            alpha, success = scipy.optimize.leastsq( errFunc, p0, args=(aValue, MAR))
            zStar, zStarHalf = 2./alpha[0]-1., np.log(2.)/alpha[0]
            reducedChi = (errFunc(alpha,aValue,MAR)**2).sum() / (MAR.size-2)
        else:
            valid, zStar, zStarHalf, reducedChi = False, -1., -1., -1.
        return (valid,zHalf,zStar,zStarHalf,reducedChi)
    
    def FormationRedshift_2(self,minNoSnapshots=5):
        """Computes the redshift of formation using different methods. Output:
            ( True/False, zHalf, zStar, zStarHalf, reducedChi )
        where:
            True/False - true if the mass accretion history is well fitted with an exponential
            zHalf - the redshift where the halo got half of its mass for the first time
            zStar - from exponential fit to the mass accretion rate (MAR) (z value where a slope of MAR is 2)
            zStarHalf - from exponential fit - when had half current mass
            reducedChi - chi divided by the number of degrees of freedom
        It uses the 'minNoSnapshots' parameter to decide betwen valid fits to the MAR history. It also checks the chi paramater to check the quality of the fit."""
        if not self.MMPComputed:
            self.MostMassiveProgenitor()
        
        # Compute the simple expression for formation time
        fractionList = [ 0.5, 0.33, 0.25, 0.1]
        res = []
        for frac in fractionList:
            indexHalfMass = ( np.arange(self.noSnap)[ self.MMP_noPart>self.MMP_noPart[-1]*frac ] )[0]
            a1, a2 = self.aFactor[indexHalfMass-1:indexHalfMass+1]      # expansion before and after halo had half of the mass
            N1, N2 = self.MMP_noPart[indexHalfMass-1:indexHalfMass+1]   # mass (number particles) before and after halo had half of the mass
            slope = (N2-N1) / (a2-a1)
            const = N1 - slope*a1
            aHalf = (self.MMP_noPart[-1]*frac - const) / slope  # linear interpolated expansion factor when halo got half of its current mass
            zHalf = 1./aHalf - 1.
            res.append(zHalf)
        
        
        # use exponential fit to MAR
        valid, zStar, zStarHalf, reducedChi = True, 0., 0., 0.
        fitFunc = lambda p, x: np.exp( p[0]*(1.-1./x) ) # the function to fit
        errFunc = lambda p, x, y: fitFunc(p,x) - y      # the error function
        p0 = [np.log(2)/res[0]]     # take redshift zHalf as initial guess
        mask   = self.MMP_noPart>0.
        aValue = self.aFactor[mask]
        MAR    = self.MMP_noPart[mask]/(self.MMP_noPart[-1]*1.)
        if MAR.size>=minNoSnapshots and valid:
            alpha, success = scipy.optimize.leastsq( errFunc, p0, args=(aValue, MAR))
            zStar, zStarHalf = 2./alpha[0]-1., np.log(2.)/alpha[0]
            reducedChi = (errFunc(alpha,aValue,MAR)**2).sum() / (MAR.size-2)
        else:
            valid, zStar, zStarHalf, reducedChi = False, -1., -1., -1.
        return (valid,res[0],res[1],res[2],res[3],zStar,zStarHalf,reducedChi)


def FileList(fileRoot,fileExtension=''):
    """Returns a list of files starting from an input root. It adds 'i' to the input root and checks if the file exists."""
    fileList = []
    noFiles = 0
    while (True):
        fileName = fileRoot + str(noFiles) + fileExtension
        if not(os.path.exists(fileName)):
            break
        noFiles += 1
        fileList.append(fileName)
    return fileList




class MMPBranchSet:
    """Class that stores the Most Massive Progenitor (MMP) branch for a set of merger trees. It has functions to write and read this data from a binary file."""
    
    def ReadMMPFile(self,fileName,VERBOSE=True):
        """Reads the binary file that stores the Most Massive Progenitor (MMP) branch for the merger trees in binary format."""
        if VERBOSE:
            print "Reading the most massive progenitor branch data from the binary file '%s' ..." % fileName
        f = open(fileName, 'rb')
        
        # Read the header information
        self.noTrees, self.noHalos, self.noSnapshots = np.fromfile( f, np.int32, 3 )
        self.redshifts = np.fromfile( f, np.float32, self.noSnapshots )
        if VERBOSE:
            print "\t found %i merger tree branches with %i halos build from %i snapshots that range from z = %.1f to %.1f ..." % (self.noTrees, self.noHalos, self.noSnapshots,self.redshifts[0],self.redshifts[-1])
        
        # Read the data position information
        self.MMP_noSnaps    = np.fromfile( f, np.int32, self.noTrees )
        self.MMP_offset     = np.fromfile( f, np.int32, self.noTrees )
        
        # Read the halo data
        self.MMP_haloID         = np.fromfile( f, np.int32, self.noHalos )
        self.MMP_haloSnaps      = np.fromfile( f, np.int32, self.noHalos )
        self.MMP_mass           = np.fromfile( f, np.float32, self.noHalos )
        self.MMP_position       = np.fromfile( f, np.float32, 3*self.noHalos ).reshape( self.noHalos, 3 )
        self.MMP_velocity       = np.fromfile( f, np.float32, 3*self.noHalos ).reshape( self.noHalos, 3 )
        self.MMP_spin           = np.fromfile( f, np.float32, 3*self.noHalos ).reshape( self.noHalos, 3 )
        self.MMP_noMergers      = np.fromfile( f, np.int32, 3*self.noHalos ).reshape( self.noHalos, 3 )
        self.MMP_velDisp        = np.fromfile( f, np.float32, self.noHalos )
        self.MMP_velMax         = np.fromfile( f, np.float32, self.noHalos )
        self.MMP_RHalf          = np.fromfile( f, np.float32, self.noHalos )
        self.MMP_mostBoundId    = np.fromfile( f, np.int64, self.noHalos )
    
    
    def WriteMMPFile(self,fileName,VERBOSE=True):
        """ Writes the Most Massive Progenitor (MMP) branch to a binary file.
        The input 'data' variable must be a class that contains the following numpy arrays:
            noTrees, noHalos, noSnaphots - number of trees, halos and snapshots
            redshifts - array giving redshift value
            MMP_noSnaps, MMP_offset - array that gives the number of snapshoots for each tree as well as the offset where each tree data starts
            MMP_haloID, MMP_haloSnaps, MMP_mass, MMP_position, MMP_velocity, MMP_spin, MMP_noMergers, MMP_velDisp, MMP_velMax, MMP_RHalf, MMP_mostBoundId - halo data for all the merger tree branches
        """
        
        if VERBOSE:
            print "Writing the results to the file '%s' (see file '%s' for informations) ..." % (outputFile,outputFile+'.info')
        
        # write a info file to know what is inside the binary file
        description = """This file specifies the format of the binary file '%s' which keeps track of the most massive progenitor (MMP) branch of the merger trees in the %s simulation. It only keeps track of the z=0 halos that have a mass of at least %.2e %s M0. The results were obtained using the command: %s
        
        The binary file has the following format:
        1) int32 - number of MMP in the file = noTrees
        2) int32 - number of Halos in the file (for all the MMP branches) = noHalos
        3) int32 - number of snapshots used to get the merger tree = noSnaphots
        4) array of float32 - size=noSnaphots - the redshift of each snapshot
        
        5) array of int32 - size=noTrees - number of snapshots in each MMP branch
        6) array of int32 - size=noTrees - offset of each MMP branch in the data to come
        
        10) array of int32   - size=noHalos   - HALO ID array
        11) array of float32 - size=noHalos   - array giving the SNAPSHOT ID in which the halo was identified
        12) array of float32 - size=noHalos   - array giving the MASS of the halo at each snapshot present in the tree
        13) array of float32 - size=3 noHalos - array giving the POSITION of each halo
        14) array of float32 - size=3 noHalos - array giving the VELOCITY of the main subhalo in the halo
        15) array of float32 - size=3 noHalos - array giving the SPIN of the main subhalo in the halo
        16) array of int32   - size=3 noHalos - array giving the NUMBER OF MERGERS of the halo at each snapshot present in the tree (1st value = major mergers, 2nd value = intermediate mergers and 3rd value = minor mergers). Major mergers = mass ration 1:3 and higher, intermediate mergers = 1:10 to 1:3 and minor mergers = 1:100 to 1:10.)
        17) array of float32 - size=noHalos   - array giving the VELOCITY DISPERSION of the main subhalo in the halo
        18) array of float32 - size=noHalos   - array giving the VELOCITY MAXIMUM of the main subhalo in the halo
        19) array of float32 - size=noHalos   - array giving the HALF RADUS of the main subhalo in the halo
        20) array of int64   - size=noHalos   - array giving the MOST BOUND PARTICLE ID of the main subhalo in the halo
        """ % (outputFile,simulationName,minMass,massUnit,programOptionsDesc)
        f = open(outputFile+'.info','w')
        f.write(description)
        f.close()
        
        f = open(outputFile,'wb')
        
        temp = np.array( [self.noTrees,self.noHalos,self.noSnaphots], np.int32)
        temp.tofile(f)
        self.redshifts.tofile(f)
        
        self.MMP_noSnaps.tofile(f)
        self.MMP_offset.tofile(f)
        
        self.MMP_haloID.tofile(f)
        self.MMP_haloSnaps.tofile(f)
        self.MMP_mass.tofile(f)
        self.MMP_position.tofile(f)
        self.MMP_velocity.tofile(f)
        self.MMP_spin.tofile(f)
        self.MMP_noMergers.tofile(f)
        self.MMP_velDisp.tofile(f)
        self.MMP_velMax.tofile(f)
        self.MMP_RHalf.tofile(f)
        self.MMP_mostBoundId.tofile(f)
        
        f.close()


class MMPBranch:
    """Selects a single Most Massive Progenitor (MMP) branch from the MMPBranchSet and computes properties for that tree branch. """
    
    def __init__(self,mmpBranchSet,treeID):
        """Copy from the MMP branch set 'mmpBranchSet' the data partaining to tree 'treeID'."""
        self.elems      = mmpBranchSet.MMP_noSnaps[treeID]
        start           = mmpBranchSet.MMP_offset[treeID]
        end             = start + self.elems
        self.start, self.end = start, end
        self.finalSnap  = mmpBranchSet.MMP_haloSnaps[end-1]
        
        self.haloID     = mmpBranchSet.MMP_haloID[start:end]
        self.redshifts  = mmpBranchSet.redshifts[ mmpBranchSet.MMP_haloSnaps[start:end] ]
        self.a          = 1. / (1. + self.redshifts)
        self.mass       = mmpBranchSet.MMP_mass[start:end]
        self.position   = mmpBranchSet.MMP_position[start:end]
        self.velocity   = mmpBranchSet.MMP_velocity[start:end]
        self.spin       = mmpBranchSet.MMP_spin[start:end]
        self.noMergers  = mmpBranchSet.MMP_noMergers[start:end]
        self.velDisp    = mmpBranchSet.MMP_velDisp[start:end]
        self.velMax     = mmpBranchSet.MMP_velMax[start:end]
        self.Rhalf      = mmpBranchSet.MMP_RHalf[start:end]
        self.mostBoundId= mmpBranchSet.MMP_mostBoundId[start:end]
    
    def Final_MMFEnvironment(self,mmfData,dx,offset):
        """Computes the final environment of the halo."""
        tempPos = ( (self.position[-1] - offset) / dx ).astype(np.int32)
        return mmfData[tempPos[0],tempPos[1],tempPos[2]]
    
    def MMFEnvironment(self,mmfData,dx,offset):
        """Returns the MMF environment at each snapshot - using the same environment for each snapshot."""
        tempPos = ( (self.position - offset) / dx ).astype(np.int32)
        return mmfData[tempPos[:,0],tempPos[:,1],tempPos[:,2]]
    
    def FormationRedshift(self,minNoSnapshots=5):
        """Computes the redshift of formation using different methods. Output:
            ( True/False, z1/2, z1/3, z1/4, z1/10, zStar, zStarHalf, r2 )
        where:
            True/False - true if the mass accretion history is well fitted with an exponential
            z1/2, z1/3,... - the redshift where the halo got 1/2, 1/3,... of its mass for the first time
            zStar - from exponential fit to the mass accretion rate (MAR) (z value where a slope of MAR is 2)
            zStarHalf - from exponential fit - when had half current mass
            r2 - square of the correlation coefficient for the linear least square fit
        It uses the 'minNoSnapshots' parameter to decide betwen valid fits to the MAR history. It also checks the chi paramater to check the quality of the fit."""
        
        # Compute the simple expression for formation time
        fractionList = [ 0.5, 0.33, 0.25, 0.1]
        res = []
        for frac in fractionList:
            index = ( np.arange(self.elems)[ self.mass>self.mass[-1]*frac ] )[0]
            a1, a2 = None, self.a[index]        # expansion before and after halo had 'fraction' of the mass
            m1, m2 = 0., self.mass[index]       # mass before and after halo had 'fraction' of the mass
            if index!=0:
                a1 = self.a[index-1]
                m1 = self.mass[index-1]
            else:
                a1 = 0.99*a2
            slope = (m2-m1) / (a2-a1)
            const = m1 - slope*a1
            aFrac = (self.mass[-1]*frac - const) / slope    # linear interpolated expansion factor when halo got half of its current mass
            zFrac = 1./aFrac - 1.
            res.append(zFrac)
        
        # use exponential fit to MAR - exp(alpha (1-1/a)) - => use linear fit ln(MAR) = - alpha 1/a + alpha
        valid, zStar, zStarHalf, r2 = True, 0., 0., 0.
        if self.mass.size>=minNoSnapshots and valid:
            mask   = self.mass>0.
            X, Y   = self.a[mask], self.mass[mask]/self.mass[-1]
            aValue = 1./X                                   # 1/a
            MAR    = np.log( Y )                            # ln( MAR )
            alpha, constant, r2 = analysis.LeastSquare_linearFit( aValue, MAR )
            alpha *= -1.        # alpha is the negative of the slope
            fitfunc = lambda p, x: np.exp(p[0]*(1.-1./x[0]))    # the fitting function
            errfunc = lambda p, x, y: fitfunc(p, x) - y         # Distance to the target function
            alpha2, success = scipy.optimize.leastsq(errfunc, (alpha,), args=(X, Y))
            #alpha2, discard, success = analysis.LeastSquare_nonlinearFit_singleParameter( X, Y, lambda x,p: np.exp(p*(1.-1./x)), alpha )
            if success: alpha = alpha2
            zStar, zStarHalf = 2./alpha-1., np.log(2.)/alpha
        else:
            valid, zStar, zStarHalf, r2 = False, -1., -1., -1.
        
        return (valid,res[0],res[1],res[2],res[3],zStar,zStarHalf,r2)
    
    def MassAccretionRate(self,noLastSnapshots=5,snapshotCenter=None):
        """Computes the mass accreation rate today using the information from the last 'noLastSnapshots' snapshots."""
        if self.elems<noLastSnapshots:
            return 0.
        if snapshotCenter==None:
            aValue = self.a[-noLastSnapshots:]
            MAR = self.mass[-noLastSnapshots:]/self.mass[-1]
            return analysis.LeastSquare_linearFit( aValue, MAR )[0]
        
        minOffset = snapshotCenter - noLastSnapshots/2 - self.finalSnap + self.elems
        maxOffset = minOffset + noLastSnapshots
        if minOffset<0 or maxOffset>=self.elems:
            return 0.
        aValue = self.a[minOffset:maxOffset]
        MAR    = self.mass[minOffset:maxOffset] / self.mass[-1]
        return analysis.LeastSquare_linearFit( aValue, MAR )[0]
        
    def TotalMergers(self):
        """Returns a 3-tuple giving the total number of major, intermediate and minor mergers."""
        return self.noMergers.sum(axis=0)
    
    
















