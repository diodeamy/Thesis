
import numpy as np
import os.path
from scipy.weave import inline, converters
from miscellaneous import throwError, throwWarning, dtype2ctype
from gadget import GadgetHeader
import analysis
try:
    import sqlite3 as sqlite
    sqlite3ImportError = None
except:
    sqlite3ImportError = "Could not find the 'sqlite3' module. Cannot read/write to a database."





noColumnsAHFTextFile = 40
noDataColumns = 41

class AHFHalos:
    """ Class used to store the properties of AHF haloes. It contains an array for each AHF halo property present in the AHF halo finder output. """
    
    def __init__(self,data=None,header=None,numberColumns=None):
        self.noHalos = np.uint64( 0 )
        self.hasData = False
        self.hasFileHeader = False
        self.noColumns = noDataColumns
        if data!=None: self.noColumns = data.shape[1]
        
        self.name = { 0:'id', 1:'npart', 2:'nvpart', 3:'Xc', 4:'Yc', 5:'Zc', 6:'VXc', 7:'VYc', 8:'VZc', 9:'Mvir', 10:'Rvir', 11:'Vmax', 12:'Rmax', 13:'sigV', 14:'lambda', 15:'Lx', 16:'Ly', 17:'Lz', 18:'a', 19:'Eax', 20:'Eay', 21:'Eaz', 22:'b', 23:'Ebx', 24:'Eby', 25:'Ebz', 26:'c', 27:'Ecx', 28:'Ecy', 29:'Ecz', 30:'ovdens', 31:'Redge', 32:'nbins', 33:'Ekin', 34:'Epot', 35:'mbp_offset', 36:'com_offset', 37:'r2', 38:'lambdaE', 39:'v_esc', 40:'Phi0' }
        self.array = { 'id':None, 'npart':None, 'nvpart':None, 'Xc':None, 'Yc':None, 'Zc':None, 'VXc':None, 'VYc':None, 'VZc':None, 'Mvir':None, 'Rvir':None, 'Vmax':None, 'Rmax':None, 'sigV':None, 'lambda':None, 'Lx':None, 'Ly':None, 'Lz':None, 'a':None, 'Eax':None, 'Eay':None, 'Eaz':None, 'b':None, 'Ebx':None, 'Eby':None, 'Ebz':None, 'c':None, 'Ecx':None, 'Ecy':None, 'Ecz':None, 'ovdens':None, 'Redge':None, 'nbins':None, 'Ekin':None, 'Epot':None, 'mbp_offset':None, 'com_offset':None, 'r2':None, 'lambdaE':None, 'v_esc':None, 'Phi0':None }
        
        if numberColumns==43:
            self.noColumns = numberColumns
            self.name = { 0:'id', 1:'hostHalo', 2:'numSubStruct', 3:'Mvir', 4:'npart', 5:'Xc', 6:'Yc', 7:'Zc', 8:'VXc', 9:'VYc', 10:'VZc', 11:'Rvir', 12:'Rmax', 13:'r2', 14:'mbp_offset', 15:'com_offset', 16:'Vmax', 17:'v_esc', 18:'sigV', 19:'lambda', 20:'lambdaE', 21:'Lx', 22:'Ly', 23:'Lz', 24:'b', 25:'c', 26:'Eax', 27:'Eay', 28:'Eaz', 29:'Ebx', 30:'Eby', 31:'Ebz', 32:'Ecx', 33:'Ecy', 34:'Ecz', 35:'ovdens', 36:'nbins', 37:'fMhires', 38:'Ekin', 39:'Epot',  40:'SurfP', 41:'Phi0', 42:'cNFW' }
            self.array = { 'id':None, 'hostHalo':None, 'numSubStruct':None, 'Mvir':None, 'npart':None, 'Xc':None, 'Yc':None, 'Zc':None, 'VXc':None, 'VYc':None, 'VZc':None, 'Rvir':None, 'Rmax':None, 'r2':None, 'mbp_offset':None, 'com_offset':None, 'Vmax':None, 'v_esc':None, 'sigV':None, 'lambda':None, 'lambdaE':None, 'Lx':None, 'Ly':None, 'Lz':None, 'b':None, 'c':None, 'Eax':None, 'Eay':None, 'Eaz':None, 'Ebx':None, 'Eby':None, 'Ebz':None, 'Ecx':None, 'Ecy':None, 'Ecz':None, 'ovdens':None, 'nbins':None, 'fMhires':None, 'Ekin':None, 'Epot':None, 'SurfP':None, 'Phi0':None, 'cNFW':None }
        elif numberColumns!=noDataColumns and numberColumns!=None:
            self.noColumns = numberColumns
            throwWarning( "Match problem between the number of input columns and the column names stored in the 'AHFHalos' class." )
        self.nameList = self.name.values()
        if data is not None:
            self.AddData(data)
        if header is not None:
            self.AddFileHeader(header)
    
    def AddData(self,Data,hasId=True):
        self.hasData = True
        self.noHalos = Data.shape[0]
        self.data = None
        if hasId:
            if Data.shape[1]!=self.noColumns: throwError( "The argument 'Data' in function 'AHFHalos:AddData' must have %i columns but it has %i columns." % (len(self.name),Data.shape[1]) )
            self.data = Data.astype(np.float32)
        else:
            if Data.shape[1]!=self.noColumns-1: throwError( "The argument 'Data' in function 'AHFHalos:AddData' must have %i columns but it has %i columns." % (len(self.name)-1,Data.shape[1]) )
            __id = np.arange( self.noHalos ).astype(np.uint32)
            self.data = np.column_stack( (__id,Data) ).astype(np.float32).copy()
        self.shape = self.data.shape
        for i in range(len(self.name)):
            self.array[ self.name[i] ] = self.data[:,i]
        self.CleanColumns()
    
    def CreateVectors(self):
        if (self.array['Xc'] is not None) and (self.array['Yc'] is not None) and (self.array['Zc'] is not None):
            self.pos = np.column_stack( (self.array['Xc'],self.array['Yc'],self.array['Zc']) )
            self.array['pos'] = self.pos
        if (self.array['VXc'] is not None) and (self.array['VYc'] is not None) and (self.array['VZc'] is not None):
            self.vel = np.column_stack( (self.array['VXc'],self.array['VYc'],self.array['VZc']) )
            self.array['vel'] = self.vel
        if (self.array['Lx'] is not None) and (self.array['Ly'] is not None) and (self.array['Lz'] is not None):
            self.L = np.column_stack( (self.array['Lx'],self.array['Ly'],self.array['Lz']) )
            self.array['L'] = self.L
    
    def CleanColumns(self):
        self.array['id'] = self.data[:,0].astype(np.uint32)
        self.array['npart'] = self.data[:,1].astype(np.uint32)
        self.array['nbins'] = self.data[:,33].astype(np.uint32)
        self.CreateVectors()
    
    def AddFileHeader(self,header):
        self.hasFileHeader = True
        self.header = header
    
    def Description(self):
        __description = [ ('id',"the halo id (unsigned integer)"),
        ('npart',"the number of particles in the halo (unsigned integer)"),
        ('nvpart',"the number of virtual particles ?"),
        ('Xc',"the x position of the halo center (Mpc/h)"),
        ('Yc',"the y position of the halo center (Mpc/h)"),
        ('Zc',"the z position of the halo center (Mpc/h)"),
        ('VXc',"the peculiar velocity of the halo center along the x direction (km/s)"),
        ('VYc',"the peculiar velocity of the halo center along the y direction (km/s)"),
        ('VZc',"the peculiar velocity of the halo center along the z direction (km/s)"),
        ('Mvir',"the virial mass of the halo (M_0/h)"),
        ('Rvir',"the virial radius of the halo (kpc/h)"),
        ('Vmax',"maximum of the rotation curve of the halo (km/s)"),
        ('Rmax',"position of velocity maximum of the halo (kpc/h)"),
        ('sigV',"3D velocity dispersion (km/s)"),
        ('lambda',"spin parameter (Bullock+ 2001 for definition)"),
        ('Lx',"orientation along x of angular momentum (for |L|=1)"),
        ('Ly',"orientation along y of angular momentum (for |L|=1)"),
        ('Lz',"orientation along z of angular momentum (for |L|=1)"),
        ('a',"largest axis (derived from inertia tensor, normalized to unity)"),
        ('Eax',"orientation of axis a along x (for |a|=1)"),
        ('Eay',"orientation of axis a along y (for |a|=1)"),
        ('Eaz',"orientation of axis a along z (for |a|=1)"),
        ('b',"second largest axis (value here = b/a)"),
        ('Ebx',"orientation of axis b along x (for |b|=1)"),
        ('Eby',"orientation of axis b along y (for |b|=1)"),
        ('Ebz',"orientation of axis b along z (for |b|=1)"),
        ('c',"third largest axis (value here = c/a)"),
        ('Ecx',"orientation of axis c along x (for |c|=1)"),
        ('Ecy',"orientation of axis c along y (for |c|=1)"),
        ('Ecz',"orientation of axis c along z (for |c|=1)"),
        ('ovdens',"overdensity at virial radius"),
        ('Redge',"actual edge of the halo (not implemented)"),
        ('nbins',"number of bins used for the halo profile computation"),
        ('Ekin',"kinetic energy (M_0/h (km/s)^2)"),
        ('Epot',"potential energy (M_0/h (km/s)^2)"),
        ('mbp_offset',"offset between most bound particle and halo center (kpc/h)"),
        ('com_offset',"offset between center of mass and halo center (kpc/h)"),
        ('r2',"position where ('rho r^2' peaks (kpc/h)"),
        ('lambdaE',"spin parameter - Peebles' definition"),
        ('v_esc',"the escape velocity at Rvir (km/s)"),
        ('Phi0',"phi_0 from the gravitational unbinding procedure ((km/s)^2)") ]
        return __description




def numberOfLines(filename):
    """ Counts the number of lines in a text file. """
    lines = 0
    for line in open(filename):
        lines += 1
    return lines


def AHFMultipleFiles(splitFilename,fileNumber,fileIndex):
    temp = str(fileNumber)[1:]+'.'
    if len(splitFilename)==2:
        return splitFilename[0] + temp + splitFilename[1]
    elif len(splitFilename)==1 and fileIndex==0:
        return splitFilename[0]
    else:
        return 'UnknowFileName'


def numberOfAHFHalos_modified(splitFilename,largeNumber,maxFileCount):
    if maxFileCount>largeNumber:
        throwError( "Error in function 'numberOfAHFHalos_modified'. The expected number of files is larger than the root used to compose the individual files." )
    noHalosInFile = {}
    noTotalHalos = 0
    noValidFiles = 0
    while 1:
        fileNumber = largeNumber + noValidFiles
        tempName = AHFMultipleFiles(splitFilename,fileNumber,noValidFiles)
        if os.path.isfile(tempName):
            noHalosInFile[noValidFiles] = numberOfLines(tempName)
            if noValidFiles is 0: noHalosInFile[noValidFiles] -= 1
            noTotalHalos += noHalosInFile[noValidFiles]
        elif noValidFiles<maxFileCount:
            noHalosInFile[noValidFiles] = 0
        noValidFiles += 1
        if noValidFiles>=maxFileCount:
            break
        
    return [noValidFiles,noTotalHalos,noHalosInFile]
    
def numberOfAHFHalos(splitFilename,largeNumber,maxFileCount=None):
    if maxFileCount:
        return numberOfAHFHalos_modified(splitFilename,largeNumber,maxFileCount)
    noHalosInFile = {}
    noTotalHalos = 0
    noValidFiles = 0
    while 1:
        fileNumber = largeNumber + noValidFiles
        tempName = AHFMultipleFiles(splitFilename,fileNumber,noValidFiles)
        if os.path.isfile(tempName):
            noHalosInFile[noValidFiles] = numberOfLines(tempName)
            if noValidFiles is 0: noHalosInFile[noValidFiles] -= 1
            noTotalHalos += noHalosInFile[noValidFiles]
        else: break
        noValidFiles += 1
    return [noValidFiles,noTotalHalos,noHalosInFile]



def readAHFHalos_textFile(file,numberZeros=3,VERBOSE=True,maxFileCount=None):
    """ Reads the AHF halo properties from the output AHF file '.halos'.
    If there are multiple output files than the first file should be 000, the next 001 and so on. One can choose how many such 000's there are using the function parameter 'numberZeros'.
    The argument 'file' is always the filename for the first file in multiple files (i.e. 000)."""
    
    # create the string used to detect the results stored in multiple files 
    largeNumber = 10**numberZeros
    zeroStr = str(largeNumber)[1:] #now this is '0..0' numberZeros
    zeroStr += '.'
    splitFilename = file.split(zeroStr)
    if len(splitFilename)==1:
        #throwWarning( "Couldn't split the supplied AHF '_halos' filename to be able to read multiple files (if that is the case). Please check the argument 'numberZeros'=%i of function 'readAHFHalos' and see if it has the correct value."  % numberZeros)
        print "Single AHF file is: '%s'" % splitFilename[0]
    elif len(splitFilename)!=2: throwError( "Split the AHF '_halos' filename using the '%s' string in %i parts. There was an error somewhere since the program expected only 2 parts." % (zeroStr,len(splitFilename)) )
    
    # count how many lines there are in the AHF files, this tells the program the number of halos
    [noValidFiles,noTotalHalos,noHalosInFile] = numberOfAHFHalos(splitFilename,largeNumber,maxFileCount)
    if VERBOSE: print "Found %i valid AHF '_halos' file(s) with %i halos." % (noValidFiles,noTotalHalos)
    
    # read the halos properties in each file
    noColumns = noColumnsAHFTextFile
    columnsDescription = None
    data = None
    for i in range(noValidFiles):
        if noHalosInFile[i]<=0: continue
        tempName = AHFMultipleFiles( splitFilename, largeNumber+i, i )
        if VERBOSE: print "Reading file %i of %i files with name '%s' ... " % (i+1,noValidFiles,tempName),
        f = open( tempName, 'r' )
        if i is 0:
            columnsDescription = f.readline()    #skip first line which contains the names of the columns (only for the first file)
            noColumns = len( columnsDescription.split() )
            print "\t found that the data has %i columns ..." % noColumns
            data = np.append( [], np.fromfile(f,np.float64,noColumns*noHalosInFile[i],sep=" ") )
        else:
            data = np.append( data, np.fromfile(f,np.float64,noColumns*noHalosInFile[i],sep=" ") )
        f.close()
        if VERBOSE: print "Done."
    
    # copy the data into the 'AHFHalos' class
    data.shape = (data.size/noColumns,noColumns)
    data = sortHalos( data, sortColumn=data[:,0],order='descending',VERBOSE=VERBOSE)
    hasID = [False, True][ int(maxFileCount!=None) ]
    if hasID:
        data[:,0] = np.arange( data.shape[0] )
    noColumns2 = [noColumns+1,noColumns][ int(hasID) ]
    halos = AHFHalos( numberColumns=noColumns2 )
    halos.AddData( data, hasId=hasID )
    halos.AddFileHeader(columnsDescription)
    
    return halos


def writeAHFHalos_textFile(file,halos,description=None,VERBOSE=True):
    """ Writes the AHF halo properties to a single text file. """
    
    if VERBOSE: print "Writting the AHF halos information to the file '%s' ... " % file,
    
    #use C++ to write the following data in a text filename
    ahfOutputCode = """
    #line 1000 "ahfOutputCode"
    std::fstream  f;
    f.open( file.c_str(), std::ios::out );
    f << description;
    for (int i=0; i<noHalos; ++i)
    {
        f << std::setw(10) << std::setprecision(10)  << data(i,0);
        for (int j=1; j<noColumns; ++j)
            f << "  " << std::setw(12) << std::setprecision(10) << data(i,j);
        f << "\\n";
    }
    f.close();
    """
    
    if halos.hasFileHeader and description==None:
        description = halos.header
    data = halos.data[1:,:] .copy()    #do not write the ids
    (noHalos,noColumns) = data.shape
    inline( ahfOutputCode, ['file','description','noHalos','noColumns','data'], type_converters=converters.blitz, headers=['<fstream>','<iomanip>'] )
    if VERBOSE: print "Done."


sqliteTableNames = [ "description", "halos" ] # 'description' stores a description of the data, while 'halos' stores the halo data
def readAHFHalos_sqlite(file,VERBOSE=True):
    """ Reads the AHF halo properties from a sqlite database. """
    if sqlite3ImportError is not None:
        throwError( sqlite3ImportError )
    
    if VERBOSE: print "Reading the AHF halos from the sqlite database '%s':" % (file)
    # open the database
    connection = sqlite.connect(file)
    cursor = connection.cursor()
    
    # query the database to get the halo data
    query = "SELECT %s FROM %s" % ( ",".join(AHFHalos().nameList), sqliteTableNames[1] )
    cursor.execute( query )
    data = np.array( cursor.fetchall() )
    cursor.execute( "SELECT variable,description FROM %s" % sqliteTableNames[0] )
    header = cursor.fetchall()
    cursor.close()
    connection.close()
    
    # copy the data into the 'AHFHalos' class
    if VERBOSE: print "\tFound %i halos." % data.shape[0]
    #data = sortHalos( data, sortColumn=data[:,0],order='descending',VERBOSE=VERBOSE)
    halos = AHFHalos()
    halos.AddData( data )
    halos.AddFileHeader( header )
    
    return halos

def writeAHFHalos_sqlite(file,halos,description=None,VERBOSE=True):
    """ Writes the AHF halo properties to a sqlite database. """
    if sqlite3ImportError is not None:
        throwError( sqlite3ImportError )
    
    if VERBOSE: print "Writing the AHF halos to the sqlite database '%s':" % (file)
    # check that the data is fine
    haloData = None
    hasValidHeader = False
    if halos.__class__ is AHFHalos().__class__:
        haloData = halos.data
        if halos.hasFileHeader and isinstance(halos.header,list): hasValidHeader = True
    else:
        haloData = halos
    __columns = AHFHalos().nameList
    if haloData.shape[1]!=len(__columns):
        throwError( "The 'halo' input argument to the function 'writeAHFHalos_sqlite' must be a 2D numpy array with %i columns. The input array has %i columns." %(len(__columns),haloData.shape[1]) )
    
    # open the database
    connection = sqlite.connect(file)
    cursor = connection.cursor()
    
    # query the database to create the tables that will store the data
    descTable = "CREATE TABLE IF NOT EXISTS %s (variable TEXT PRIMARY KEY, description TEXT)" % sqliteTableNames[0]
    haloTableColumns = "id INT PRIMARY KEY"
    for i in range(1,len(__columns)): haloTableColumns += ", %s REAL" % __columns[i]
    haloTable = "CREATE TABLE IF NOT EXISTS %s (%s)" % (sqliteTableNames[1],haloTableColumns)
    cursor.execute(descTable)
    cursor.execute(haloTable)
    
    # write the data into the description table
    descInsert = "INSERT OR REPLACE INTO %s VALUES (?,?)" % sqliteTableNames[0]
    if description!=None and isinstance(description,list):
        try: cursor.executemany(descInsert,description)
        except: throwError( "When inserting the 'description' argument of function 'writeAHFHalos_sqlite' into the database. The 'description' argument must be a 2D matrix, with 2 columns, giving the parameter name and its description. You need python2.6 or newer for this to work." )
    if hasValidHeader: cursor.executemany(descInsert,halos.header)
    cursor.executemany(descInsert,halos.Description())
    
    # write the halo data
    haloTableColumns = "?"
    for i in range(1,len(__columns)): haloTableColumns += ", ?"
    haloInsert = "INSERT OR REPLACE INTO %s VALUES (%s)" % (sqliteTableNames[1],haloTableColumns)
    cursor.executemany(haloInsert,haloData)
    
    # close the database
    connection.commit()
    cursor.close()
    connection.close()
    if VERBOSE: print "\tDone"


def readAHFHalos_binary(file,VERBOSE=True):
    """ Reads the AHF halo data from a binary file. All the data is dumped in a binary file that has the same header as the gadget file.
    The 'npart' field of the Gadget Header gives the number of halos using: [ halosInFile, 0, 0, 0, 0, 0]  """
    
    if VERBOSE: print "Reading the AHF halo data from the binary file '%s' ... " % (file)
    # read the header
    header = GadgetHeader()
    f = open(file, 'rb')
    header.fromfile( f )
    
    # read the data
    noHalos = header.npart[0]
    noColumns = noDataColumns
    if header.npart[1]>0: noColumns = header.npart[1]
    __buffer1 = np.fromfile( f, np.dtype('i4'), 1 )[0]
    if noHalos!=0: data = np.fromfile( f, np.dtype('f4'), noColumns*noHalos )
    else: data = np.zeros( 0, np.dtype('f4') )
    __buffer2 = np.fromfile( f, np.dtype('i4'), 1 )[0]
    f.close()
    if __buffer1!=__buffer2: throwError( "While reading the AHF halo data from the binary file '%s'" % file )
    
    # write the data to the AHFHalos structure
    data.shape = -1, noColumns
    halos = AHFHalos( data=data, header=header, numberColumns=noColumns )
    return halos

def writeAHFHalos_binary(file,halos,VERBOSE=True):
    """ Writes the AHF halo data to a binary file. All the data is dumped in a binary file that has the same header as the gadget file.
    The 'npart' field of the Gadget Header gives the number of halos using: [ halosInFile, 0, 0, 0, 0, 0]  """
    
    if VERBOSE: print "Writing the AHF halo data to the binary file '%s' ... " % (file)
    # write the header
    header = GadgetHeader()
    data = halos
    if halos.__class__ is AHFHalos().__class__:
        if halos.header.__class__ is GadgetHeader().__class__: header = halos.header
        data = halos.data
    if data.shape[1]!=halos.noColumns:
        throwError( "The AHF halo data supplied to the 'writeAHFHalos_binary' function must be a 2D numpy matrix with %i columns. But the data has %i columns." % (noDataColumns,data.shape[1]) )
    if data.dtype!=np.float32:
        throwError( "The type of the AHF halo data must be %s, but the input data is of type %s." %('np.float32',data.dtype) )
    header.npart[:] = data.shape[0], data.shape[1], 0, 0, 0, 0
    f = open(file, 'wb')
    header.tofile( f )
    
    # write the data
    noBytes = data.nbytes
    __buffer = np.array( [noBytes], dtype=np.uint32 )
    __buffer.tofile( f )
    if noBytes!=0: data.tofile( f )
    __buffer.tofile( f )
    f.close()



def readAHFHalos(file,fileType='binary',numberZeros=3,VERBOSE=True,maxFileCount=None):
    """ Reads the AHF halos from the input file. """
    fileFormats = [ 'binary', 'text_file', 'sqlite' ]
    if file.split('.')[-1] in ['AHF_halos'] and fileType!='text_file':
        fileType = 'text_file'
        print "Detected that the input file is a 'AHF_halos' text file. The program will use the text file reader to read in the AHF halo data."
    if fileType not in fileFormats:
        throwError( "When reading the AHF halos from the input file in function 'readAHFHalos'. Unknow value for argument 'fileType'. Allowed values are '%s', while the input value is '%s'" % (str(fileFormats),fileType) )
    elif fileType is 'binary':
        return readAHFHalos_binary(file,VERBOSE)
    elif fileType is 'text_file':
        return readAHFHalos_textFile(file,numberZeros,VERBOSE,maxFileCount)
    elif fileType is 'sqlite':
        return readAHFHalos_sqlite(file,VERBOSE)

def writeAHFHalos(file,halos,fileType='binary',description=None,VERBOSE=True):
    """ Writes the AHF halo properties to a single file. """
    fileFormats = [ 'binary', 'text_file', 'sqlite' ]
    if fileType not in fileFormats:
        throwError( "When writing the AHF halos to the output file in function 'writeAHFHalos'. Unknow value for argument 'fileType'. Allowed values are '%s', while the input value is '%s'" % (str(fileFormats),fileType) )
    elif fileType is 'binary':
        return writeAHFHalos_binary(file,halos,VERBOSE=VERBOSE)
    elif fileType is 'text_file':
        return writeAHFHalos_textFile(file,halos,description=description,VERBOSE=VERBOSE)
    elif fileType is 'sqlite':
        return writeAHFHalos_sqlite(file,halos,description=description,VERBOSE=VERBOSE)


#
#  Functions that operate on AHF halo data
#

def sortHalos(Data,sortColumn,order='ascending',VERBOSE=True):
    """ Function that sorts the rows of a matrix according to the values given in the 'sortColumn'. The order can be 'ascending' or 'descending'.
    The result is written to the input array. """
    if VERBOSE: print "Sorting the data in %s order ... " % order,
    
    # check the variables for consistency
    data = None
    if Data.__class__ is AHFHalos().__class__:
        data = Data.data
    elif Data.__class__ is np.zeros(1,'i').__class__:
        data = Data
    else: throwError( "Unknown type for argument 'Data' in function 'sortHalos'. Recognized types are: 'AHFHalos' and 'numpy arrays'. The type of the input data is '%s'." % str(Data.__class__) )
    orderValues = ['ascending','descending']
    if order not in orderValues: throwError( "The argument 'order' of function 'sortHalos' can take only the values %s." % str(orderValues) )
    if sortColumn.ndim is not 1: throwError( "The argument 'sortColumn' of function 'sortHalos' muts be a 1D numpy array whose values are used for sorting the halo properties." )
    dataIs1D = False
    if data.ndim is 1:
        data.shape = (-1,1)
        dataIs1D = True
    elif data.ndim is not 2:
        throwError( "The data array supplied to function 'sortHalos' must be a 1D or 2D numpy array." )
    noHalos, noColumns = data.shape
    if sortColumn.ndim is not 1:
        throwError( "The 'sortColumn' argument of function 'sortHalos' must be a 1D numpy array giving the values in the sorting column corresponding to the data in the same order." )
    if noHalos!=sortColumn.size:
        throwWarning( "There are %i values in the 'sortColumn' while there are %i rows in tha data matrix. Will sort only the minimum of the two values." % (sortColumn.size,noHalos) )
        noHalos = min(sortColumn.size,noHalos)
    
    # use C++ to sort the halos
    sortHalosCode = """
    #line 1000 "sortHalosCode"
    // sort the entries according to the values in 'sortColumn'
    std::vector< std::pair<double,int> > tempSort;
    tempSort.reserve(noHalos);
    for (int i=0; i<noHalos; ++i)
        tempSort.push_back( std::make_pair(sortColumn(i),i) );
    std::sort( tempSort.begin(), tempSort.end() );
    if (order==1) std::reverse( tempSort.begin(), tempSort.end() );
    
    // rewrite the entries in the 'data' matrix according to the new order
    std::vector<int> newPosition( noHalos, -1 );
    for (size_t i=0; i<tempSort.size(); ++i)
        newPosition[ tempSort[i].second ] = i;  //the row at 'tempSort[i].second' needs to go at position 'i'
    // create temporary matrix to store the values in the correct order
    double **temp;
    temp = new double*[noHalos];
    for (int i=0; i<noHalos; ++i)
        temp[i] = new double[noColumns];
    for (int i=0; i<noHalos; ++i)
        for (int j=0; j<noColumns; ++j)
            temp[newPosition[i]][j] = data(i,j);
    // now copy back the values to the 'data' array
    for (int i=0; i<noHalos; ++i)
        for (int j=0; j<noColumns; ++j)
            data(i,j) = temp[i][j];
    
    // free the memory
    for (int i=0; i<noHalos; ++i)
        delete[] temp[i];
    delete[] temp;
    """
    if order is 'ascending': order = 0
    elif order is 'descending': order = 1
    
    inline( sortHalosCode, ['noHalos','noColumns','data','sortColumn', 'order'], type_converters=converters.blitz, headers=['<new>','<vector>','<utility>','<algorithm>'] )
    if VERBOSE: print "Done."
    
    #return the results
    if Data.__class__ is AHFHalos().__class__:
        Data.AddData(data)
        return Data
    elif Data.__class__ is np.zeros(1,'i').__class__:
        return data


def halosInMassRange(Data,massColumn,minMass,maxMass,VERBOSE=True):
    """ Returns the halos with masses in the range 'minMass' to 'maxMass'. The masses are given by the argument 'massColumn'.
    The result is copied to a new array. """
    if VERBOSE: print "Selecting the halos with masses in the interval [%e,%e] ... " % (minMass,maxMass)
    
    #do some error checking
    data = None
    if Data.__class__ is AHFHalos().__class__:
        data = Data.data
    elif Data.__class__ is np.zeros(1,'i').__class__:
        data = Data
    else: throwError( "Unknown type for argument 'Data' in function 'halosInMassRange'. Recognized types are: 'AHFHalos' and 'numpy arrays'. The type of the input data is '%s'." % str(Data.__class__) )
    if data.ndim!=1 and data.ndim!=2:
        throwError( "The data array supplied to function 'halosInMassRange' must be a 1D or 2D numpy array." )
    if massColumn.ndim!=1:
        throwError( "The 'massColumn' array supplied to function 'halosInMassRange' must be a 1D numpy array." )
    if massColumn.shape[0]!=data.shape[0]:
        throwError( "The 'massColumn' and 'data' arrays supplied to function 'halosInMassRange' must have the same length of the first dimension." )
    
    # get the haloes in the mass range
    validEntries = (massColumn>=minMass) * (massColumn<=maxMass)
    result = data[validEntries,:]
    if VERBOSE: print "Done."
    
    #return the results
    if Data.__class__ is AHFHalos().__class__:
        returnHalos = AHFHalos(data=result,header=Data.header)
        return returnHalos
    elif Data.__class__ is np.zeros(1,'i').__class__:
        return result



def halosInMask(Data,pos,boxLength,mask,VERBOSE=True):
    """Returns two new arrays which contains the halos in the mask and the ones outside the mask. Valid regions of the mask are those with 'values>=1'.
            Data - full halo data
            pos - aray giving halo position in the same order as in Data
            boxLength - 6 coordinates giving xMin, xMax, yMin, yMax, zMin, zMax for the mask box
            mask - integer array with values >0 for the valid voxels for which to select the halos.
    return halosInMask, halosOutsideMask """
    if VERBOSE: print "Finding the halos in the mask ..."
    
    #do some error checking
    data = None
    if Data.__class__ is AHFHalos().__class__:
        data = Data.data
    elif Data.__class__ is np.zeros(1,'i').__class__:
        data = Data
    else: throwError( "Unknown type for argument 'Data' in function 'halosInMask'. Recognized types are: 'AHFHalos' and 'numpy arrays'. The type of the input data is '%s'." % str(Data.__class__) )
    if data.ndim>2:
        throwError( "The data array supplied to function 'halosInMask' must be a 1D or 2D numpy array." )
    if not (pos.ndim==2 and pos.shape[1]==3):
        throwError( "The 'pos' array supplied to function 'halosInMask' must be a 2D numpy array giving the position of the halos in the same order as in the data array. The second dimension of the 'pos' array must be 3." )
    if pos.shape[0]!=data.shape[0]:
        throwError( "The 'pos' and 'data' arrays supplied to function 'halosInMask' must have the same length of the first dimension." )
    if mask.ndim!=3:
        throwError( "The 'mask' array supplied to function 'halosInMask' must be a 3D integer numpy array giving the mask from where to select the valid halos." )
    
    
    # Use C++ code to rapidly find the halos in the mask and outside the mask
    halosInMaskCode = """
    #line 1000 "halosInMaskCode"
    int n[] = { mask.extent(0), mask.extent(1), mask.extent(2) };
    float dx[3];
    for (int i=0; i<3; ++i) dx[i] = (boxLength(2*i+1)-boxLength(2*i)) / n[i];
    int noHalos = data.extent(0), noColumns = data.extent(1);
    int noInMask = 0, noOutsideMask = 0;
    for (int i=0; i<noHalos; ++i)
    {
        bool validHalo = true;
        int temp[3];
        for (int j=0; j<3; ++j)
        {
            temp[j] = int(std::floor( (pos(i,j)-boxLength(2*j)) / dx[j] ));
            if (temp[j]<0 or temp[j]>=n[j])
            {
                validHalo = false;
                break;
            }
        }
        if ( validHalo and mask(temp[0],temp[1],temp[2]) )
        {
            for (int j=0; j<noColumns; ++j)
                inmask(noInMask,j) = data(i,j);
            ++noInMask;
        }
        else
        {
            for (int j=0; j<noColumns; ++j)
                outsidemask(noOutsideMask,j) = data(i,j);
            ++noOutsideMask;
        }
    }
    py::tuple results(2);
    results[0] = noInMask;
    results[1] = noOutsideMask;
    return_val = results;
    """
    
    inmask = np.zeros( data.shape, data.dtype )
    outsidemask = np.zeros( data.shape, data.dtype )
    noInMask, noOutsideMask = inline( halosInMaskCode, ['data','pos','boxLength','mask','inmask', 'outsidemask'], type_converters=converters.blitz, headers=['<cmath>'] )
    inmask.resize( (noInMask,data.shape[1]), refcheck=False )
    outsidemask.resize( (noOutsideMask,data.shape[1]), refcheck=False )
    if VERBOSE:
        noHalos = data.shape[0]
        print "\t found %i (%.1f%%) inside the mask and %i (%.1f%%) outside the mask." % (noInMask,float(noInMask)/noHalos*100.,noOutsideMask,float(noOutsideMask)/noHalos*100.)
    
    #return the results
    if Data.__class__ is AHFHalos().__class__:
        halosInMask = AHFHalos(data=inmask,header=Data.header)
        halosOutsideMask = AHFHalos(data=outsidemask,header=Data.header)
        return halosInMask, halosOutsideMask
    elif Data.__class__ is np.zeros(1,'i').__class__:
        return inmask, outsidemask



def halosSubhalos(Data,pos,radius,boxLength,VERBOSE=True):
    """ This function returns two different arrays with the halos and subhalos.
        Data - full data array
        pos - halo position matrix (same order as Data array)
        radius - halo radius matrix (same order as Data array)
        boxLength - array giving the periodic box length along each direction
    """
    if VERBOSE: print "Finding the AHF halos and subhalos ..."
    
    #do some error checking
    data = None
    if Data.__class__ is AHFHalos().__class__:
        data = Data.data
    elif Data.__class__ is np.zeros(1,'i').__class__:
        data = Data
    else: throwError( "Unknown type for argument 'Data' in function 'halosSubhalos'. Recognized types are: 'AHFHalos' and 'numpy arrays'. The type of the input data is '%s'." % str(Data.__class__) )
    if data.ndim!=2:
        throwError( "The 'Data' array supplied to function 'halosSubhalos' must be a 2D numpy array." )
    if pos.ndim!=2 or pos.shape[1]!=3:
        throwError( "The 'pos' array supplied to function 'halosSubhalos' must be a 2D numpy array giving the position of the halos in the same order as in the data array. The second dimension of the 'pos' array must be 3." )
    if pos.shape[0]!=data.shape[0]:
        throwError( "The 'pos' and 'data' arrays supplied to function 'halosSubhalos' must have the same length of the first dimension." )
    if radius.ndim!=1:
        throwError( "The 'radius' array supplied to function 'halosSubhalos' must be a 1D numpy array giving the radius of the halos." )
    if radius.shape[0]!=data.shape[0]:
        throwError( "The 'radius' and 'data' arrays supplied to function 'halosSubhalos' must have the same length of the first dimension." )
    if boxLength.size!=3:
        throwError( "The 'boxLength' argument supplied to function 'halosSubhalos' must be a 1D array with 3 elements giving the periodic box lengths along each direction." )
    
    # Use C++ code to rapidly find the halos in the mask and outside the mask
    halosSubhalosSupportCode ="""
    #line 1000 "halosSubhalosSupportCode"
    float periodicDistance(float *pos1, float *pos2, float *length)
    {
        float temp, dist = 0.;
        for (int i=0; i<3; ++i)
        {
            temp = pos1[i] - pos2[i];
            if ( std::fabs(temp)>std::fabs(temp+length[i]) ) temp += length[i];
            else if ( std::fabs(temp)>std::fabs(temp-length[i]) ) temp -= length[i];
            dist += temp*temp;
        }
        return dist;
    }
    """
    
    halosSubhalosCode = """
    #line 1000 "halosSubhalosCode"
    int const noObjects = data.extent(0);
    int const noColumns = data.extent(1);
    int noHalos = 0, noSubhalos = 0;
    bool isSubhalo[noObjects];
    for (int i=0; i<noObjects; ++i) isSubhalo[i] = false;
    float length[] = {boxLength(0), boxLength(1), boxLength(2)};
    
    // find which are halos and subhalos by checking if distance < sum of radia (more massive one is the halo)
    for (int i=0; i<noObjects; ++i)
    {
        float pos1[] = {pos(i,0), pos(i,1), pos(i,2)};
        for (int j=0; j<i; ++j)
        {
            float pos2[] = {pos(j,0), pos(j,1), pos(j,2)};
            float d = periodicDistance( pos1, pos2, length );
            float d2 = radius(i)+radius(j);
            d2 = d2*d2;
            if ( d2>d )
            {
                isSubhalo[i] = true;
                break;
            }
        }
    }
    
    // copy the halo and subhalo data to the two arrays
    for (int i=0; i<noObjects; ++i)
    {
        if ( isSubhalo[i] )
        {
            for (int j=0; j<noColumns; ++j)
                subhalos(noSubhalos,j) = data(i,j);
            ++noSubhalos;
        }
        else
        {
            for (int j=0; j<noColumns; ++j)
                halos(noHalos,j) = data(i,j);
            ++noHalos;
        }
    }
    
    py::tuple results(2);
    results[0] = noHalos;
    results[1] = noSubhalos;
    return_val = results;
    """
    
    
    halos, subhalos = None, None
    noHalos, noSubhalos = None, None
    if Data.__class__ is AHFHalos().__class__ and 'hostHalo' in Data.nameList:    #the halo properties keep track of which are halo and which subhalos
        if VERBOSE: print "\t finding the halos and subhalos using the host halo information from the halo properties ..."
        select = Data.array['hostHalo']==-1
        halos  = data[select]
        subhalos = data[-select]
        noHalos, noSubhalos = halos.shape[0], subhalos.shape[0]
    else:
        if VERBOSE: print "\t finding the halos and subhalos by searching all the objects within the same radius ..."
        halos = np.zeros( data.shape, data.dtype )
        subhalos = np.zeros( data.shape, data.dtype )
        noHalos, noSubhalos = inline( halosSubhalosCode, ['data','pos','radius','halos', 'subhalos', 'boxLength'], type_converters=converters.blitz, support_code=halosSubhalosSupportCode )
        halos.resize( (noHalos,data.shape[1]), refcheck=False )
        subhalos.resize( (noSubhalos,data.shape[1]), refcheck=False )
    
    if VERBOSE:
        noObjects = data.shape[0]
        print "\t out of the total %i objects: %i (%.1f%%) are halos and %i (%.1f%%) are subhalos" % (noObjects,noHalos,float(noHalos)/noObjects*100.,noSubhalos,float(noSubhalos)/noObjects*100.)
    
    #return the results
    if Data.__class__ is AHFHalos().__class__:
        Halos = AHFHalos(data=halos,header=Data.header)
        Subhalos = AHFHalos(data=subhalos,header=Data.header)
        return Halos, Subhalos
    elif Data.__class__ is np.zeros(1,'i').__class__:
        return halos, subhalos
    



def computeStatistics(data,massBins,massBinValues,noBootstrapSamples=100):
    """Computes the statistics for a given distribution of halo properties. The data is stored in columns where the first column is the mass of the halo, while the rest are halo properties for which the program computes the median, error associated to the median, the 16th and 84th percentile, the average and error of the average."""
    noData = data.shape[1] - 1
        
    averages = np.zeros( (massBins.size-1,6*noData+2), np.float64 )
    averages[:,0] = massBinValues
    for i in range(1,massBins.size):    # loop over the mass bins and get median formation redshif for each bin
        selection = (data[:,0]>=massBins[i-1]) * (data[:,0]<massBins[i])    #select only the halos with mass in the given mass bin
        averages[i-1,1] = np.sum( selection )
        if averages[i-1,1]<=0: continue
        for j in range(0,noData):
            percentile = analysis.Percentile( data[selection,j+1], percentile=(16.,84.), N=1 )
            averages[i-1,6*j+2:6*j+4] = analysis.Median( data[selection,j+1], noBootstrapSamples )
            averages[i-1,6*j+4:6*j+6] = percentile[0][0], percentile[1][0]
            averages[i-1,6*j+6:6*j+8] = analysis.Average( data[selection,j+1], noBootstrapSamples )
    return averages


def computeDistributions(data, massBins,massBinValues, distBins,distBinValues,bin_type='linear'):
    """Computes the distribution in halo properties for the given halo data. The data needs to be in two columns, with the first giving the halo mass and the second the halo property. It outputs in columns: the 'distBinValues' values, the total number of halos for each 'distBins' bin, the PDF associated to it. Than, for each mass bin it outputs the total number of halos and the PDF for the halos in that mass bin as a function of 'distBins' bin. """
    distribution = np.zeros( (distBins.size-1,2*massBins.size+1), np.float32 )
    distribution[:,0] = distBinValues
    binNormalization = None
    if bin_type=='linear': binNormalization = distBins[1:] - distBins[:-1]
    elif bin_type=='logarithm': binNormalization = np.log10(distBins[1:]) - np.log10(distBins[:-1])
    else: throwError( "Unrecognized value for the 'bin_type' parameter in function 'computeDistributions'. The inserted value is '%s' while the allowed values are '%s'!" % (bin_type,['linear','logarithm']) )
    
    # get the distribution for all masses
    distribution[:,1], discard = np.histogram( data[:,1], bins=distBins )
    if distribution[:,1].sum()>0:
        distribution[:,2] = distribution[:,1] / ( binNormalization * np.sum(distribution[:,1]) )
    for i in range(1,massBins.size):    # loop over the mass bins
        selection = (data[:,0]>=massBins[i-1]) * (data[:,0]<massBins[i])    #select only the halos with mass in the given mass bin
        if np.sum( selection )<=0: continue
        distribution[:,2*i+1], discard = np.histogram( data[selection,1], bins=distBins )
        if distribution[:,2*i+1].sum()>0:
            distribution[:,2*i+2] = distribution[:,2*i+1] / ( binNormalization * np.sum(distribution[:,1]) )
    return distribution


def computeDistributions_direction(data, massBins,massBinValues, noAngleBins=10,testDistribution=None, noRepeats=100):
    """Computes the distribution in halo properties for the given halo data. The data needs to be in two columns, with the first giving the halo mass and the second the cosinus of the angle between the halo vector and the environment. It outputs in columns: the cosAngle values, the total number of halos for each angle bin, the PDF associated to it. Than, for each mass bin it outputs the total number of halos and the PDF for the halos in that mass bin as a function of angle value. """
    distribution = np.zeros( (noAngleBins,3*massBins.size+1), np.float32 )
    cosBins = np.linspace( 0., 1., noAngleBins+1 )
    distribution[:,0] = (cosBins[1:] + cosBins[:-1]) /2.
    if testDistribution==None: 
        testDistribution = analysis.randomCorrelation_direction(1000000,noAngleBins,100)[0] / 1000000
    
    # get the distribution for all masses
    distribution[:,1], randomPoints, stdRandom = analysis.dataCorrelation_direction( data[:,1], noAngleBins, noRepeats=noRepeats )
    noPoints = data.shape[0]
    distribution[:,2] = distribution[:,1] / (noPoints * testDistribution)
    distribution[:,3] = stdRandom / (noPoints * testDistribution)
    for i in range(1,massBins.size):    # loop over the mass bins
        selection = (data[:,0]>=massBins[i-1]) * (data[:,0]<massBins[i])    #select only the halos with mass in the given mass bin
        noPoints = np.sum( selection )
        if noPoints<=0: continue
        distribution[:,3*i+1], random, stdRandom = analysis.dataCorrelation_direction( data[selection,1], noAngleBins, noRepeats=noRepeats )
        distribution[:,3*i+2] = distribution[:,3*i+1] / (noPoints * testDistribution)
        distribution[:,3*i+3] = stdRandom / (noPoints * testDistribution)
    return distribution
    
