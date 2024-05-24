
import numpy as np
import os.path
from scipy.weave import inline, converters
from miscellaneous import throwError, throwWarning


def readTextFile(file,noDescriptionLines,noRows,noColumns,dataType,VERBOSE=True):
    """ Used to read data from a text file.
        noDescriptionLines - number of description lines.
        noRows - number of rows to read.
        noColumns - number of columns on each row (must give all the columns in the file).
    It returns [description,data], whit data in a numpy matrix of noRows x noColumns dimension.
    """
    if not os.path.isfile(file): throwError( "Could not find the file '%s' for reading (in function 'readTextFile')." % file )
    if VERBOSE: print "Reading the data from the ASCII file '%s' ... " % file,
    f = open( file, 'r' )
    description = ''
    for i in range(noDescriptionLines):
        description += f.readline()
    dataSize = noRows*noColumns
    data = np.fromfile( f, dataType, dataSize, ' ' )
    data.shape = (noRows,noColumns)
    if VERBOSE: print "Done"
    
    return [description, data]


def writeTextFile2(file,description,data,VERBOSE=True):
    if VERBOSE: print "Writing the data to the ASCII file '%s' ... " % file
    f = open( file, 'w' )
    f.write( description )
    data.tofile(f,'   ', '%12.7g')
    f.close()

def writeTextFile(file,description,data,VERBOSE=True):
    """ Used to write a matrix to a text file.
        file - the name of the output file
        description - description lines in file
        data - the matrix to be written (a numpy array of 1 or 2 dimensions)
    """
    noRows = 1
    noColumns = 1
    if data.ndim is 1:
        (noRows,) = data.shape
        data.shape = (noRows,1)
    elif data.ndim is 2:
        (noRows,noColumns) = data.shape
    else: throwError( "The data argument supplied to function 'writeTextFile' must be a 1D or 2D numpy array." )
    
    if VERBOSE: print "Writing the data to the ASCII file '%s' ... " % file,
    # Use C++ to write the data to the file
    writeTextFileCode ="""
    #line 1000 "writeTextFileCode"
    std::fstream  f;
    f.open( file.c_str(), std::ios::out );
    f << description;
    for (int i=0; i<noRows; ++i)
    {
        f << std::setw(10) << std::setprecision(10)  << data(i,0);
        for (int j=1; j<noColumns; ++j)
            f << "  " << std::setw(12) << std::setprecision(10) << data(i,j);
        f << "\\n";
    }
    f.close();
    """
    inline( writeTextFileCode, ['file','description','noRows','noColumns','data'], type_converters=converters.blitz, headers=['<fstream>','<iomanip>'] )
    if VERBOSE: print "Done"


def writeTextFile_gnuplot3D(file,description,data,VERBOSE=True):
    """ Used to write a 3D data set to a text file.
        file - the name of the output file
        description - description lines in file
        data - the matrix to be written (a numpy array of 1 or 2 dimensions)
    """
    xSize = 1
    ySize = 1
    noColumns = 1
    if data.ndim is 3:
        xSize, ySize, noColumns = data.shape
    else: throwError( "The data argument supplied to function 'writeTextFile_gnuplot3D' must be a 3D numpy array." )
    
    if VERBOSE: print "Writing the data to the ASCII file '%s' ... " % file,
    # Use C++ to write the data to the file
    writeTextFile_gnuplot3DCode ="""
    #line 1000 "writeTextFile_gnuplot3DCode"
    std::fstream  f;
    f.open( file.c_str(), std::ios::out );
    f << description;
    for (int i=0; i<xSize; ++i)
    {
        for (int j=0; j<ySize; ++j)
        {
            f << std::setw(10) << std::setprecision(10)  << data(i,j,0);
            for (int k=1; k<noColumns; ++k)
                f << "  " << std::setw(12) << std::setprecision(10) << data(i,j,k);
            f << "\\n";
        }
        f << "\\n";
    }
    f.close();
    """
    inline( writeTextFile_gnuplot3DCode, ['file','description','xSize','ySize','noColumns','data'], type_converters=converters.blitz, headers=['<fstream>','<iomanip>'] )
    if VERBOSE: print "Done"

