
import numpy as np
import numpy.fft as fft
import math
from scipy.weave import inline, converters
from miscellaneous import throwError, throwWarning

ownCIncludePath = [  "/net/plato/data/users/cautun/Programs/stow/include/", "/net/plato/data/users/cautun/Code/MMF_C/classes/", "/net/plato/data/users/cautun/Code/MMF_C/miscellaneous/" ]
ownCLibraryPath = [ "/net/plato/data/users/cautun/Programs/stow/lib/" ]
fftwLibraryPath = "/net/plato/data/users/cautun/Programs/stow/lib/" 



def MultidimensionalArray(noDims,arrays,operation='multiply',VERBOSE=True):
    """This function takes 'noDims' 1D arrays and than uses 'operation' to rewrite the 1D arrays into a 'noDims'-D array given by:
            result(i,j,k) = arrays[0](i) 'operation' arrays[1](j) 'operation' arrays[2](k)
    where:
            noDims = gives the number of dimensions
            arrays = a list of 'noDims' 1-D arrays
            operation = ['add','multiply'] to add the elements or multiply them
    """
    funcName = 'MultidimensionalArray'
    
    if noDims>4 or noDims<1:
        throwError( "The argument 'noDims' in the '%s' function must in the interval 1 to 4. You inserted 'noDim'=%s" % (funcName,noDims) )
    if noDims!=len(arrays):
        throwError( "The length of the input list 'arrays' in the '%s' function must be the same as the number of dimensions 'noDims=%i'. The list 'arrays' has the length %i." % (funcName,noDims,len(arrays)) )
    for i in range(noDims):
        if arrays[i].ndim!=1: throwError( "All the elements of the input list arrays' in the '%s' function must be 1D numpy arrays. Element %i is a %i-D numpy array." % (funcName,i,arrays[i].ndim) )
    
    operationList = {'add':'+', 'multiply':'*'}
    if operation not in operationList.keys(): 
        throwError( "The argument 'operation' in function '%s' must have one of the values in the list '%s'. You inserted for the argument the value '%s'." % (funcName,str(operationList.keys()),operation) )
    operationText = operationList[ operation ]
    
    if noDims==1:
        return arrays[0]
    elif noDims==2:
        return eval( 'arrays[0][:,np.newaxis] %s arrays[1][np.newaxis,:]' % operationText )
    elif noDims==3:
        return eval( 'arrays[0][:,np.newaxis,np.newaxis] %s arrays[1][np.newaxis,:,np.newaxis] %s arrays[2][np.newaxis,np.newaxis,:]' % (operationText,operationText) )
    elif noDims==4:
        return eval( 'arrays[0][:,np.newaxis,np.newaxis,np.newaxis] %s arrays[1][np.newaxis,:,np.newaxis,np.newaxis] %s arrays[2][np.newaxis,np.newaxis,:,np.newaxis] %s arrays[3][np.newaxis,np.newaxis,np.newaxis,:]' % (operationText,operationText,operationText) )
    
    throwError( "Unrecognized value for the argument 'noDims' in the '%s' function must in the interval 1 to 4. You inserted 'noDim'=%s" % (funcName,noDims) )
    return None



class KSpace:
    """Class used to return the k-space values for the momenta in discrete space. The momenta values return obey the same order as for the DFT (Discrete Fourier Transform)."""
    def __init__(self,kShape,boxSize):
        if len(kShape)!=len(boxSize): throwError( "The 'KSpace' class arguments 'kShape' and 'boxSize' must have the same number of elements." )
        self.shape = np.array(kShape)
        self.boxLength = np.array(boxSize)
        self.dx = self.boxLength / self.shape
        self.scaling = 2.*np.pi / self.shape
        
        self.ndim = len(self.shape)
        self.totalGridSize = 1
        for e in self.shape: self.totalGridSize *= e
    
    def nn(self,axis):
        if axis>=self.ndim: throwError( "In class 'KSpace'. You are trying to return the n-values for axis %i when there are only %i axes." % (axis,self.ndim) )
        N = self.shape[axis]
        indices = np.arange(N)
        return np.where( indices>N/2, indices-N, indices )
    def nx(self):
        return self.nn(0)
    def ny(self):
        return self.nn(1)
    def nz(self):
        return self.nn(2)
    def nSquare(self):
        nArrays = []
        for i in range(self.ndim):
            nArrays.append( self.nn(i)**2 )
        return MultidimensionalArray(self.ndim,nArrays,operation='add',VERBOSE=False)
    def n(self):
        return np.sqrt( self.nSquare() )
    
    def kn_FT(self,axis):
        return self.nn(axis) * 2*np.pi / self.boxLength[axis]
    def kx_FT(self):
        return self.kn_FT(0)
    def ky_FT(self):
        return self.kn_FT(1)
    def kz_FT(self):
        return self.kn_FT(2)
    def kSquare_FT(self):
        nArrays = []
        for i in range(self.ndim):
            nArrays.append( self.kn_FT(i)**2 )
        return MultidimensionalArray(self.ndim,nArrays,operation='add',VERBOSE=False)
    def k_FT(self):
        return np.sqrt( self.kSquare_FT() )
    
    def kn(self,axis):
        return np.sin( self.scaling[axis] * self.nn(axis) ) / self.dx[axis]
    def kx(self):
        return self.kn(0)
    def ky(self):
        return self.kn(1)
    def kz(self):
        return self.kn(2)
    def kSquare(self):
        nArrays = []
        for i in range(self.ndim):
            nArrays.append( (np.sin( self.scaling[i]/2. * self.nn(i) ) *2./ self.dx[i])**2 )
        return MultidimensionalArray(self.ndim,nArrays,operation='add',VERBOSE=False)
    def knSquare(self,axis):
        return ( np.sin( self.scaling[axis]/2. * self.nn(axis) ) *2./ self.dx[axis])**2
    def kxSquare(self):
        return self.knSquare(0)
    def kySquare(self):
        return self.knSquare(1)
    def kzSquare(self):
        return self.knSquare(2)
    def k(self):
        return np.sqrt( self.kSquare() )


def BoxLength(noDims,boxCoordinates):
    """Computes the box length starting from the box coordinates 'boxCoordinates'=[xMin, xMax, yMin, yMax, ...]. 'noDims' gives the number of dimensions of the space. """
    boxLength = np.zeros( noDims )
    boxLength[:] = boxCoordinates[1:2*noDims:2] - boxCoordinates[0:2*noDims:2]
    return boxLength


def NumpyArray(size,elements,functionName,variableName):
    """This functions takes a tuple/list of 1 or 'size' elements and returns the corresponding numpy array of 'size' elements:
            size = the size of the numpy array returned
            elements = the elements put in the numpy array
            functionName = the name of the function calling this function - used only for error reporting purposes
            variableName = the variable name given under 'elements' - used only for error reporting purposes
    """
    result = np.zeros(size)
    if len(elements)!=size:
        if len(elements)==1: result[:] = elements[0]
        else: throwError( "The argument '%s' in function '%s' must have 1 or %i elements. But the argument has %i elements." % (variableName,functionName,size,len(elements)) )
    else:
        result[:] = elements[:]
    return result


def SolvePoissonEquation(Data,multiplicationFactor=1.,boxLength=[1.,],nMin=0.1,VERBOSE=True):
    """Solves the Poisson equation in 1D, 2D and 3D. It takes the following arguments:
            data = the input data to the solver (the right side of the equation)
            multiplicationFactor = if a multiplication factor is needed for the right side
            boxLength = 1 or n entries array giving the box size along each coordinate
            nMin = tell  how many of the low k modes to discard - (nMin>0 since divison by 0 is Nan) - give these values in terms of the sqrt(i1^2+i2^2+i3^3) modes that need to be discarded.
    """
    funcName = "SolvePoissonEquation"
    noDims = Data.ndim
    if noDims not in [1,2,3]: throwError( "The input array to the '%s' function must be a 1D, 2D or 3D numpy array. You inserted a %iD array." % (funcName,noDims) )
    boxLength = NumpyArray(noDims,boxLength,funcName,'boxLength')
    shape = Data.shape
    
    if VERBOSE: print "Solving the Poisson equation in %iD for an array of %s size ..." % (noDims,str(shape))
    momenta =  KSpace( Data.shape, boxLength )
    tempK = momenta.kSquare()
    select = momenta.nSquare() < nMin
    tempK[select] = 1.
    
    data = fft.fftn(Data)
    data *= -1.*multiplicationFactor/tempK
    data[select] = 0.
    result = fft.ifftn(data).astype(Data.dtype).real
    del data
    return result
    
    
#def HessianMatrix(Data,


def EigenvaluesEigenvectors(Data,sortOrder=None,VERBOSE=True):
    """Computes the eigenvalues and eigenvectors for an array of matrices. NOTE: the matrices need to be Hermitian. The arguments are:
            Data - the input data, at least a 2D numpy array for a single matrix or an nD (n>=3) numpy array for a multidimension array of matrices.
            sortOrder - if to sort the eigenvalues and eigenvectors. Available values: ['ascending', 'abs_ascending', 'descending', 'abs_denscending'] where the 'abs_*' methods use the absolute values of the eigenvalues for sorting.
    """
    funcName = 'EigenvaluesEigenvectors'
    shape = Data.shape
    if Data.ndim<=1: throwError( "The 'Data' argument of function '%s' needs to be at least a 2D numpy array. The dimension of 'Data' is %i." % (funcName,Data.ndim) )
    elif Data.ndim==2:
        Data.shape = 1, shape[0], shape[1]
    else:
        Data.shape = -1, shape[-2], shape[-1]   #flatten the array to a 1D array of matrices
    if shape[-2]!=shape[-1]:
        throwError( "The last two dimensions of the 'Data' array need to have the same value since the matrices need to be square." )
    
    sortingOptionsDesc = { 'ascending':'ascending eigenvalues', 'abs_ascending':'ascending eigenvalues absolute value', 'descending':'descending eigenvalues', 'abs_denscending':'descending eigenvalues absolute value' }
    sortingOptions     = { 'ascending': 'eigenVals.argsort()', 'abs_ascending':'np.abs(eigenVals).argsort()', 'descending':'(-1.*eigenVals).argsort()', 'abs_denscending':'(-1.*np.abs(eigenVals)).argsort()' }
    sorting = None
    if sortOrder!=None and sortOrder not in sortingOptions.keys():
        throwError( "Unrecognized values for the 'sortOrder' parameter in function %s. Allowed values are: %s" % (funcName,sortingOptions.keys()) )
    elif sortOrder!=None and sortOrder in sortingOptions.keys():
        sorting = sortingOptions[sortOrder]
    
    if VERBOSE:
        if sortOrder!=None: print "Selected sorting order '%s' = %s  -- in function '%s'." % (sortOrder,sortingOptionsDesc[sortOrder],funcName)
        print "Computing the eigenvalues and eigenvectors for %i matrices of dimension %ix%i ..." % (Data.shape[0],Data.shape[1],Data.shape[2])
    
    eigenValues  = np.empty( Data.shape[:-1], np.float64 )
    eigenVectors = np.empty( Data.shape, np.float64 )
    for i in range(Data.shape[0]):
        eigenVals, eigenVecs = np.linalg.eigh( Data[i,:,:] )
        if sorting is not None:
            select = eval( sorting )
            eigenVals = eigenVals[select]
            eigenVecs = eigenVecs[:,select]
        # assign the results to the output array
        eigenValues[i,:] = eigenVals
        eigenVectors[i,:,:] = eigenVecs
    
    Data.shape         = shape
    eigenValues.shape  = shape[:-1]
    eigenVectors.shape = shape
    return eigenValues, eigenVectors
    


def MultiplyArrayOnAxes(Data,axesArrays,VERBOSE=True):
    """ Takes a n-dimensional array (n=1,2 or 3) and multiplies those values with n-arrays alone for each axis.
            data = the input n-dim array
            axesArrays = a tuple of n elements giving the multiplication values 
    For example, in 2D, the code computes:
    
                result(i,j) = data(i,j) * axesArrays(0)(i) * axesArrays(1)(j)
    And returns the result.
    """
    funcName = 'MultiplyArrayOnAxes'
    noDims = Data.ndim
    if noDims not in [1,2,3]: throwError( "The input array to the 'MultiplyArrayOnAxes' function must be a 1D, 2D or 3D numpy array. You inserted a %iD array." % noDims )
    shape = Data.shape
    X1, X2, X3 = 0, 0, 0
    if len(axesArrays)!=noDims: throwError( "The input argument 'axesArrays' in function 'MultiplyArrayOnAxes' must be a tuple with %i elements (this is the dimension of the input data array). You supplied a tuple with %i elements." % (noDims,len(axesArrays) ) )
    if noDims>=1:
        X1 = axesArrays[0]
        if len(X1)!=shape[0]: throwError( "The length of the array to multiply along the x-axis is %i which is different from the input data array shape %s." % (len(X1),str(shape)) )
    if noDims>=2:
        X2 = axesArrays[1]
        if len(X2)!=shape[1]: throwError( "The length of the array to multiply along the y-axis is %i which is different from the input data array shape %s." % (len(X2),str(shape)) )
    if noDims>=3:
        X3 = axesArrays[2]
        if len(X3)!=shape[2]: throwError( "The length of the array to multiply along the z-axis is %i which is different from the input data array shape %s." % (len(X3),str(shape)) )
    
    
    if VERBOSE: print "Multiply a %iD array of shape %s with %i additional arrays, one for each axis..." % (noDims,str(shape),noDims)
    if noDims==1:
        return Data * X1
    elif noDims==2:
        temp = MultidimensionalArray( noDims, [X1,X2], 'multiply' )
        temp = Data * temp
        return temp
    elif noDims==3:
        temp = MultidimensionalArray( noDims, [X1,X2,X3], 'multiply' )
        temp = Data * temp
        return temp
    
    throwError( "Unrecognized value for the argument 'noDims' in the '%s' function must in the interval 1 to 4. You inserted 'noDim'=%s" % (funcName,noDims) )
    return None
    


def Filter(Data,filterRadius,boxLength,filterType='gaussian',VERBOSE=True):
    """ Applies a filter to the input data. The function takes the following parameters:
            Data = the input data - a n-D array (with n=1,2 or 3)
            filterRadius = a tuple/list giving the filter radius in physical units
            boxLength = a tuple/list giving the periodic box length in physical units
            filter_type = what type of smoothing to apply: 'gaussian', 'top-hat'
    The function works only for periodic data.
    """
    funcName = 'Filter'
    
    noDims = Data.ndim
    shape = Data.shape
    if noDims not in [1,2,3]: throwError( "The input array to the '%s' function must be a 1D, 2D or 3D numpy array. You inserted a %iD array." % (funcName,noDims) )
    filterRadius = NumpyArray(noDims,filterRadius,funcName,'filterRadius')
    boxLength = NumpyArray(noDims,boxLength,funcName,'boxLength')
    filterTypeList = ['gaussian','top-hat']
    if filterType not in filterTypeList: throwError( "The argument 'filter_type' in function '%s' must have one of the values in the list '%s'. You inserted for the argument the value '%s'." % (funcName,str(filterTypeList),filterType) )
    
    if VERBOSE: print "Computing the '%s' filter of radius %s applied to the input data of dimension %s ..." % (filterType,str(filterRadius),str(shape))
    filterFunction = []
    momenta = KSpace(shape,boxLength)
    if filterType is 'gaussian':
        for i in range(noDims):
            tempK = momenta.kn_FT(i)
            filterFunction.append( np.exp(-0.5 * tempK**2 * filterRadius[i]**2) )
    elif filterType is 'top-hat':
        for i in range(noDims):
            tempK = momenta.kn_FT(i)
            filterFunction.append( np.sin(np.pi * tempK * filterRadius[i]) / (np.pi*tempK*filterRadius[i]) )
            filterFunction[i][0] = 1.
    
    fftData = fft.fftn(Data)
    fftData = MultiplyArrayOnAxes(fftData,filterFunction,VERBOSE=False)
    result = fft.ifftn(fftData).real.astype(Data.dtype)
    return result


def FourierTransformSpectrum(Data,boxLength=[1.],noBins=100,bin_type='linear',VERBOSE=True):
    """ This function computes the spectrum ofthe input data by taking the Fourier Transform and bining the data according to k-values.
            Data = the input data - up to a 3D array
            boxLength = the box dimensions along each axis (1 or 2 values)
            noBins = the number of bins in the intervale kMin (=2\pi/L) to KMax (=2\pi/L N/2)
            bin_type = 'linear' or 'logarithm' specifies how to compute the k-bins used to get the spectrum
    The function returns a 2D array with 3 columns giving: k-bin value; spectrum of abs( f(k) ) and spectrum of real( f(k) )
    """
    funcName = 'FourierTransformSpectrum'
    
    noDims = Data.ndim
    if noDims not in [1,2,3]: throwError( "The input array to the '%s' function must be a 1D, 2D or 3D numpy array. You inserted a %iD array." % (funcName,noDims) )
    boxLength = NumpyArray(noDims,boxLength,funcName,'boxLength')
    shape = Data.shape
    if bin_type not in ['linear','logarithm']: throwError( "The argument 'bin_type' in function '%s' must have one of the values in the list '%s'. You inserted for the argument the value '%s'." % (funcName,str(['linear','logarithm']),bin_type) )
    
    if VERBOSE: print "Computing the FT spectrum of a %iD input data of shape %s ..." % (noDims,str(shape))
    dxValue_lin, indexValue_lin, kValue_lin = "(kMax-kMin)/noBins", "int( floor( (tempK-kMin)/dx ) )", "kMin + (i+0.5)*dx"
    dxValue_log, indexValue_log, kValue_log = "log10(kMax/kMin)/noBins", "int( floor( log10(tempK/kMin) /dx ) )", "kMin * pow( 10, (i+0.5)*dx )"
    dxValue, indexValue, kValue = None, None, None
    if bin_type is 'linear': dxValue, indexValue, kValue = dxValue_lin, indexValue_lin, kValue_lin
    elif bin_type is 'logarithm': dxValue, indexValue, kValue = dxValue_log, indexValue_log, kValue_log
    
    
    # C++ code that does the computation
    FourierTransformSpectrumSupportCode = """
    #line 1000 "FourierTransformSpectrumSupportCode"
    #include <cmath>
    typedef float   Real;
    #include <k_space.h>
    using namespace std;
    """
    
    FourierTransformSpectrumCode_1D = """
    #line 1000 "FourierTransformSpectrumCode_1D"
    int const noDims = 1, noBins = int(parameters(0)); int grid[] = {data.extent(0)};
    Real boxLength[] = {parameters(1)};
    K_space<Real,noDims> k(grid, noDims, boxLength);
    Real kMin = k.k(1), kMax = k.k(grid[0]/2);
    blitz::Array<int,noDims> count(noBins); count = 0;
    
    // computes dk for the given bin type
    Real dx = %s;
    
    // loop over all the entries and find the corresponding k-bin where to put the data
    result = 0.;
    for (int i1=0; i1<grid[0]; ++i1)
    {
        Real tempK = k.k( i1 );
        int index = %s;
        if (index>=0 and index<noBins)
        {
            result(index,1) += real( data(i1) );
            result(index,2) += abs( data(i1) );
            ++count(index);
        }
    }
    for (int i=0; i<noBins; ++i)
    {
        result(i,0) = %s;
        if (count(i)==0) continue;
        result(i,1) /= count(i);
        result(i,2) /= count(i);
    }
    """ % (dxValue,indexValue,kValue)
    FourierTransformSpectrumCode_2D = """
    #line 1000 "FourierTransformSpectrumCode_2D"
    int const noDims = 2, noBins = int(parameters(0)); int grid[] = {data.extent(0),data.extent(1)};
    Real boxLength[] = {parameters(1),parameters(2)};
    K_space<Real,noDims> k(grid, noDims, boxLength);
    Real kMin = k.k(1,0), kMax = k.k(grid[0]/2,grid[1]/2);
    blitz::Array<int,noDims> count(noBins); count = 0;
    
    // computes dk for the given bin type
    Real dx = %s;
    
    // loop over all the entries and find the corresponding k-bin where to put the data
    result = 0.;
    for (int i1=0; i1<grid[0]; ++i1)
        for (int i2=0; i2<grid[1]; ++i2)
        {
            Real tempK = k.k( i1,i2 );
            int index = %s;
            if (index>=0 and index<noBins)
            {
                result(index,1) += real( data(i1,i2) );
                result(index,2) += abs( data(i1,i2) );
                ++count(index);
            }
        }
    for (int i=0; i<noBins; ++i)
    {
        result(i,0) = %s;
        if (count(i)==0) continue;
        result(i,1) /= count(i);
        result(i,2) /= count(i);
    }
    """ % (dxValue,indexValue,kValue)
    FourierTransformSpectrumCode_3D = """
    #line 1000 "FourierTransformSpectrumCode_3D"
    int const noDims = 3, noBins = int(parameters(0)); int grid[] = {data.extent(0),data.extent(1),data.extent(2)};
    Real boxLength[] = {parameters(1),parameters(2),parameters(3)};
    K_space<Real,noDims> k(grid, noDims, boxLength);
    Real kMin = k.k(1,0,0), kMax = k.k(grid[0]/2,grid[1]/2,grid[2]/3);
    blitz::Array<int,noDims> count(noBins); count = 0;
    
    // computes dk for the given bin type
    Real dx = %s;
    
    // loop over all the entries and find the corresponding k-bin where to put the data
    result = 0.;
    for (int i1=0; i1<grid[0]; ++i1)
        for (int i2=0; i2<grid[1]; ++i2)
            for (int i3=0; i3<grid[2]; ++i3)
            {
                Real tempK = k.k( i1,i2,i3 );
                int index = %s;
                if (index>=0 and index<noBins)
                {
                    result(index,1) += abs( data(i1,i2,i3) );
                    result(index,2) += real( data(i1,i2,i3) );
                    ++count(index);
                }
            }
    for (int i=0; i<noBins; ++i)
    {
        result(i,0) = %s;
        if (count(i)==0) continue;
        result(i,1) /= count(i);
        result(i,2) /= count(i);
    }
    """ % (dxValue,indexValue,kValue)
    
    
    # Call the C++ code
    parameters = np.zeros( 1+noDims, boxLength.dtype )
    parameters[0], parameters[1:] = noBins, boxLength[:]
    data = fft.fftn(Data)
    result = np.zeros( (noBins,3), np.float64 )
    inline( eval("FourierTransformSpectrumCode_%iD" % noDims), ['data','parameters','result'], type_converters=converters.blitz, support_code=FourierTransformSpectrumSupportCode, include_dirs=ownCIncludePath )
    
    return result
