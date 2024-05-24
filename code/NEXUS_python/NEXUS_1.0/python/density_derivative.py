#!/usr/bin/python

import sys
import numpy as np
import numpy.fft as fft
from scipy.weave import inline, converters
from density import DensityHeader, readDensityData, writeDensityData
from mathTools import KSpace
from miscellaneous import throwError



help = """
This program is used to compute the derivatives of input field (scalar or vector field - determined from the type of the input file). The derivatives are computed in Fourier Space.

Usage: %s  inputDensityFile  outputFile  derivativeType

Where the output is also a density file and the derivative type can be:
\t 1 = compute the first derivative (OUTPUT = a 3-component vector for an input scalar/ a 9-component tensor for an input vector field).
\t 2 = compute the 2nd derivative (OUTPUT = a tensor for an input scalar/ does NOT work for a vector field).
\t gradient = computes the gradient of an input scalar field (3-components).
\t divergence = computes the divergence of a input vector field (1-component).
\t shear = computes the shear of a input vector field (6-components).
\t vorticity = computes the vorticity of a input vector field (3-components).
\n""" % sys.argv[0]
allowedFiles = ['Density']
allowedFields = ['scalar','vector']
dataTypes = { 1:'scalar', 3:'vector', -1:'unknown', 5:'unknown', 9:'unknown'}
derivativeTypes = ['1','2','gradient','divergence','shear','vorticity']



# read the arguments supplied by the user
if len(sys.argv)!=4:
    print help
    sys.exit(1)

inputFile = sys.argv[1]
outputFile = sys.argv[2]
derivativeType = sys.argv[3]
programOptions = sys.argv[0].rsplit('/')[-1] + ' ' + ' '.join(sys.argv[1:])
if derivativeType not in derivativeTypes:
    throwError( "Unrecognized derivative type in the program '%s'. Allowed derivatives types are '%s', while the requested derivative is '%s'." % (sys.argv[0],str(derivativeTypes),derivativeType) )


#define some functions to return the k-space 1st derivative
def derivative_1(fftShape,K,dX):
    momenta = KSpace(fftShape,dX*fftShape)
    
    mx = 1.j/dX[0] * np.sin(K[0]*dX[0])
    my = 1.j/dX[1] * np.sin(K[1]*dX[1])
    mz = 1.j/dX[2] * np.sin(K[2]*dX[2])
    ax, ay, az = np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128')
    for i in range(fftShape[0]):
        ax[i,:,:] = mx[i]
        for j in range(fftShape[1]):
            ay[i,j,:] = my[j]
            az[i,j,:] = mz
    return (ax,ay,az)

def derivative_2(fftShape,K,dX):
    secondDerivativeCode = """
    #line 1000 "firstDerivativeCode"
    int const nx = fftShape(0), ny = fftShape(1), nz = fftShape(2);
    complex<double> mx[nx], my[ny], mz[nz], mxx[nx], myy[ny], mzz[nz];
    for (int i=0; i<nx; ++i)
    {
        mx[i] = complex<double>( 0., -1./dx * std::sin(kx(i)*dx) );
        mxx[i] = complex<double>( -4./(dx*dx) * std::sin(kx(i)*dx/2) * std::sin(kx(i)*dx/2), 0. );
    }
    for (int i=0; i<ny; ++i)
    {
        my[i] = complex<double>( 0., -1./dy * std::sin(ky(i)*dy) );
        myy[i] = complex<double>( -4./(dy*dy) * std::sin(ky(i)*dy/2) * std::sin(ky(i)*dy/2), 0. );
    }
    for (int i=0; i<nz; ++i)
    {
        mz[i] = complex<double>( 0., -1./dz * std::sin(kz(i)*dz) );
        mzz[i] = complex<double>( -4./(dz*dz) * std::sin(kz(i)*dz/2) * std::sin(kz(i)*dz/2), 0. );
    }
    
    for (int i=0; i<nx; ++i)
        for (int j=0; j<ny; ++j)
            for (int k=0; k<nz; ++k)
            {
                axx(i,j,k) = mxx[i];
                axy(i,j,k) = mx[i] * my[j];
                axz(i,j,k) = mx[i] * mz[k];
                ayy(i,j,k) = myy[j];
                ayz(i,j,k) = my[j] * mz[k];
                azz(i,j,k) = mzz[k];
            }
    """
    
    dx, dy, dz = dX
    kx, ky, kz = K
    axx, axy, axz, ayy, ayz, azz = np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128'), np.zeros(fftShape, dtype='complex128')
    inline( secondDerivativeCode, ['fftShape','dx','dy','dz','kx', 'ky', 'kz', 'axx', 'axy', 'axz', 'ayy', 'ayz', 'azz'], type_converters=converters.blitz, headers=['<complex>','<cmath>'] )
    return (axx,axy,axz,ayy,ayz,azz)



# read the density file
headerType = getHeaderType(inputFile)
if headerType not in allowedFiles:
    throwError( "Unrecognized input file type in the program '%s'. Allowed file types are '%s', while the input file is of type: '%s'." % (sys.argv[0],str(allowedFiles),headerType) )

header, data = None, None
if headerType is allowedFiles[0]:
    header, data = readDensityData(inputFile)
data.shape = header.gridSize
length = header.box[1::2] - header.box[0::2] #the box length along each dimensions in Mpc
fftShape = header.gridSize
dataType = dataTypes[header.DataComponents()]
if dataType not in allowedFields:
    throwError( "Unrecognized input field type in the program '%s'. Allowed input field types are '%s', while the input field is of type: '%s'." % (sys.argv[0],str(allowedFields),dataType) )


# Get the Fourier Space momenta
kx = fft.fftfreq( fftShape[0] ) * 2*np.pi*fftShape[0]/length[0] # now kx contains the value of the momentum along the x-direction (kx(n) = 2\pi/L * n)
ky = fft.fftfreq( fftShape[1] ) * 2*np.pi*fftShape[1]/length[1]
kz = fft.fftfreq( fftShape[2] ) * 2*np.pi*fftShape[2]/length[2]
dx, dy, dz = length[0]/fftShape[0], length[1]/fftShape[1], length[2]/fftShape[2]


# computations for scalar fields
result = None
if dataType is 'scalar':
    tempField = 'scalar'
    tempDerivatives = ['1', '2', 'gradient']
    if derivativeType not in tempDerivatives:
        throwError( "Unrecognized derivative type '%s' for %s fields. Allowed derivatives for %s fields are '%s', while the required derivative is '%s'." % (derivativeType,tempField,tempField,str(tempDerivatives),derivativeType) )
    
    if derivativeType in ['1','gradient']:
        print "Computing the first derivative (gradient) of the scalar field:"
        
        print "\t computing the Fourier transform of the data ..."
        FT = fft.fftn(data)
        
        print "\t computing the derivative in Fourier space ..."
        A = derivative_1( fftShape, (kx,ky,kz), (dx,dy,dz) )
        FTx = FT * A[0]
        FTy = FT * A[1]
        FTz = FT * A[2]
        
        print "\t transforming from Fourier Space to real space ..."
        dataX = fft.ifftn(FTx).astype(data.dtype).real
        dataY = fft.ifftn(FTy).astype(data.dtype).real
        dataZ = fft.ifftn(FTz).astype(data.dtype).real
        dataX.shape, dataY.shape, dataZ.shape = -1, -1, -1
        result = np.column_stack( dataX, dataY, dataZ )





