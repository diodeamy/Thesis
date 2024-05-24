#include <iostream>
#include <cmath>
#include <omp.h>

#include "additionalComputations.h"
#include <array.h>
#include <k_space.h>
#include <FFTW.h>
#include <MMF_file.h>
using namespace std;

typedef K_space<Real,NO_DIM>    Momenum3D;




/* This function takes the logarithm of the density field.

"minValue" sets the minimal value of the density for which to take the logarithm - all density values with values < 'minValue' are taken as = 'minValue'. */
void densityLogarithm(ArrayReal3D *result,
                      Real const minValue,
                      bool VERBOSE)
{
    if ( not minValue>Real(0.) ) throwError( "The argument 'minValue' of function 'densityLogarithm' must be always greater than 0." );
    Real minLog = log10( minValue );
    
    if (VERBOSE) cout << "Computing the logarithm of the density field ... ";
    for (size_t i=0; i<result->size(); ++i)
    {
        if ( (*result)[i]>minValue )
            (*result)[i] = Real( log10( (*result)[i] ) );
        else
            (*result)[i] = minLog;
    }
    if (VERBOSE) cout << "Done.\n";
}

void densityLogarithm(ArrayReal3D &input,
                      Real const minValue,
                      ArrayReal3D *result,
                      bool VERBOSE)
{
    if ( not minValue>Real(0.) ) throwError( "The argument 'minValue' of function 'densityLogarithm' must be always greater than 0." );
    Real minLog = log10( minValue );
    
    if (VERBOSE) cout << "Computing the logarithm of the density field ... ";
    for (size_t i=0; i<input.size(); ++i)
    {
        if ( input[i]>minValue )
            (*result)[i] = Real( log10( input[i] ) );
        else
            (*result)[i] = minLog;
    }
    if (VERBOSE) cout << "Done.\n";
}
/* Takes the exponential of the logarithm of the density. */
void densityExponential(ArrayReal3D *result,
                        Real const factor,
                        bool VERBOSE)
{
    if (VERBOSE) cout << "Computing the exponential of the logarithm of the density field ... ";
    for (size_t i=0; i<result->size(); ++i)
        (*result)[i] = pow( 10., (*result)[i]*factor );
    if (VERBOSE) cout << "Done.\n";
}



/* Computes the Fourier Transform of the gravitational field from the FT of the density field by multiplication with 1/k^2 in k-space. */
void computeFT_gravitationalField(ArrayComplex3D *result,
                                  size_t *grid,
                                  Real const boxLength[],
                                  Real const factor)
{
    Real kMin = Real(.1);    // discard the large scale modes
    kMin = kMin * kMin;
    intervalCheck( NO_DIM, 2, 3, "'NO_DIM' in function 'computeFT_gravitationalField'" );
    
    Array<Int,1> gridSize( NO_DIM );
    gridSize = result->getSize<Int>();
    Momenum3D k( grid, NO_DIM, boxLength );
    Real multFactor = Real(-1. * factor);
#if NO_DIM==2
    for (Int i1=0; i1<gridSize[0]; ++i1)
        for (Int i2=0; i2<gridSize[1]; ++i2)
        {
            (*result)(i1,i2) *= 1./k.kSquare( i1, i2 ) * multFactor;
            if ( k.nSquare(i1,i2)<=kMin ) (*result)(i1,i2) = Real(0.);
        }
#elif NO_DIM==3
    for (Int i1=0; i1<gridSize[0]; ++i1)
        for (Int i2=0; i2<gridSize[1]; ++i2)
            for (Int i3=0; i3<gridSize[2]; ++i3)
            {
                (*result)(i1,i2,i3) *= 1./k.kSquare( i1, i2, i3 ) * multFactor;
                if ( k.nSquare(i1,i2,i3)<=kMin ) (*result)(i1,i2,i3) = Real(0.);
            }
#endif
}


/* This function solves the Poisson equation. 

'factor' gives the normalization factor - see function 'computeFT_gravitationalField' for more details. */
void solvePoisson(ArrayReal3D *result,
                  Real boxLength[],
                  Real const factor,
                  string variableName)
{
    cout << "Computing the " << variableName << " using the Poisson equation ... " << flush;
    Array<size_t,1> grid = result->axisSize<size_t>();
    Real totalSize = Real( result->size() );
    
    // first due a FT
    ArrayComplex3D denComplex( grid[0], grid[1], grid[2]/2+1 );
    computeFFTW( *result, denComplex, FFTW_class::R2C, omp_get_max_threads() );
    
    // solve Poissoin in Fourier Space - the 1/totalSize factor comes from the normalization of the FT
    computeFT_gravitationalField( &denComplex, grid.ptrData(), boxLength, factor / totalSize );
    
    // translate back to real space via another FT
    computeFFTW( denComplex, *result, FFTW_class::C2R, omp_get_max_threads() );
    
    cout << "Done.\n";
}



/* This function takes the logarithm of the absolute value of the divergence. */
void velocityDivergenceLogarithm(ArrayReal3D *result,
                                 Real const minValue)
{
    if ( not minValue>Real(0.) ) throwError( "The argument 'minValue' of function 'velocityDivergenceLogarithm' must be always greater than 0." );
    Real minLog = log10( minValue );
    
    cout << "Computing the logarithm of the absolute value of the velocity divergence field ... ";
    for (size_t i=0; i<result->size(); ++i)
    {
        Real temp = fabs( (*result)[i] );
        if ( temp>minValue )
            (*result)[i] = Real( log10( temp ) );
        else
            (*result)[i] = minLog;
    }
    cout << "Done.\n";
}

/* This function takes the logarithm of the positive and negative values of the divergence and shifts the two parts by a certain 'offset'. */
void velocityDivergenceLogarithm(ArrayReal3D *result,
                                 Real const minValue,
                                 Real const hubbleValue)
{
    if ( not minValue>Real(0.) ) throwError( "The argument 'minValue' of function 'velocityDivergenceLogarithm' must be always greater than 0." );
    Real minLog = log10( minValue );
    
    cout << "Computing the logarithm of the modified velocity divergence field ... ";
    for (size_t i=0; i<result->size(); ++i)
    {
        if ( (*result)[i]<Real(0.) ) (*result)[i] *= Real(-1.);
        else if ( (*result)[i]>=Real(0.) and (*result)[i]<hubbleValue ) (*result)[i] *= Real(-1.);
        
        (*result)[i] += hubbleValue;
        
        Real temp = fabs( (*result)[i] );
        if ( temp>minValue )
            (*result)[i] = Real( log10( temp ) );
        else
            (*result)[i] = minLog;
    }
    cout << "Done.\n";
}



/* This function reads a clean response file and uses this result to mask the input field. 

'maskValue' is the value with which the field is replaced at the mask pixels. */
void maskInputField(ArrayReal3D *result,
                    MMF_header &mmfHeader,
                    Real maskValue,
                    string maskFilename,
                    string maskName)
{
    cout << "\nReading the " << maskName << " clean MMF response:" << flush;
    MMF_header tempHeader;
    readMMF_Header( &tempHeader, maskFilename );
    mmfHeader.compatible( tempHeader );// check that the two MMF files are compatible (i.e. have the same grid dimensions)
    if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE )
        throwError( "The input file '" + maskFilename + "' is not a 'clean' MMF response file." );

    Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
    readMMF_cleanResponse( &tempResponse, maskFilename );
        
    //discard the values of 'result' in the pixels where 'tempResponse' is a valid mask
    cout << "Masking the input field using the " << maskName << " clean MMF response ... " << flush;
    if ( result->size()!=tempResponse.size() )
        throwError( "You are trying to do operations with array of different sizes. Error in function 'maskInputField'." );
    
    for (size_t i=0; i<result->size(); ++i)
        if ( tempResponse[i]>=1 )
            (*result)[i] = maskValue;
    cout << "Done.\n";
}




