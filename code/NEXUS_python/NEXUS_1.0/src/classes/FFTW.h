
#ifndef FFTW_HEADER
#define FFTW_HEADER

#include <defines.h>
#include <iostream>
using std::cout;


#ifdef DOUBLE
    #include "FFTW_double.h"
    typedef FFTW_double    FFTW_class;
#else
    #include "FFTW_float.h"
    typedef FFTW_float     FFTW_class;
#endif


/* Computes the DFT:
        FFTW_type: FFTW_class::R2C, FFTW_class::C2R, FFTW_class::C2C_x2k, FFTW_class::C2C_k2x. */
template <typename Array1, typename Array2 >
void computeFFTW(Array1 & input,
                 Array2 & output,
                 int const FFTW_type,
                 int const noThreads = 1,
                 bool const VERBOSE = true)
{
    // define some constants
    int const dimension = input.dimensions();
    int gridSize[dimension];
    for (int i=0; i<dimension; ++i)
        gridSize[i] = input.getSize(i);
    if ( sizeof(output[0])==sizeof(Real) )
        for (int i=0; i<dimension; ++i)
            gridSize[i] = output.getSize(i);
    
    if (VERBOSE)
    {
        cout << "Computing a " << dimension << "D FFTW of size " << gridSize[0];
        for (int i=1; i<dimension; ++i)
            cout << " * " << gridSize[i];
        if ( noThreads!=1 ) cout << " using " << noThreads << " threads";
        cout << " ... " << std::flush;
    }
    
    if ( not( FFTW_type==FFTW_class::R2C or FFTW_type==FFTW_class::C2R or FFTW_type==FFTW_class::C2C_x2k or FFTW_type==FFTW_class::C2C_k2x ) )
        throwError( "Unknown value for argument 'FFTW_type' of function 'computeFFTW' in header 'FFTW.h'. This argument gives the type of the FFTW tranform to be computed. Allowed values: FFTW_class::R2C, FFTW_class::C2R, FFTW_class::C2C_x2k and FFTW_class::C2C_k2x." );
    
    if ( sizeof(input[0])==sizeof(Real) and FFTW_type!=FFTW_class::R2C )
        throwError( "The first argument of function 'computeFFTW' is an array of reals so the only allowed FFTW type is FFTW_class::R2C." );
    if ( sizeof(output[0])==sizeof(Real) and FFTW_type!=FFTW_class::C2R )
        throwError( "The second argument of function 'computeFFTW' is an array of reals so the only allowed FFTW type is FFTW_class::C2R." );
    if ( sizeof(input[0])!=2*sizeof(Real) and sizeof(output[0])!=2*sizeof(Real) )
        throwError( "At least one of the two arrays in function 'computeFFTW' must have a complex type." );
    if ( (sizeof(input[0])!=2*sizeof(Real) or sizeof(output[0])!=2*sizeof(Real)) and (FFTW_type==FFTW_class::C2C_x2k or FFTW_type==FFTW_class::C2C_k2x) )
        throwError( "Both arrays supplied to function 'computeFFTW' must have a complex type since you requested for a C2C FFTW transform." );
    
    FFTW_class fftwClass( gridSize, dimension, FFTW_type, noThreads );
    if ( sizeof(output[0])==sizeof(Real) )
        fftwClass.plan( output.ptrData(), input.ptrData() );
    else 
        fftwClass.plan( input.ptrData(), output.ptrData() );
    fftwClass.FFTW();
    
    if (VERBOSE)
        cout << "Done!\n";
}




#endif
