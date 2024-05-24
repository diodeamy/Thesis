#include <fftw3.h>
#include <iostream>
#include <string>
#include <sstream>


#include "FFTW_float.h"
#include <miscellaneous.h>
using namespace std;


/*! These can be used by the user to specify to the class what kind of FT it wants. */
int const FFTW_float::R2C = 1;
int const FFTW_float::C2R = 2;
int const FFTW_float::C2C_x2k = 11;
int const FFTW_float::C2C_k2x = 12;



/* Class constructors. In the following: r2c = real to complex FT, c2r = complex to real FT and c2c = complex to complex FT
Values allowed for typeFFTW:
       R2C (=1) - r2c -- forward FT
       C2R (=2) - c2r -- backward FT 
       C2C_x2k (=11) - complex to complex -- from coord. space to k-space
       C2C_k2x (=12) - complex to complex -- from k-space to coord. space
*/
FFTW_float::FFTW_float(int dimensions[],
                       int const noDimensions,
                       int const typeFFTW,
                       int const noThreads)
{
    this-> setClassVariables( dimensions, noDimensions, typeFFTW, noThreads );
    fftwPlan = NULL;
}
FFTW_float::FFTW_float(int const dim_1,
                       int const typeFFTW,
                       int const noThreads)
{
    int const noDimensions = 1;
    int dimensions[noDimensions] = { dim_1 };
    
    this-> setClassVariables( dimensions, noDimensions, typeFFTW, noThreads );
    fftwPlan = NULL;
}
FFTW_float::FFTW_float(int const dim_1,
                       int const dim_2,
                       int const typeFFTW,
                       int const noThreads)
{
    int const noDimensions = 2;
    int dimensions[noDimensions] = { dim_1, dim_2 };
    
    this-> setClassVariables( dimensions, noDimensions, typeFFTW, noThreads );
    fftwPlan = NULL;
}
FFTW_float::FFTW_float(int const dim_1,
                       int const dim_2,
                       int const dim_3,
                       int const typeFFTW,
                       int const noThreads)
{
    int const noDimensions = 3;
    int dimensions[noDimensions] = { dim_1, dim_2, dim_3 };
    
    this-> setClassVariables( dimensions, noDimensions, typeFFTW, noThreads );
    fftwPlan = NULL;
}
/* Class destructor. */
FFTW_float::~FFTW_float()
{
    fftwf_destroy_plan( fftwPlan );
    fftwf_cleanup_threads();
    delete[] dim;
    
    if ( not userPointers[0] and fftwType==FFTW_float::R2C )
        fftwf_free( realArray );
    else if ( not userPointers[0] )
        fftwf_free( reinterpret_cast<FFTW_COMPLEX*>(complexArray_in) );
    
    if ( not userPointers[1] and fftwType==FFTW_float::C2R )
        fftwf_free( realArray );
    else if ( not userPointers[1] )
        fftwf_free( reinterpret_cast<FFTW_COMPLEX*>(complexArray_out) );
}




//! Initialize the FFTW plans
void FFTW_float::plan()
{
    userPointers[0] = false;   //no FFTW 'in' pointers supplied
    userPointers[1] = false;   //no FFTW 'out' pointers supplied
    this->ptrInitialization();                 // reserve memory for the 'in' and 'out' pointers if not supplied by the user
    this->planInitialization();                // give a value to the FFTW plan
}
void FFTW_float::plan(FFTW_REAL *ptrReal,
                      FFTW_COMPLEX *ptrComplex)
{
    this->ptrAssigment( ptrReal, ptrComplex ); // check the user supplied pointers and assign the address they point to to the internal class pointers
    this->ptrInitialization();                 // reserve memory for the 'in' and 'out' pointers if not supplied by the user
    this->planInitialization();                // give a value to the FFTW plan
}
void FFTW_float::plan(FFTW_COMPLEX *ptrComplex_in,
                      FFTW_COMPLEX *ptrComplex_out)
{
    this->ptrAssigment( ptrComplex_in, ptrComplex_out ); // check the user supplied pointers and assign the address they point to to the internal class pointers
    this->ptrInitialization();                 // reserve memory for the 'in' and 'out' pointers if not supplied by the user
    this->planInitialization();                // give a value to the FFTW plan
}
void FFTW_float::plan(FFTW_REAL *ptrReal)
{
    this->ptrAssigment( ptrReal );             // check the user supplied pointers and assign the address they point to to the internal class pointers
    this->ptrInitialization();                 // reserve memory for the 'in' and 'out' pointers if not supplied by the user
    this->planInitialization();                // give a value to the FFTW plan
}
void FFTW_float::plan(FFTW_COMPLEX *ptrComplex,
                      bool const in)
{
    this->ptrAssigment( ptrComplex, in );      // check the user supplied pointers and assign the address they point to to the internal class pointers
    this->ptrInitialization();                 // reserve memory for the 'in' and 'out' pointers if not supplied by the user
    this->planInitialization();                // give a value to the FFTW plan
}




/* Execute the Fourier transform.
*/
void FFTW_float::FFTW()
{
    if ( not fftwData )
        throwError( "You need to supply input values for the input FFTW data before doing the Fourier transform." );
    
    fftwf_execute( fftwPlan );
    fftwExecuted = true;
}


/* Return the type of the Fourier transform
*/
int FFTW_float::FFTW_type() const
{
    return fftwType;
}


/* Output to standard terminal the type of the Fourier tranform that is being executed.
*/
void FFTW_float::outputFTTW_type() const
{
    ostringstream buffer1, buffer2;
    buffer1 << dim[0];
    for (int i=1; i<noDim; ++i)
        buffer1 << " * " << dim[i];
    
    if ( userPointers[0] and userPointers[1] )
        buffer2 << "Both input and ouput pointers are supplied by the main program.";
    else if ( userPointers[0] )
        buffer2 << "Only input pointer is supplied by the main program.";
    else if ( userPointers[1] )
        buffer2 << "Only output pointer is supplied by the main program.";
    else
        buffer2 << "No pointers were supplied by the main program.";
    
    switch (fftwType)
    {
        case FFTW_float::R2C:
            cout << "\nReal to complex Fourier transform of dimension  " << buffer1.str() << ". " << buffer2.str() << "\n";
            break;
        case FFTW_float::C2R:
            cout << "\nComplex to real Fourier transform of dimension  " << buffer1.str() << ". " << buffer2.str() << "\n";
            break;
        case FFTW_float::C2C_x2k:
            cout << "\nComplex to complex coordinate space to k-space Fourier transform of dimension  " << buffer1.str() << ". " << buffer2.str() << "\n";
            break;
        case FFTW_float::C2C_k2x:
            cout << "\nComplex to complex k-space to coordinate space Fourier transform of dimension  " << buffer1.str() << ". " << buffer2.str() << "\n";
            break;
    }
}




//! Pointers to the different input/output arrays
/* Pointer to input complex array. */
FFTW_COMPLEX * FFTW_float::ptrComplexInput()
{
    if ( not this->complexInput() )
        throwError( "Trying to access the input complex pointer in function 'FFTW_float::ptrComplexInput' but the FFTW wasn't initialized for a C2R (real to complex) or C2C Fourier transform. There is no FFTW complex input pointer to be accessed." );
    if ( userPointers[0] )
        throwWarning( "Trying to access an input pointer supplied by the user. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    fftwData = true;  // must be set to true since now the user the insert FFTW input data outside the class
    return complexArray_in;
}


/* Pointer to output complex array. */
FFTW_COMPLEX * FFTW_float::ptrComplexOutput() const
{
    if ( not this->complexOutput() )
        throwError( "Trying to access the output complex pointer in function 'FFTW_float::ptrComplexOutput' but the FFTW wasn't initialized for a R2C (complex to real) or C2C Fourier transform. There is no FFTW complex output pointer to be accessed." );
    if ( userPointers[1] )
        throwWarning( "Trying to access an output pointer supplied by the user. This operation may be valid, but it is suspicious since the user supplied a pointer for the output data for the FFTW.");
    
    return complexArray_out;
}


/* Pointer to input complex array. */
FFTW_REAL * FFTW_float::ptrRealInput()
{
    if ( not this->realInput() )
        throwError( "Trying to access the input real pointer in function 'FFTW_float::ptrRealInput' but the FFTW wasn't initialized for a R2C (real to complex) Fourier transform. There is no FFTW real input pointer to be accessed." );
    if ( userPointers[0] )
        throwWarning( "Trying to access an input pointer supplied by the user. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    fftwData = true;  // must be set to true since now the user the insert FFTW input data outside the class
    return realArray;
}


/* Pointer to output complex array. */
FFTW_REAL * FFTW_float::ptrRealOutput() const
{
    if ( not this->realOutput() )
        throwError( "Trying to access the output real pointer in function 'FFTW_float::ptrRealOutput' but the FFTW wasn't initialized for a C2R (complex to real) Fourier transform. There is no FFTW real output pointer to be accessed." );
    if ( userPointers[1] )
        throwWarning( "Trying to access an output pointer supplied by the user. This operation may be valid, but it is suspicious since the user supplied a pointer for the output data for the FFTW.");
    
    return realArray;
}






//! Copy input
/* Copy the values for the FFTW input array. It copies to the input array of the FFTW the entries at 'location[0]',  'location[offset]', 'location[2*offset]' ...
*/
void FFTW_float::copyInput(FFTW_REAL *location,
                           int const offset)
{
    if ( not this->realInput() )
       throwError( "Trying to copy a real input in function 'FFTW_float::copyInput' but the FFTW wasn't initialized for a R2C (real to complex) Fourier transform. There is no FFTW real input array where to copy the input to function 'FFTW_float::copyInput'." );
    if ( userPointers[0] )
       throwWarning( "Trying to write input values into the user supplied input FFTW array. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    for (int i=0; i<realSize; ++i)
        realArray[i] = location[i*offset];
    
    fftwData = true; //update that indeed there is valid input data for the FFTW
}


/* Copy the values for the FFTW input array. It copies to the input array of the FFTW the entries at 'location[0]',  'location[offset]', 'location[2*offset]' ...
*/
void FFTW_float::copyInput(FFTW_COMPLEX *location,
                           int const offset)
{
    if ( not this->complexInput() )
       throwError( "Trying to copy a complex input in function 'FFTW_float::copyInput' but the FFTW wasn't initialized for a C2R (real to complex) or C2C Fourier transform. There is no FFTW complex input array where to copy the input to function 'FFTW_float::copyInput'." );
    if ( userPointers[0] )
       throwWarning( "Trying to write input values into the user supplied input FFTW array. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    for (int i=0; i<complexSize; ++i)
        complexArray_in[i] = location[i*offset];
    
    fftwData = true; //update that indeed there is valid input data for the FFTW
}


/* Copy the values for the FFTW input array.
*/
void FFTW_float::copyInput(FFTW_REAL *location)
{
    if ( not this->realInput() )
       throwError( "Trying to copy a real input in function 'FFTW_float::copyInput' but the FFTW wasn't initialized for a R2C (real to complex) Fourier transform. There is no FFTW real input array where to copy the input to function 'FFTW_float::copyInput'." );
    if ( userPointers[0] )
       throwWarning( "Trying to write input values into the user supplied input FFTW array. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    for (int i=0; i<realSize; ++i)
        realArray[i] = location[i];
    
    fftwData = true; //update that indeed there is valid input data for the FFTW
}


/* Copy the values for the FFTW input array. It copies to the input array of the FFTW the entries at 'location[0]',  'location[offset]', 'location[2*offset]' ...
*/
void FFTW_float::copyInput(FFTW_COMPLEX *location)
{
    if ( not this->complexInput() )
       throwError( "Trying to copy a complex input in function 'FFTW_float::copyInput' but the FFTW wasn't initialized for a C2R (real to complex) or C2C Fourier transform. There is no FFTW complex input array where to copy the input to function 'FFTW_float::copyInput'." );
    if ( userPointers[0] )
       throwWarning( "Trying to write input values into the user supplied input FFTW array. This operation may be valid, but it is suspicious since the user supplied a pointer for the input data for the FFTW.");
    
    for (int i=0; i<complexSize; ++i)
        complexArray_in[i] = location[i];
    
    fftwData = true; //update that indeed there is valid input data for the FFTW
}




//! Copy output
/* Copy the output to a new array. It copies the output entries at 'location[0]',  'location[offset]', 'location[2*offset]' ...
*/
void FFTW_float::copyOutput(FFTW_REAL *location,
                            int const offset) const
{
    if ( not this->realOutput() )
       throwError( "Trying to copy the real output of the FFTW in function 'FFTW_float::copyOutput' but the FFTW wasn't initialized for a C2R (complex to real) Fourier transform. There is no real output that can be copied by this function." );
    if ( not fftwExecuted and not userPointers[1] )
       throwError( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). There are no results in the output pointer of FFTW." );
    if ( not fftwExecuted and userPointers[1] )
       throwWarning( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). Since the output pointer was supplied by the user, the operation of copying the data may be a valid one and hence will be executed.");
    
    for (int i=0; i<realSize; ++i)
        location[i*offset] = realArray[i];
}


/* Copy the output to a new array. It copies the output entries at 'location[0]',  'location[offset]', 'location[2*offset]' ...
*/
void FFTW_float::copyOutput(FFTW_COMPLEX *location,
                            int const offset) const
{
    if ( not this->complexOutput() )
       throwError( "Trying to copy the real output of the FFTW in function 'FFTW_float::copyOutput' but the FFTW wasn't initialized for a R2C (complex to real) or C2C Fourier transform. There is no complex output that can be copied by this function." );
    if ( not fftwExecuted and not userPointers[1] )
       throwError( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). There are no results in the output pointer of FFTW." );
    if ( not fftwExecuted and userPointers[1] )
       throwWarning( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). Since the output pointer was supplied by the user, the operation of copying the data may be a valid one and hence will be executed.");
    
    for (int i=0; i<complexSize; ++i)
        location[i*offset] = complexArray_out[i];
}



/* Copy the output to a new array.
*/
void FFTW_float::copyOutput(FFTW_REAL *location) const
{
    if ( not this->realOutput() )
       throwError( "Trying to copy the real output of the FFTW in function 'FFTW_float::copyOutput' but the FFTW wasn't initialized for a C2R (complex to real) Fourier transform. There is no real output that can be copied by this function." );
    if ( not fftwExecuted and not userPointers[1] )
       throwError( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). There are no results in the output pointer of FFTW." );
    if ( not fftwExecuted and userPointers[1] )
       throwWarning( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). Since the output pointer was supplied by the user, the operation of copying the data may be a valid one and hence will be executed.");
    
    for (int i=0; i<realSize; ++i)
        location[i] = realArray[i];
}


/* Copy the output to a new array.
*/
void FFTW_float::copyOutput(FFTW_COMPLEX *location) const
{
    if ( not this->complexOutput() )
       throwError( "Trying to copy the real output of the FFTW in function 'FFTW_float::copyOutput' but the FFTW wasn't initialized for a R2C (complex to real) or C2C Fourier transform. There is no complex output that can be copied by this function." );
    if ( not fftwExecuted and not userPointers[1] )
       throwError( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). There are no results in the output pointer of FFTW." );
    if ( not fftwExecuted and userPointers[1] )
       throwWarning( "Trying to copy the results of the FFTW before doing the actual Fourier transform (in function 'FFTW_float::copyOutput'). Since the output pointer was supplied by the user, the operation of copying the data may be a valid one and hence will be executed.");
    
    for (int i=0; i<complexSize; ++i)
        location[i] = complexArray_out[i];
}



//! Check type of input and output
/* Checks if the FFTW input is a real array. */
bool FFTW_float::realInput() const
{
     if ( fftwType==FFTW_float::R2C )
        return true;
     else
         return false;
}

/* Checks if the FFTW input is a complex array. */
bool FFTW_float::complexInput() const
{
     if ( fftwType==FFTW_float::C2R or fftwType==FFTW_float::C2C_x2k or fftwType==FFTW_float::C2C_k2x )
        return true;
     else
         return false;
}

/* Checks if the FFTW output is a real array. */
bool FFTW_float::realOutput() const
{
     if ( fftwType==FFTW_float::C2R )
        return true;
     else
         return false;
}

/* Checks if the FFTW output is a complex array. */
bool FFTW_float::complexOutput() const
{
     if ( fftwType==FFTW_float::R2C or fftwType==FFTW_float::C2C_x2k or fftwType==FFTW_float::C2C_k2x )
        return true;
     else
         return false;
}





//! Private functions
/* Assign values to the fftwPlan for a R2C transform. */
void FFTW_float::initialize_R2C_plan()
{
    fftwPlan = fftwf_plan_dft_r2c( noDim, dim, realArray, reinterpret_cast<fftwf_complex*>(complexArray_out), FFTW_ESTIMATE);
}

/* Assign values to the fftwPlan for a C2R transform. */
void FFTW_float::initialize_C2R_plan()
{
    fftwPlan = fftwf_plan_dft_c2r( noDim, dim, reinterpret_cast<fftwf_complex*>(complexArray_in), realArray, FFTW_ESTIMATE);
}

/* Assign values to the fftwPlan for a C2R transform. */
void FFTW_float::initialize_C2C_plan()
{
    if ( C2C_x2k )
         fftwPlan = fftwf_plan_dft( noDim, dim, reinterpret_cast<fftwf_complex*>(complexArray_in), reinterpret_cast<fftwf_complex*>(complexArray_out), FFTW_BACKWARD, FFTW_ESTIMATE);
    else if ( C2C_k2x )
         fftwPlan = fftwf_plan_dft( noDim, dim, reinterpret_cast<fftwf_complex*>(complexArray_in), reinterpret_cast<fftwf_complex*>(complexArray_out), FFTW_FORWARD, FFTW_ESTIMATE);
}

/* Find the size of the arrays needed for the FFTW transform. */
void FFTW_float::findSize()
{
    // store in some temporary variables the size of the memory that needs to be allocated
    int r2cRealSize = 1;       // size of the real array in r2c and c2r FT
    int r2cComplexSize = 1;    // size of the complex array in r2c and c2r FT
    int c2cSize = 1;           // size of the complex arrays in c2c FT
    for (int i=0; i<noDim-1; ++i)
        r2cRealSize *= dim[i];
    r2cComplexSize = r2cRealSize * (dim[noDim-1]/2+1);
    r2cRealSize *= dim[noDim-1];
    c2cSize = r2cRealSize;
    
    if ( fftwType==R2C or fftwType==C2R )
    {
         realSize = r2cRealSize;
         complexSize = r2cComplexSize;
    }
    else if ( fftwType==C2C_x2k or fftwType==C2C_k2x )
    {
         realSize = 0;
         complexSize = c2cSize;
    }
    else throwError( "Invalid value for variable 'fftwType' in function 'FFTW_float::findSize'." );
}

/* Assign memory and values to the *dim pointer of the array. */
void FFTW_float::initializeDimensions(int dimensions[], int const noDimensions)
{
    lowerBoundCheck( noDimensions, 1, "The argument 'noDimensions' in function 'FFTW_float::initializeDimensions'" );
    noDim = noDimensions;
    
    dim = new (nothrow) int[noDim];
    memoryAllocationCheck( dim, "Error while allocating memory for the 'dim' variable in the function 'FFTW_float::initializeDimensions'." );
    for (int i=0; i<noDim; ++i)
        dim[i] = dimensions[i];
}


/* Reserve memory for the FFTW arrays. Call this function only after first calling 'FFTW_float::findSize'.*/
void FFTW_float::ptrInitialization()
{
    if ( not userPointers[0] )      //reserve memory for input array if not supplied by the user
    {
        if ( fftwType==R2C )
        {
            realArray = reinterpret_cast<FFTW_REAL*>( fftwf_malloc( realSize * sizeof(FFTW_REAL)) );
            memoryAllocationCheck( realArray, "Error while allocating memory for the 'realArray' variable in the constructor of class 'FFTW_float'." );
        }
        else if ( fftwType==C2R or fftwType==C2C_x2k or fftwType==C2C_k2x )
        {
            complexArray_in = reinterpret_cast<FFTW_COMPLEX*>( fftwf_malloc( complexSize * sizeof(FFTW_COMPLEX)) );
            memoryAllocationCheck( complexArray_in, "Error while allocating memory for the 'complexArray_in' variable in the constructor of class 'FFTW_float'." );
        }
        else throwError( "Invalid value for variable 'fftwType' in function 'FFTW_float::ptrInitialization'." );
    }
    
    if ( not userPointers[1] )      // reserve memory for the output array if not supplied by the user
    {
        if ( fftwType==C2R )
        {
            realArray = reinterpret_cast<FFTW_REAL*>( fftwf_malloc( realSize * sizeof(FFTW_REAL)) );
            memoryAllocationCheck( realArray, "Error while allocating memory for the 'realArray' variable in the constructor of class 'FFTW_float'." );
        }
        else if ( fftwType==R2C or fftwType==C2C_x2k or fftwType==C2C_k2x )
        {
            complexArray_out = reinterpret_cast<FFTW_COMPLEX*>( fftwf_malloc( complexSize * sizeof(FFTW_COMPLEX)) );
            memoryAllocationCheck( complexArray_out, "Error while allocating memory for the 'complexArray_out' variable in the constructor of class 'FFTW_float'." );
        }
        else throwError( "Invalid value for variable 'fftwType' in function 'FFTW_float::ptrInitialization'." );
    }
}

/* Set the value of 'fftwPlan'. */
void FFTW_float::planInitialization()
{
    if ( userPointers[0] )
        fftwData = true;
    else fftwData = false;
    
    switch( fftwType )
    {
        case FFTW_float::R2C:
            this->initialize_R2C_plan();
            break;   
        case FFTW_float::C2R:
            this->initialize_C2R_plan();
            break;
        case FFTW_float::C2C_x2k:
            this->initialize_C2C_plan();
            break;
        case FFTW_float::C2C_k2x:
            this->initialize_C2C_plan();
            break;
        default:
             throwError( "The class 'FFTW_float' allows for only 4 different values for the argument 'typeFFTW': 'FFTW_float::R2C', 'FFTW_float::C2R', 'FFTW_float::C2C_x2k' and 'FFTW_float::C2C_k2x'. You supplied a nonrecognized value for 'typeFFTW'." );
    }
}



/* General part of the constructor - deals with the general tasks found in all the class constructors. 
It sets the variables, reserves memory if needed and defines the plans. */
void FFTW_float::setClassVariables(int dimensions[],
                                   int const noDimensions,
                                   int const typeFFTW,
                                   int const noThreads)
{
    fftwType = typeFFTW;
    if ( not (fftwType==R2C or fftwType==C2R or fftwType==C2C_x2k or fftwType==C2C_k2x) )
        throwError( "The class 'FFTW_float' allows for only 4 different values for the argument 'typeFFTW': 'FFTW_float::R2C', 'FFTW_float::C2R', 'FFTW_float::C2C_x2k' and 'FFTW_float::C2C_k2x'. You supplied a nonrecognized value for 'typeFFTW'." );
    if ( noDimensions<=0 ) throwError( "Cannot compute the Fourier transform with a dimension smaller than 1." );
    
    userPointers[0] = false;   //no FFTW 'in' pointers supplied
    userPointers[1] = false;   //no FFTW 'out' pointers supplied
    fftwExecuted = false;      // the FFTW was not executed
    _noThreads = noThreads;
    
    this->initializeDimensions( dimensions, noDimensions );   // read the number of dimensions and their size
    this->findSize();          // store in some temporary variables the size of the memory that needs to be allocated
    
    
    // set FFTW characteristics
    fftwf_init_threads();
    fftwf_plan_with_nthreads( noThreads );
}


/* Functions that assigns the user supplied pointers to internal FFTW pointers. */
void FFTW_float::ptrAssigment(FFTW_REAL *ptrReal,
                              FFTW_COMPLEX *ptrComplex)
{
    realArray = ptrReal;
    if ( fftwType==R2C )
        complexArray_out = ptrComplex;
    else if ( fftwType==C2R ) 
        complexArray_in = ptrComplex;
    else throwError( "When supplying a real and a complex pointer the FFTW type can be only R2C or C2R (error in function 'FFTW_float::ptrAssigment')." );
    
    userPointers[0] = true; // FFTW 'in' pointers supplied
    userPointers[1] = true; // FFTW 'out' pointers supplied
}
void FFTW_float::ptrAssigment(FFTW_REAL *ptrReal)
{
    realArray = ptrReal;
    if ( fftwType==R2C )
    {
        userPointers[0] = true; // FFTW 'in' pointers supplied
        userPointers[1] = false; // no FFTW 'out' pointers supplied
    }
    else if ( fftwType==C2R ) 
    {
        userPointers[0] = false; // FFTW 'in' pointers supplied
        userPointers[1] = true; // no FFTW 'out' pointers supplied
    }
    else throwError( "When supplying a real pointer the FFTW type can be only R2C or C2R (error in function 'FFTW_float::ptrAssigment')." );
}
void FFTW_float::ptrAssigment(FFTW_COMPLEX *ptrComplex,
                              bool in)
{
    bool complexBool = false;
    if ( fftwType==C2C_x2k or fftwType==C2C_k2x )
        complexBool = true;
    
    if ( fftwType==R2C or (complexBool and not in) )
    {
        complexArray_out = ptrComplex;
        userPointers[0] = false; // FFTW 'in' pointers supplied
        userPointers[1] = true; // no FFTW 'out' pointers supplied
    }
    else if ( fftwType==C2R or (complexBool and in) )
    {
        complexArray_in = ptrComplex;
        userPointers[0] = true; // FFTW 'in' pointers supplied
        userPointers[1] = false; // no FFTW 'out' pointers supplied
    }
    else throwError( "When supplying a complex pointer the FFTW type can be only R2C, C2R or C2C (error in function 'FFTW_float::ptrAssigment')." );
}
void FFTW_float::ptrAssigment(FFTW_COMPLEX *ptrComplex_in,
                              FFTW_COMPLEX *ptrComplex_out)
{
    complexArray_in = ptrComplex_in;
    complexArray_out = ptrComplex_out;
    
    userPointers[0] = true; // FFTW 'in' pointers supplied
    userPointers[1] = true; // FFTW 'out' pointers supplied
    
    if ( not (fftwType==C2C_x2k or fftwType==C2C_k2x) )
        throwError( "When supplying two complex pointers the FFTW type can be only C2C_x2k or C2C_k2x (error in function 'FFTW_float::ptrAssigment')." );
}




