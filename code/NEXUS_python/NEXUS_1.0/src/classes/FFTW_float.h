/*

Class to implement the FFTW Fourier transform in an easily accesible form.
To use the class the following steps must be taken:
    1) Call class constructor: FFTW_float A(...);
    2) Set the plan:           A.plan(... can supply pointers to in/output arrays...)
    3) If didn't supply pointers, set the input data for the FFTW - the FFTW will not execute otherwise.
    4) Execute the FFTW:       A.FFTW();
    5) Access the resulting data using the out pointer supplied or through the class functions.
    6) Can execute the FFTW multiple times.

*/

#ifndef CLASS_FFTW_FLOAT 
#define CLASS_FFTW_FLOAT

#include <complex>
#include <fftw3.h>


typedef float                   FFTW_REAL;
typedef std::complex<float>     FFTW_COMPLEX;



class FFTW_float
{
    bool fftwData;      // true if the arrays have been initialized with data, or the plans are initialized with arrays defined outside the class 
    int fftwType;       // the type of the FT (r2c, c2r, c2c forward or backward)
    bool userPointers[2];  // true if the user supplied FFTW 'in' or/and 'out' pointers to store the data - used for the destructor
    bool fftwExecuted;  // true after calling function fftw();
    int noDim, *dim;
    int realSize, complexSize;     //keep track of the number of elements in the real and compplex arrays respectively
    int _noThreads;     // number of threads
    fftwf_plan fftwPlan;
    
    
    public:
    static int const R2C, C2R, C2C_x2k, C2C_k2x;
    FFTW_COMPLEX *complexArray_in, *complexArray_out;
    FFTW_REAL *realArray;
    
    
    // class constructors
    FFTW_float(int dimensions[],
               int const noDimensions,
               int const typeFFTW,
               int const noThreads = 1);
    FFTW_float(int const dim_1,
               int const typeFFTW,
               int const noThreads = 1);
    FFTW_float(int const dim_1,
               int const dim_2,
               int const typeFFTW,
               int const noThreads = 1);
    FFTW_float(int const dim_1,
               int const dim_2,
               int const dim_3,
               int const typeFFTW,
               int const noThreads = 1);
    
    ~FFTW_float();
    
    void plan();                      // construct FFTW plan using internal class pointers  
    void plan(FFTW_REAL *ptrReal,
              FFTW_COMPLEX *ptrComplex); // construct plan using user supplied pointers
    void plan(FFTW_COMPLEX *ptrComplex_in,
              FFTW_COMPLEX *ptrComplex_out);
    void plan(FFTW_REAL *ptrReal);
    void plan(FFTW_COMPLEX *ptrComplex,
              bool const in);         // in tells the function if the complex pointer is for the in ('in'=true) or out ('in'=false) data pointer  
    
    void FFTW();                      // execute plan (forward or backward depending on value of fftwType)
    int FFTW_type() const;            // write to stdin the type and dimension of the FFTW
    void outputFTTW_type() const;     // write to stdin the type and dimension of the FFTW
    
    FFTW_COMPLEX * ptrComplexInput();             // return a pointer to the complex input array
    FFTW_COMPLEX * ptrComplexOutput() const;      // return a pointer to the complex output array
    FFTW_REAL * ptrRealInput();                   // return a pointer to the float input array
    FFTW_REAL * ptrRealOutput() const;            // return a pointer to the float output array
    
    void copyInput(FFTW_REAL *location, int const offset);      // copy real input from 'location' pointer at every 'offset' memory location starting with first
    void copyInput(FFTW_COMPLEX *location, int const offset);   // copy complex input from 'location' pointer at every 'offset' memory location starting with first
    void copyInput(FFTW_REAL *location);                        // copy real input from 'location' pointer
    void copyInput(FFTW_COMPLEX *location);                     // copy complex input from 'location' pointer
    
    void copyOutput(FFTW_REAL *location, int const offset) const;   // copy real output to 'location' pointer at every 'offset' memory location starting with first
    void copyOutput(FFTW_COMPLEX *location, int const offset) const;// copy complex output to 'location' pointer at every 'offset' memory location starting with first
    void copyOutput(FFTW_REAL *location) const;                     // copy real output to 'location' pointer
    void copyOutput(FFTW_COMPLEX *location) const;                  // copy complex output to 'location' pointer
    
    bool realInput() const;      // check if the input array is real
    bool complexInput() const;   // check if the input array is complex
    bool realOutput() const;     // check if the output array is real
    bool complexOutput() const;  // check if the output array is complex
    
    
    
    private:
    void initialize_R2C_plan();  // initialize the respective plans
    void initialize_C2R_plan();
    void initialize_C2C_plan();
    
    void findSize();    // computes the size of the array needed for the FFTW
    void initializeDimensions(int dimensions[], int const noDimensions);// reserve memory for the *dim pointer and assign values
    void ptrInitialization();                     // reserve memory for FFTW variables
    void planInitialization();                    // set the value of the FFTW plan
    void setClassVariables(int dimensions[],
                           int const noDimensions,
                           int const typeFFTW,
                           int const noThreads);  // set all the variables of the FFTW class
    
    // assign the user supplied pointers to the internal ones
    void ptrAssigment(FFTW_REAL *ptrReal,
                      FFTW_COMPLEX *ptrComplex);
    void ptrAssigment(FFTW_REAL *ptrReal);
    void ptrAssigment(FFTW_COMPLEX *ptrComplex,
                      bool in); //variable 'in' says if this a 'in/out' complex pointer for a C2C transform
    void ptrAssigment(FFTW_COMPLEX *ptrComplex_in,
                      FFTW_COMPLEX *ptrComplex_out);
};


#endif
