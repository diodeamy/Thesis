#ifndef DEFINES_HEADER
#define DEFINES_HEADER

#include <cstddef>




// number of OPENMP threads
#define ENABLE_OPENMP
#ifndef NO_OPENMP_THREADS
#define NO_OPENMP_THREADS 1
#endif

// define conversion factor to get Gadget length units in Mpc
#ifndef MPC_UNIT
#define MPC_UNIT 1000.
#endif


// choose number of dimensions; not all functions work for 2D
#ifdef _2D
    #define NO_DIM 2
    #define NO_DIM2 4
#else
    #define NO_DIM 3
    #define NO_DIM2 9
#endif


// define some shorthand notations
typedef short int                      shortInt;
typedef size_t                         Int;
#ifdef DOUBLE
    typedef double                         Real;
#else
    typedef float                          Real;
#endif



// physical constants
// the critical density in units of (M0/h) / (Mpc/h)^3
#define RHO_CRITICAL     27.7538e+10    



#endif
