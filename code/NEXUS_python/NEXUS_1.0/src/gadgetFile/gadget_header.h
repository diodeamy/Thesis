#ifndef GADGET_HEADER
#define GADGET_HEADER

#include <algorithm>




//see Gadget2 documentation for meaning of each variable
struct Gadget_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    int      flag_stellarage;
    int      flag_metals;
    int      num_total_particles_hw[6]; /* High word of total particle number */
    int      flag_entropy_instead_u;    /*< flags that IC-file contains entropy instead of u */
    int      flag_doubleprecision;      /*< flags that snapshot contains double-precision instead of single precision */
    int      flag_ic_info;              /*< flag to inform whether IC files are generated with Zeldovich approximation or 2nd order lagrangian perturbation theory initial conditions. */\
    float    lpt_scalingfactor;         /*< scaling factor for 2lpt initial conditions */
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 4*12];  /* fills to 256 Bytes */
    
    Gadget_header();   // class constructor
    void print();      // print the header values to the standard output
    void swapBytes();  // swaps the endianess of the gadget header values
    
    int totalParticles();        // get total number of particles in a given header
    int totalTotalParticles();   // get total number of particles in the snapshot
    int gasDMParticles();        // return the total number of gas and DM particles
    
    bool differentMass();        // true if at least one species of particle has particles with different masses
    bool differentMass(int const n);        // true if particle species 'n' has particles with different masses
    int totalParticlesDifferentMass();      // get the number of particles with different masses
    void gasOnly();              // keep in the header only gas relevant information
    void darkMatterOnly();       // keep in the header only DM relevant information
    
    //internal functions
    void keepSpeciesOnly(int const n);   // keep in the header only species 'n' relevant information
};



// Function used to swap bytes between different endianness
inline void ByteSwap(unsigned char * b, int n)
{
    register int i = 0;
    register int j = n-1;
    while (i<j)
    {
        std::swap(b[i], b[j]);
        i++, j--;
    }
}
template <typename T>
void ByteSwapArray(T *x, size_t const elements)
{
    int size = sizeof(x[0]);
    for (size_t i=0; i<elements; ++i)
        ByteSwap( (unsigned char *) &(x[i]), size );
}

#define BYTESWAP(x) ByteSwap( (unsigned char *) &x, sizeof(x) )

#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3)           { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }





#include "gadgetParticles.h"


// class for the cosmological parameters
struct Cosmos
{
    double omegaBaryons;    //baryonic density
    double omegaMatter;     //matter density
    double omegaLambda;     // lambda density
    double h0;              // Hubble parameter h0 where the Hubble constat H0 = 100 km/s/Mpc * h0
    double R8;              // the radius at which sigma8 is computed (usually 8 h^-1 Mpc)
    double sigma8;          // the value of sigma8 at present
    double n;               // index of the primordial power fluctuations
    double amplitude;       // amplitude of the power spectrum
    
    Cosmos( void )
    {
        omegaBaryons = 0.; omegaMatter = 0.; omegaLambda = 0.;  h0 = 0.; R8 = 0.; sigma8 = 0.; n = 0.; amplitude = 0.;
    }
    Cosmos(Gadget_header const & gadgetHeader)
    {
        omegaBaryons = 0.;
        omegaMatter = gadgetHeader.Omega0;
        omegaLambda = gadgetHeader.OmegaLambda;
        h0 = gadgetHeader.HubbleParam;
        R8 = 0.;
        sigma8 = 0.;
        n = 0.;
        amplitude = 0.;
    }
    Cosmos(double const omegaB, Gadget_header const & gadgetHeader, double const R, double const sigma, double const nPrimordial, double const amplPower = 0. )
    {
        omegaBaryons = omegaB;
        omegaMatter = gadgetHeader.Omega0;
        omegaLambda = gadgetHeader.OmegaLambda;
        h0 = gadgetHeader.HubbleParam;
        R8 = R;
        sigma8 = sigma;
        n = nPrimordial;
        amplitude = amplPower;
    }
    Cosmos(double const omegaB, double const omegaM, double const omegaL, double const hubble0, double const R, double const sigma, double const nPrimordial, double const amplPower = 0. )
    {
        omegaBaryons = omegaB;
        omegaMatter = omegaM;
        omegaLambda = omegaL;
        h0 = hubble0;
        R8 = R;
        sigma8 = sigma;
        n = nPrimordial;
        amplitude = amplPower;
    }
};

#endif

