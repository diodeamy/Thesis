#include "gadget_header.h"
#include <miscellaneous.h>
using namespace std;



/* Class constructor */
Gadget_header::Gadget_header()
{
    for (int i=0; i<6; ++i)
    {
        npart[i] = 0;
        mass[i] = 0.;
        npartTotal[i] = 0;
        num_total_particles_hw[i] = 0;
    }
    time = -1.; redshift = -1.;
    flag_sfr = 0; flag_feedback = 0; flag_cooling = 0;
    num_files = 1;
    BoxSize = -1; Omega0 = -1.; OmegaLambda = -1.; HubbleParam = -1.;
    flag_stellarage = 0; flag_metals = 0; flag_entropy_instead_u = 0; flag_doubleprecision = 0; flag_ic_info = 0;
    lpt_scalingfactor = -1.;
    
    initializeCharArray( fill, sizeof(fill) );
}


/* Function that prints the Gadget header.
*/
void Gadget_header::print()
{
    cout << "The header of the Gadget file contains the following info:\n" 
            << "npart[6]     =  " << npart[0] << "  " << npart[1] << "  " << npart[2] << "  " << npart[3] << "  " <<  npart[4] << "  " <<  npart[5] << "\n"
            << "mass[6]      =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time         =  " << time << "\n"
            << "redshift     =  " << redshift << "\n"
            << "flag_sfr     =  " << flag_sfr << "\n"
            << "flag_feedback=  " << flag_feedback << "\n"
            << "npartTotal[6]=  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "  " << "\n"
            << "flag_cooling =  " << flag_cooling << "\n"
            << "num_files    =  " << num_files << "\n"
            << "BoxSize      =  " << BoxSize << "\n"
            << "Omega0       =  " << Omega0 << "\n"
            << "OmegaLambda  =  " << OmegaLambda << "\n"
            << "h            =  " << HubbleParam << "\n\n";
}


/* Change the endianness. */
void Gadget_header::swapBytes()
{
    ByteSwapArray( npart, 6 );
    ByteSwapArray( mass, 6 );
    BYTESWAP( time );
    BYTESWAP( redshift );
    BYTESWAP( flag_sfr );
    BYTESWAP( flag_feedback );
    ByteSwapArray( npartTotal, 6 );
    BYTESWAP( flag_cooling );
    BYTESWAP( num_files );
    BYTESWAP( BoxSize );
    BYTESWAP( Omega0 );
    BYTESWAP( OmegaLambda );
    BYTESWAP( HubbleParam );
    BYTESWAP( flag_stellarage );
    BYTESWAP( flag_metals );
    ByteSwapArray( num_total_particles_hw, 6 );
    BYTESWAP( flag_entropy_instead_u );
    BYTESWAP( flag_doubleprecision );
    BYTESWAP( flag_ic_info );
    BYTESWAP( lpt_scalingfactor );
}



// Compute total number of particles in a given header.
int Gadget_header::totalParticles()
{
    int temp = 0;
    for (int i=0; i<6; ++i)
        temp += npart[i];
    return temp;
}


// Compute total number of particles in the snapshot.
int Gadget_header::totalTotalParticles()
{
    int temp = 0;
    for (int i=0; i<6; ++i)
        temp += npartTotal[i];
    return temp;
}


// Compute total number of particles in the snapshot.
int Gadget_header::gasDMParticles()
{
    int temp = 0;
    for (int i=0; i<2; ++i)
        temp += npartTotal[i];
    return temp;
}



// Check if there are particles of the same species with different masses. Return true if that is the case.
bool Gadget_header::differentMass()
{
    for (int i=0; i<6; ++i)
        if (mass[i]==0. and npart[i]!=0)
            return true;
    return false;
}


// Check if there are particles of species 'n' with different masses. Return true if that is the case.
bool Gadget_header::differentMass(int const n)
{
    intervalCheck( n, 0, 5, "'n' in function 'Gadget_header::differentMass'" ); //check the "n" value
    if (mass[n]==0. and npart[n]!=0)
        return true;
    return false;
}


// Gives the total number of particles in a header that have different masses.
int Gadget_header::totalParticlesDifferentMass()
{
    int temp = 0;
    for (int i=0; i<6; ++i)
        if (mass[i]==0. and npart[i]!=0)
            temp += npart[i];
    return temp;
}




/* Keep only information relevant for dark matter (DM) particles only. Suitable to use when outputing in a file only the DM particles.
*/
void Gadget_header::darkMatterOnly()
{
    this->keepSpeciesOnly( 1 );
}


/* Keep only information relevant for gas particles only. Suitable to use when outputing in a file only the gas particles.
*/
void Gadget_header::gasOnly()
{
    this->keepSpeciesOnly( 0 );
}


/* Keep only the relevant header information for species "n".
*/
void Gadget_header::keepSpeciesOnly(int const n)
{
    intervalCheck( n, 0, 5, "'n' in function 'Gadget_header::speciesOnly'" ); //check the "n" value
    for (int i=0; i<n; ++i)
    {
        npart[i] = 0;
        mass[i] = 0.;
        npartTotal[i] = 0;
    }
    for (int i=n+1; i<6; ++i)
    {
        npart[i] = 0;
        mass[i] = 0.;
        npartTotal[i] = 0;
    }
}


