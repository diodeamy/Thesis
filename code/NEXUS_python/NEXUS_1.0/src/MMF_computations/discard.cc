#include <list>

#include <defines.h>
#include <gadget_header.h>
#include <miscellaneous.h>
#include <gadget_header.h>
#include <array.h>
//#include <hopGroup.h>
using namespace std;


/* Discards the particles inside a voxel mask which covers all the simulation box. 'gridSize[3]' gives the dimmensions of the mask along the 3 directions.
*/
template <typename T>
void discardParticles(Array<T,1> &p,
                      Gadget_header *gadgetHeader,
                      Array<shortInt,NO_DIM> &mask,
                      shortInt threshold,
                      Real *box)
{
    if ( NO_DIM!=3 ) throwError( "Function 'discardParticles' gives correct results only for 3 dimensions." );
    cout << "Discarding the particles inside the MMF response mask ... " << flush;
    
    // grid spacing along the 3 directions
    Int const noParticles = p.totalSize();
    Array<Int,1> gridSize = mask.getSize<Int>();
    for (int i=0; i<6; ++i)
        box[i] *= 1000.;
    Real const dx[3] = { (box[1]-box[0]) / Real(gridSize[0]), (box[3]-box[2]) / Real(gridSize[1]), (box[5]-box[4]) / Real(gridSize[2]) };
    // some additional constants
    Int const nZ = gridSize[2];
    Int const nYnZ = gridSize[1] * nZ;
    
    // define array which keeps track if a particle should be discarded or not
    Array<bool,1> discard( noParticles );
    
    // find which particle should be discarded and which not
    for (Int i=0; i<noParticles; ++i)
    {
        bool validCell = true;
        int temp[3];
        for (int j=0; j<3; ++j)
        {
            temp[j] = int( floor( (p[i].pos[j]-box[2*j]) / dx[j] ) );
            if ( temp[j]<0 or temp[j]>=gridSize[j] )
                validCell = false;
        }
        if ( not validCell ) continue;
        
         
        int gridPos = temp[0] * nYnZ + temp[1] * nZ + temp[2];  //particle position in the mask's grid
        
        if ( mask[gridPos]>=threshold )
            discard[i] = true;
        else
            discard[i] = false;
    }
    
    // rearrange the particles in the particle array and discard the rest
    Int writePos = 0, temp = 0;
    Int newParticleNo[6];
    for (int j=0; j<6; ++j)     //now loop over each particle species
    {
        for (Int i=temp; i<temp+gadgetHeader->npartTotal[j]; ++i)
            if ( not discard[i] )
            {
                p[writePos] = p[i];
                ++writePos;
            }
        newParticleNo[j] = writePos;
        temp += gadgetHeader->npartTotal[j];
    }
    
    // update the gadget header to contain the new particle information
    gadgetHeader->npart[0] = newParticleNo[0];
    for (int i=1; i<6; ++i)
        gadgetHeader->npart[i] = newParticleNo[i]-newParticleNo[i-1];
    for (int i=0; i<6; ++i)
        gadgetHeader->npartTotal[i] = gadgetHeader->npart[i];
    
    cout << "Done.\n";
    Int tempNo = noParticles - (newParticleNo[0]+newParticleNo[1]);
    cout << "Discarded " << tempNo << " which represent " << tempNo/float(noParticles)*100. << "\% of the total number of particles.\n";
}




/*
// Returns the haloes inside a given mask
void hopHaloesInMask(Array<HopGroup,1> h,
                     Real const boxSize,
                     Array<bool,NO_DIM> mask,
                     list<HopGroup> *result)
{
    if ( NO_DIM!=3 ) throwError( "Function 'discardParticles' gives correct results only for 3 dimensions." );
    cout << "Searching for the haloes inside the MMF response mask ... " << flush;
    
    // grid spacing along the 3 directions
    Int const noHaloes = h.totalSize();
    Array<Int,1> gridSize = mask.getSize<Int>();
    Real const dx[3] = { boxSize / Real(gridSize[0]), boxSize / Real(gridSize[1]), boxSize / Real(gridSize[2]) };
    // some additional constants
    Int const nZ = gridSize[2];
    Int const nYnZ = gridSize[1] * nZ;
    
    for (Int i=0; i<noHaloes; ++i)
    {
        Int temp[NO_DIM];
        for (int j=0; j<NO_DIM; ++j)
            temp[j] = Int( h[i].pos[j] / dx[j] );
        
        Int gridPos = temp[0] * nYnZ + temp[1] * nZ + temp[2];  //particle position in the mask's grid
        
        if ( mask[gridPos] )
            result->push_back( h.at(i) );
    }
    cout << "Done.\n";
    cout << "Found " << result->size() << " haloes inside the mask. This represents " << float(result->size())/noHaloes*100. << "\% of the total number of " << noHaloes << " haloes.\n";
}*/




