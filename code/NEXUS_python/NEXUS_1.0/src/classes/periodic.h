/*

This is a header file that implements several functions that are very usefull when dealing with data in a periodic box.

*/


#ifndef PERIODIC_HEADER
#define PERIODIC_HEADER

#include <iostream>
#include <new>
#include <miscellaneous.h>
#include <boost/assert.hpp>

#include <defines.h>


// returns the index of a periodic grid
template <typename T>
inline T periodicGrid(T index, T nGrid)
{
    if (index>=nGrid) return index-nGrid;
    else if (index<T(0)) return nGrid+index;
    else return index;
}


// function that returns the periodic offsets when translating from periodic to non-periodic boxes
template <typename T>
void periodicOffsets(T offset[][NO_DIM], T *length)
{
#if NO_DIM==2
    for (int i1=0; i1<3; ++i1)
        for (int i2=0; i2<3; ++i2)
        {
            offset[i1*3+i2][0] = T(i1-1) * length[0];
            offset[i1*3+i2][1] = T(i2-1) * length[1];
        }
#elif NO_DIM==3
    for (int i1=0; i1<3; ++i1)
        for (int i2=0; i2<3; ++i2)
            for (int i3=0; i3<3; ++i3)
            {
                offset[i1*9+i2*3+i3][0] = T(i1-1) * length[0];
                offset[i1*9+i2*3+i3][1] = T(i2-1) * length[1];
                offset[i1*9+i2*3+i3][2] = T(i3-1) * length[2];
            }
#endif
}


/* This function copies additional particle to padd the boundaries of a periodic box up ot the 'boundary' distance from the periodic box boundaries. */
template <typename T>
Int periodicBoundary(Array<T,2> &inputData,
                     T boundary,
                     Box<T,NO_DIM> &box,
                     std::vector<T> &outputData,
                     std::vector<Int> *paddedIndices = NULL,
                     bool const VERBOSE=true)
{
    if (VERBOSE)
        std::cout << "Extending the periodic box positions to a non-periodic box with boundaries of " << boundary << " Mpc/h ..." << std::flush;
    std::vector<Int> tempIntVector;
    if ( paddedIndices==NULL ) paddedIndices = &tempIntVector;
    
    T length[NO_DIM];    //stores the box length
    box.length( length );
    int const noPoints = inputData.axisSize(0);
    int reserveSize =  NO_DIM * noPoints * int( std::pow(1.+3*boundary/(length[0]+2*boundary),3) );    //allocate some memory for the particles in the boundary
    outputData.reserve( reserveSize );
    for (int i=0; i<noPoints; ++i)
        for (int j=0; j<NO_DIM; ++j)
            outputData.push_back( inputData(i,j) );
    
    
    // variables to keep track of the offsets when applying the periodic boundaries
    int const noOffsets = (NO_DIM==2 ? 9:27);
    T offset[noOffsets][NO_DIM];
    periodicOffsets<T>( offset, length );

    
    // to speed the search first compute which particles are inside the inner box => these particles do not have any boundary contributions
    Box<T,NO_DIM> innerBox( box );  // particles inside this box do not have any boundary contributions
    Box<T,NO_DIM> outerBox( box );  // the coordinates of the larger box that includes the boundaries too
    innerBox.addPadding( -boundary );
    outerBox.addPadding( boundary );
    
    Int noFinalPoints = noPoints;     //the final number of points
    T pos[NO_DIM], tempPos[NO_DIM];
    for (int i=0; i<noPoints; ++i)
    {
        for (int j=0; j<NO_DIM; ++j)
            pos[j] = inputData(i,j);
        if ( innerBox.isPointInBox(pos) )   //this point will not have boundary contributions
            continue;
        for (int i1=0; i1<noOffsets; ++i1)
        {
            for (int j=0; j<NO_DIM; ++j)
                tempPos[j] = pos[j] + offset[i1][j];
            if ( not outerBox.isPointInBox(tempPos) )
                continue;
            ++noFinalPoints;
            for (int j=0; j<NO_DIM; ++j)
                outputData.push_back( tempPos[j] );
            paddedIndices->push_back( i );
        }
    }
    
    if (VERBOSE)
        std::cout << " Done.\n" << "\tAdded " << noFinalPoints-noPoints << " additional points (" << std::setprecision(4) << double( 100.*(noFinalPoints-noPoints)/noPoints ) << "\%) due to the periodic boundaries.\n";
    return noPoints;
}



/* This function checks that the points are inside the periodic box, otherwise is offsets them by the box length to fold them back inside the box. */
template <typename T>
void periodicFolding(Array<T,2> &inputData,
                     Box<T,NO_DIM> &box,
                     bool const VERBOSE=true)
{
    if (VERBOSE)
        std::cout << "Folding back the points into the periodic box  ..." << std::flush;
    
    T length[NO_DIM];    //stores the box length
    box.length( length );
    size_t const noPoints = inputData.axisSize(0);
    for (size_t i=0; i<noPoints; ++i)
    {
        if ( box.isPointInBox( &(inputData(i,0)) ) )
            continue;
        for (int j=0; j<NO_DIM; ++j)
        {
            if ( inputData(i,j)<box[2*j] )
                inputData(i,j) += length[j];
            else if ( inputData(i,j)>box[2*j+1] )
                inputData(i,j) -= length[j];
        }
    }
    
    if (VERBOSE)
        std::cout << " Done.\n";
}





#endif
