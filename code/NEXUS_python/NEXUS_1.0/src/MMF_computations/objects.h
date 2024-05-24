
#ifndef OBJECTS_HEADER
#define OBJECTS_HEADER


#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <limits>

//#include "MMF_computations.h"
#include "../defines/defines.h"
#include "../classes/array.h"
#include "../miscellaneous/misc.h"




//! Functions in 'objects.cc'
// identify the distinct objects
void identifyDistinctObjects(Array<int,3> *mask,
                             int const neighborFindingMethod,
                             std::vector<Int> *objectSize);
void identifyDistinctObjects(Array<int,3> *mask,
                             int const noOldObjects,
                             vector<Int> &newCells,
                             int const neighborFindingMethod,
                             std::vector<Int> *objectSize);

/* compute the distinct objects starting from:
        1 = a response and threshold
        2 = from a short int mask with 0=empty and 1=taken
        3 = from an int mask with -1=empty and 0=taken - returns the results in the same mask
NOTE: it does not order the object labels according to their size!
*/
Int compactObjects(Array<Real,3> &response,
                   Real const threshold,
                   int const neighborFindingMethod,
                   Array<int,3> *mask,
                   std::vector<Int> *objectSize);
Int compactObjects(Array<shortInt,3> &mask,
                   int const neighborFindingMethod,
                   Array<int,3> *outMask,
                   std::vector<Int> *objectSize);
Int compactObjects(Array<int,3> *mask,
                   int const neighborFindingMethod,
                   std::vector<Int> *objectSize);

/* Computes compact objects as above, but on top of that:
        - orders object labels in decreasing order of their size (large obejcts have small label number and small objects have large labels)
        - discards the objects below a certain size
        - the size is given by the object volume (in voxel counts)
*/
Int significantObjects(Array<Real,3> &response,
                       Real const threshold,
                       int const neighborFindingMethod,
                       Array<int,3> *mask,
                       std::vector<Int> *objectSize,
                       Int const minObjectSize);
Int significantObjects(Array<shortInt,3> &mask,
                       int const neighborFindingMethod,
                       Array<int,3> *outMask,
                       std::vector<Int> *objectSize,
                       Int const minObjectSize);
Int significantObjects(Array<int,3> *mask,
                       int const neighborFindingMethod,
                       std::vector<Int> *objectSize,
                       Int const minObjectSize);

/* The following need a mass array (array giving the mass in each cell), so their size will be ordered according to their mass and not volume. */
Int significantObjects(Array<Real,3> &response,
                       Real const threshold,
                       int const neighborFindingMethod,
                       Array<int,3> *mask,
                       Array<Real,3> &mass,
                       std::vector<double> *objectSize,
                       double const minObjectSize);
Int significantObjects(Array<shortInt,3> &mask,
                       int const neighborFindingMethod,
                       Array<int,3> *outMask,
                       Array<Real,3> &mass,
                       std::vector<double> *objectSize,
                       double const minObjectSize);
Int significantObjects(Array<int,3> *mask,
                       int const neighborFindingMethod,
                       Array<Real,3> &mass,
                       std::vector<double> *objectSize,
                       double const minObjectSize);


/* Returns the clean response map for the given input parameters: for the given maximum response map with values >= 'threshold' and keeping only the objects with sizes >= 'minSize'. */
void cleanResponse(Array<Real,3> &response,
                   Real const threshold,
                   int const neighborFindingMethod,
                   Int const minObjectSize,
                   Array<shortInt,3> *mask,
                   bool const VERBOSE = true);
void cleanResponse(Array<Real,3> &response,
                   Real const threshold,
                   int const neighborFindingMethod,
                   Array<Real,3> &mass,
                   double const minObjectSize,
                   Array<shortInt,3> *mask,
                   bool const VERBOSE = true);





/* Relabels the Cosmic Web objects according to their size (mass or volume): larger objects first, smaller latter. */
template <typename T>
void relabelObjects(Array<int,NO_DIM> *mask,
                    std::vector<T> *objectSize,
                    bool const VERBOSE = true)
{
    if (VERBOSE) std::cout << "Relabeling the Cosmic Web objects according to their size (mass or volume) ... " << std::flush;
    
    std::vector< std::pair<T,int> > orderObjects;
    orderObjects.reserve( objectSize->size() );
    for (size_t i=0; i<objectSize->size(); ++i)
        orderObjects.push_back( std::make_pair((*objectSize)[i],i) );
    std::sort( orderObjects.begin(), orderObjects.end() );   //sorts the objects according to their volume
    std::reverse( orderObjects.begin(), orderObjects.end() );
    
    // find the new labels of each MMF group and update their new size in 'objectSize'
    Array<int,1> newLabel( objectSize->size() );
    std::vector<T> oldSize = *objectSize;
    for (size_t i=0; i<orderObjects.size(); ++i)
    {
        newLabel[ orderObjects[i].second ] = i;
        (*objectSize)[i] = oldSize[ orderObjects[i].second ];
    }
    
    // relabel the valid lattice points in 'mask' according to the new labels
    for (size_t i=0; i<mask->size(); ++i)
        if ( (*mask)[i]>=0 )
            (*mask)[i] = newLabel[ (*mask)[i] ];
    
    if (VERBOSE) std::cout << "Done.\n";
}

/* Relabels the objects according to their size, like the function above. It discards any objects smaller than the size threshold. */
template <typename T>
void relabelSignificantObjects(Array<int,NO_DIM> *mask,
                               std::vector<T> *objectSize,
                               T const minSize,
                               bool const VERBOSE = true)
{
    if (VERBOSE) std::cout << "Relabeling the Cosmic Web objects according to their size (mass or volume). Keeping only the significant objects, with sizes larger than " << minSize << " ... " << std::flush;
    
    std::vector< std::pair<T,int> > orderObjects;
    orderObjects.reserve( objectSize->size() );
    for (size_t i=0; i<objectSize->size(); ++i)
        orderObjects.push_back( std::make_pair((*objectSize)[i],i) );
    std::sort( orderObjects.begin(), orderObjects.end() );   //sorts the objects according to their volume
    std::reverse( orderObjects.begin(), orderObjects.end() );
    
    // find the new labels of each MMF group and update their new size in 'objectSize'
    Array<int,1> newLabel( objectSize->size() );
    std::vector<T> oldSize = *objectSize;
    int noValidObjects = 0;
    for (size_t i=0; i<orderObjects.size(); ++i)
    {
        (*objectSize)[i] = oldSize[ orderObjects[i].second ];
        if ( (*objectSize)[i]>=minSize )        //if object significant, keep it
        {
            newLabel[ orderObjects[i].second ] = i;
            ++noValidObjects;
        }
        else                                    //if object insignificant, discard it
            newLabel[ orderObjects[i].second ] = -1;
    }
    objectSize->resize( noValidObjects );
    
    // relabel the valid grid points in 'mask' according to the new labels
    for (size_t i=0; i<mask->size(); ++i)
        if ( (*mask)[i]>=0 )
            (*mask)[i] = newLabel[ (*mask)[i] ];
    
    if (VERBOSE) std::cout << "Done.\n";
}



/* Returns the selection with all elements '>=' than 'minThreshold' and '<' than 'maxThreshold' as valid (value given by 'valid') and everything else as invalid. */
template <typename T1, typename T2>
void selection(Array<T1,NO_DIM> &inputData,
               T1 const minThreshold,
               T1 const maxThreshold,
               Array<T2,NO_DIM> *select,
               T2 const valid = T2(1),
               T2 const invalid = T2(0))
{
    if ( select->size()!=inputData.size() )
        throwError( "You are trying to do operations with two arrays of different sizes. Error in function 'maskData'." );
    if ( valid==invalid )
        throwError( "In function 'maskData'. Both the 'valid' and 'invalid' parameters have the same value. They should have different values to distinguish between valid and invalid grid cells." );
    
    Int const gridSize = inputData.size();
    for (Int i=0; i<gridSize; ++i)
        if ( inputData[i]>=minThreshold and inputData[i]<maxThreshold )
            (*select)[i] = valid;
        else
            (*select)[i] = invalid;
}
/* Returns the grid indices of all elements with 'minThreshold' <= data < 'maxThreshold'. */
template <typename T>
void selection(Array<T,NO_DIM> &inputData,
               T const minThreshold,
               T const maxThreshold,
               std::vector<Int> *select)
{
    Int const gridSize = inputData.size();
    for (Int i=0; i<gridSize; ++i)
        if ( inputData[i]>=minThreshold and inputData[i]<maxThreshold )
            select->push_back( i );
}



/* Returns the mass of all the objects. */
template <typename T, typename T2>
void objectsMass(Array<int,NO_DIM> &mask,
                 Array<T,3> &mass,
                 std::vector<T2> *objectSize,
                 int const noObjects)
{
    if ( mask.size()!=mass.size() )
        throwError( "You are trying to do operations with two arrays of different sizes. Error in function 'objectsMass'." );
    
    objectSize->assign( noObjects, T2(0.) );
    for (size_t i=0; i<mask.size(); ++i)
        if ( mask[i]>=0 )
            (*objectSize)[ mask[i] ] += T2( mass[i] );
}




#endif
