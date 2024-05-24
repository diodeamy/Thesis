
#include "objects.h"

using namespace std;




// Implements a periodic index on the grid
inline int pI(int index, int period)
{
    if (index<0)
        return period+index;
    else if (index>=period)
        return index-period;
    return index;
}



/* Returns the 2D and 3D neighbors of a grid cell. The parameters are:
        nx, ny, nz = the neighbor indices along the 3 dimensions
        neighborFindingMethod = what is a neighbor (1=the closest 6 neighbors, 2=the closest 26 neighbors)
        fullNeighbors = if true it returns all the neighbors, otherwise only half of them
*/
void neighborIndices(int nx[], int ny[], int nz[],
                     int const neighborFindingMethod = 1,
                     bool const fullNeighbors = true)
{
    int indicesFull_1[][6]  = { {-1,1,0,0,0,0}, {0,0,-1,1,0,0}, {0,0,0,0,-1,1} };
    int indicesFull_2[][26] = { {-1,-1,-1,  0, 0, 0,  1, 1, 1,   -1,-1,-1,  0, 0,  1, 1, 1,   -1,-1,-1,  0, 0, 0,  1, 1, 1 },
                                {-1, 0, 1, -1, 0, 1, -1, 0, 1,   -1, 0, 1, -1, 1, -1, 0, 1,   -1, 0, 1, -1, 0, 1, -1, 0, 1 },
                                {-1,-1,-1, -1,-1,-1, -1,-1,-1,    0, 0, 0,  0, 0,  0, 0, 0,    1, 1, 1,  1, 1, 1,  1, 1, 1 } };
    
    int indicesHalf_1[][3]  = { {-1,0,0}, {0,-1,0}, {0,0,-1} };
    int indicesHalf_2[][13] = { {-1,-1,-1,  0, 0, 0,  1, 1, 1,  -1,-1,-1,  0 },
                                {-1, 0, 1, -1, 0, 1, -1, 0, 1,  -1, 0, 1, -1 },
                                {-1,-1,-1, -1,-1,-1, -1,-1,-1,   0, 0, 0,  0 } };
    
    if ( fullNeighbors )
    {
        int const noElements = neighborFindingMethod==1 ? 6 : 26;
        if (neighborFindingMethod==1)
            for (int i1=0; i1<noElements; ++i1)
            {
                nx[i1] = indicesFull_1[0][i1];
                ny[i1] = indicesFull_1[1][i1];
                nz[i1] = indicesFull_1[2][i1];
            }
        else if (neighborFindingMethod==2)
            for (int i1=0; i1<noElements; ++i1)
            {
                nx[i1] = indicesFull_2[0][i1];
                ny[i1] = indicesFull_2[1][i1];
                nz[i1] = indicesFull_2[2][i1];
            }
    }
    else
    {
        int const noElements = neighborFindingMethod==1 ? 3 : 13;
        if (neighborFindingMethod==1)
            for (int i1=0; i1<noElements; ++i1)
            {
                nx[i1] = indicesHalf_1[0][i1];
                ny[i1] = indicesHalf_1[1][i1];
                nz[i1] = indicesHalf_1[2][i1];
            }
        else if (neighborFindingMethod==2)
            for (int i1=0; i1<noElements; ++i1)
            {
                nx[i1] = indicesHalf_2[0][i1];
                ny[i1] = indicesHalf_2[1][i1];
                nz[i1] = indicesHalf_2[2][i1];
            }
    }
}





/* This function takes an input mask with values: -1=empty voxels and 0=occupied voxels and finds all the distinct objects. Two voxels are part of the same object if they are neighbors. The argument 'neighborFindingMethod'=1/2 specifies that each voxel has 6/26 neighbors - the closest ones in the grid.
*/
void identifyDistinctObjects(Array<int,3> *mask,
                             int const neighborFindingMethod,
                             vector<Int> *objectSize)
{
    // get the neighbor array
    intervalCheck<int>( neighborFindingMethod, 0, 3, "'neighborFindingMethod' argument in function 'identifyDistinctObjects'" );
    int const neighbors = neighborFindingMethod==1 ? 3 : 13;  //use only half of the neighbors - the ones to the right
    int nx[neighbors], ny[neighbors], nz[neighbors];
    neighborIndices( nx, ny, nz, neighborFindingMethod, false );    //get the upper half neighbors 
    
    // assign to each non-empty cell an 'object'
    int count = 1;         //starting counting objects from 1
    int const NX = mask->getSize(0),
              NY = mask->getSize(1),
              NZ = mask->getSize(2);
    Int const NYZ = Int(NY*NZ);
    for (int i1=0; i1<NX; ++i1)
        for (int i2=0; i2<NY; ++i2)
            for (int i3=0; i3<NZ; ++i3)
            {
                if ( mask->at(i1,i2,i3)<0 )      // skip if cell is 'empty'
                    continue;
                Int const mainIndex = i1*NYZ + Int(i2*NZ + i3);
                
                for (int j=0; j<neighbors; ++j)      // loop over all neighbors
                {
                    Int const index = pI(i1+nx[j],NX)*NYZ + Int(pI(i2+ny[j],NY)*NZ + pI(i3+nz[j],NZ));
                    if ( mask->at(index)>0 )// if neighbor='object', compute further
                    {
                        (*mask)[mainIndex] = (*mask)[index];
                        break;
                    }
                }
                if ( mask->at(mainIndex)==0 )      // if current cell has no 'object' neighbors, assign a new 'object' to it
                    (*mask)(mainIndex) = count++;
            }
    
    // loop again over all cells and merge neighboring 'objects'
    list< pair<int,int> > mergeObjects;// variable to keep track what 'objects' to merge - always merge object with higher id to the one with lower id
    for (int i1=0; i1<NX; ++i1)
        for (int i2=0; i2<NY; ++i2)
            for (int i3=0; i3<NZ; ++i3)
            {
                if ( mask->at(i1,i2,i3)<0 )      // skip if cell is 'empty'
                    continue;
                Int const mainIndex = i1*NYZ + Int(i2*NZ + i3);
                
                for (int j=0; j<neighbors; ++j)      // loop over all neighbors
                {
                    Int const index = pI(i1+nx[j],NX)*NYZ + Int(pI(i2+ny[j],NY)*NZ + pI(i3+nz[j],NZ));
                    if ( mask->at(index)>0  and mask->at(index)!=mask->at(mainIndex) ) // if two different 'objects' are neighbor
                    {
                        int obj1 = Min<int>( mask->at(mainIndex), mask->at(index) );
                        int obj2 = Max<int>( mask->at(mainIndex), mask->at(index) );
                        if ( mergeObjects.empty() ) mergeObjects.push_back( make_pair(obj1,obj2) );
                        if ( mergeObjects.back()!=make_pair(obj1,obj2) )
                            mergeObjects.push_back( make_pair(obj1,obj2) );
                    }
                }
            }
    int noObjects = count-1;        //total number of 'objects; found - before merging adjacent objects
    
    // sort the pair array and keep only unique entries
    mergeObjects.sort();
    mergeObjects.unique();
    
    
    // decide what neighboring objects must be merged
    Array<int,1> finalObjectId( noObjects+1 ); // a non-zero value of it tells that the object 'i' must be merged to the value stored at entry 'i' 
    finalObjectId.assign( 0 );
    while ( mergeObjects.size()!=0 )        //loop while there are still objects to merge
    {
//        cout << iter << "   " << mergeObjects.size() << "\n";
        list< pair<int,int> > temp;
        for (list< pair<int,int> >::iterator it=mergeObjects.begin(); it!=mergeObjects.end() ; ++it)
        {
            if ( finalObjectId[ it->second ]==0 and finalObjectId[ it->first ]==0 )  // if both first and second object isn't merged to something else
                finalObjectId[ it->second ] = it->first;
            else if ( finalObjectId[ it->second ]==0 )          // if first object is already merged, but not second
                finalObjectId[ it->second ] = finalObjectId[ it->first ];
            else if ( finalObjectId[ it->second ]!=it->first )  // if also second object is merged, must keep track of this
            {
                if ( finalObjectId[ it->second ]<it->first )
                    temp.push_back( make_pair( finalObjectId[ it->second ], it->first ) );
                else temp.push_back( make_pair( it->first, finalObjectId[ it->second ] ) );
            }
        }
        temp.sort();
        temp.unique();
        
        mergeObjects = temp;    // the new list of objects still to be merged
    }
    
    
    // now entry 'i' of 'finalObjectId', if !=0, tells that the group 'i' must be merged with the group 'finalObjectId[i]'<'i'
    // relabel the remaining unique groups from 0 to 'number unique groups-1'
    count = 0;
    for (int i=1; i<noObjects+1; ++i)
    {
        if ( finalObjectId[i]==0 )   //if this group will not be merged
            finalObjectId[i] = count++;
        else
            finalObjectId[i] = finalObjectId[ finalObjectId[i] ];
    }
    
    
    // now relabel the mask 'non-empty' values with the updated object ID and count how many cells are in each object
    noObjects = count;                 // number of unique objects in the mask
    objectSize->assign( noObjects, Int(0) );// counts how many cells are in each unique object
    for (Int i=0; i<mask->totalSize(); ++i)
    {
        if ( (*mask)[i]<0 )
            continue;
        (*mask)[i] = finalObjectId( (*mask)[i] );
        ++(objectSize->at( (*mask)[i] ));  // increase cell count for the given object
    }
}

/* Does similar computations as the above function, but it takes an already existing mask and adds additional valid voxels. It than recomputes the distinct objects in the input mask. */
void identifyDistinctObjects(Array<int,3> *mask,
                             int const noOldObjects,
                             vector<Int> &newCells,
                             int const neighborFindingMethod,
                             vector<Int> *objectSize)
{
    // get the neighbor array
    intervalCheck<int>( neighborFindingMethod, 0, 3, "'neighborFindingMethod' argument in function 'identifyDistinctObjects'" );
    int const neighbors = neighborFindingMethod==1 ? 6 : 26;  //use all the neighbors
    int nx[neighbors], ny[neighbors], nz[neighbors];
    neighborIndices( nx, ny, nz, neighborFindingMethod, true );
    
    
    // loop over the new cells to be add and assign a valid object label to each of them
    int count = noOldObjects;         // start counting objects from previous labels
    int const NX = mask->getSize(0),
              NY = mask->getSize(1),
              NZ = mask->getSize(2);
    Int const NYZ = NY*NZ;
    for (Int i=0; i<newCells.size(); ++i)
    {
        Int mainIndex = newCells[i];
        int i1 = int( mainIndex / NYZ );
        int i2 = int( (mainIndex-i1*NYZ) / NZ );
        int i3 = int( mainIndex - i1*NYZ - i2*NZ );
        
        for (int j=0; j<neighbors; ++j)      // loop over all neighbors
        {
            Int const index = pI(i1+nx[j],NX)*NYZ + Int(pI(i2+ny[j],NY)*NZ + pI(i3+nz[j],NZ));
            if ( mask->at(index)>=0 )// if neighbor='object', compute further
            {
                (*mask)(mainIndex) = (*mask)(index);
                break;
            }
        }
        if ( mask->at(mainIndex)<0 )      // if current cell has no 'object' neighbors, assign a new 'object' to it
            (*mask)(mainIndex) = count++;
    }
//    cout << "\n\t inside : oldObjects " << noOldObjects << "    new objects " << count << "\n" <<flush;
    
    // loop over the new cells and see if they merged new objects
    list< pair<int,int> > mergeObjects;// variable to keep track what 'objects' to merge - always merge object with higher id to the one with lower id
    for (Int i=0; i<newCells.size(); ++i)
    {
        Int mainIndex = newCells[i];
        int i1 = int( mainIndex / NYZ );
        int i2 = int( (mainIndex-i1*NYZ) / NZ );
        int i3 = int( mainIndex - i1*NYZ - i2*NZ );
        
        for (int j=0; j<neighbors; ++j)      // loop over all neighbors
        {
            Int const index = pI(i1+nx[j],NX)*NYZ + Int(pI(i2+ny[j],NY)*NZ + pI(i3+nz[j],NZ));
            if ( mask->at(index)>=0  and mask->at(index)!=mask->at(mainIndex) ) // if two different 'objects' are neighbors
            {
                int obj1 = Min<int>( mask->at(mainIndex), mask->at(index) );
                int obj2 = Max<int>( mask->at(mainIndex), mask->at(index) );
                if ( mergeObjects.empty() ) mergeObjects.push_back( make_pair(obj1,obj2) );
                if ( mergeObjects.back()!=make_pair(obj1,obj2) )
                    mergeObjects.push_back( make_pair(obj1,obj2) );
            }
        }
    }
    int noObjects = count;        //total number of 'objects; found - before merging adjacent objects
    
    // sort the pair array and keep only unique entries
    mergeObjects.sort();
    mergeObjects.unique();
    
    
    // decide what neighboring objects must be merged
    Array<int,1> finalObjectId( noObjects ); // a non-negative value of it tells that the object 'i' must be merged to the value stored at entry 'i' 
    finalObjectId.assign( -1 );
    
    while ( mergeObjects.size()!=0 )        //loop while there are still objects to merge
    {
        list< pair<int,int> > temp;
        for (list< pair<int,int> >::iterator it=mergeObjects.begin(); it!=mergeObjects.end() ; ++it)
        {
            if ( finalObjectId[ it->second ]==-1 and finalObjectId[ it->first ]==-1 )  // if both first and second object isn't merged to something else
                finalObjectId[ it->second ] = it->first;
            else if ( finalObjectId[ it->second ]==-1 )          // if first object is already merged, but not second
                finalObjectId[ it->second ] = finalObjectId[ it->first ];
            else if ( finalObjectId[ it->second ]!=it->first )  // if also second object is merged, must keep track of this
            {
                if ( finalObjectId[ it->second ]<it->first )
                    temp.push_back( make_pair( finalObjectId[ it->second ], it->first ) );
                else temp.push_back( make_pair( it->first, finalObjectId[ it->second ] ) );
            }
        }
        temp.sort();
        temp.unique();
        
        mergeObjects = temp;    // the new list of objects still to be merged
    }
    
    
    // now entry 'i' of 'finalObjectId', if !=0, tells that the group 'i' must be merged with the group 'finalObjectId[i]'<'i'
    // relabel the remaining unique groups from 0 to 'number unique groups-1'
    count = 0;
    for (int i=0; i<noObjects; ++i)
    {
        if ( finalObjectId[i]==-1 )   //if this group will not be merged
            finalObjectId[i] = count++;
        else
            finalObjectId[i] = finalObjectId[ finalObjectId[i] ];
    }
    
    
    // now relabel the mask 'non-empty' values with the updated object ID and count how many cells are in each object
    noObjects = count;                 // number of unique objects in the mask
    objectSize->assign( noObjects, Int(0) );// counts how many cells are in each unique object
    for (Int i=0; i<mask->totalSize(); ++i)
    {
        if ( (*mask)[i]<0 )
            continue;
        (*mask)[i] = finalObjectId( (*mask)[i] );
        ++(objectSize->at( (*mask)[i] ));  // increase cell count for the given object
    }
}




/* This function computes how many continuos groups are in an 'int' grid with values '-1'=vacuum and '0'=stuff.
'method' tells the function what means a continuos object:
    1 = check only left and right hand only neighbors along each axis
    2 = check all 26 neighbors (in 3D) which share a side or edge with the given cell
*/
Int computeCompactObjects(Array<Real,3> &response,
                          Real const threshold,
                          int const neighborFindingMethod,
                          Array<int,3> *mask,
                          vector<Int> *objectSize)
{
    // compute selection of cells which have response value larger than threshold
    selection( response, threshold, numeric_limits<Real>::max(), mask, 0, -1 );
    
    // compute distinct objects
    identifyDistinctObjects( mask, neighborFindingMethod, objectSize );
    return objectSize->size();
}
Int computeCompactObjects(Array<shortInt,3> &mask,
                          int const neighborFindingMethod,
                          Array<int,3> *outMask,
                          vector<Int> *objectSize)
{
    // compute selection of valid cells (i.e. a value of 1)
    selection( mask, shortInt(1), numeric_limits<shortInt>::max(), outMask, 0, -1 );
    
    // compute distinct objects
    identifyDistinctObjects( outMask, neighborFindingMethod, objectSize );
    return objectSize->size();
}
Int computeCompactObjects(Array<int,3> *mask,
                          int const neighborFindingMethod,
                          vector<Int> *objectSize)
{
    identifyDistinctObjects( mask, neighborFindingMethod, objectSize );
    return objectSize->size();
}


/* Computes the distinct objects, labels them in decreasing order of size and keeps only the most significant ones (the ones above a certain size).
NOTE: the size is given by the volume of each object. */
Int significantObjects(Array<Real,3> &response,
                       Real const threshold,
                       int const neighborFindingMethod,
                       Array<int,3> *mask,
                       vector<Int> *objectSize,
                       Int const minObjectSize)
{
    // get the distinct objects
    computeCompactObjects( response, threshold, neighborFindingMethod, mask, objectSize );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( mask, objectSize, minObjectSize, false );
    return objectSize->size();
}
Int significantObjects(Array<shortInt,3> &mask,
                       int const neighborFindingMethod,
                       Array<int,3> *outMask,
                       vector<Int> *objectSize,
                       Int const minObjectSize)
{
    // get the distinct objects
    computeCompactObjects( mask, neighborFindingMethod, outMask, objectSize );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( outMask, objectSize, minObjectSize, false );
    return objectSize->size();
}
Int significantObjects(Array<int,3> *mask,
                       int const neighborFindingMethod,
                       vector<Int> *objectSize,
                       Int const minObjectSize)
{
    // get the distinct objects
    computeCompactObjects( mask, neighborFindingMethod, objectSize );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( mask, objectSize, minObjectSize, false );
    return objectSize->size();
}


/* Computes the distinct objects, labels them in decreasing order of size and keeps only the most significant ones (the ones above a certain size).
NOTE: the size is given by the mass of each object. */
Int significantObjects(Array<Real,3> &response,
                       Real const threshold,
                       int const neighborFindingMethod,
                       Array<int,3> *mask,
                       Array<Real,3> &mass,
                       vector<double> *objectSize,
                       double const minObjectSize)
{
    // get the distinct objects
    vector<Int> objectVolume;
    computeCompactObjects( response, threshold, neighborFindingMethod, mask, &(objectVolume) );
    
    // compute the mass of each object
    objectsMass( *mask, mass, objectSize, objectVolume.size() );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( mask, objectSize, minObjectSize, false );
    return objectSize->size();
}
Int significantObjects(Array<shortInt,3> &mask,
                       int const neighborFindingMethod,
                       Array<int,3> *outMask,
                       Array<Real,3> &mass,
                       vector<double> *objectSize,
                       double const minObjectSize)
{
    // get the distinct objects
    vector<Int> objectVolume;
    computeCompactObjects( mask, neighborFindingMethod, outMask, &(objectVolume) );
    
    // compute the mass of each object
    objectsMass( *outMask, mass, objectSize, objectVolume.size() );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( outMask, objectSize, minObjectSize, false );
    return objectSize->size();
}
Int significantObjects(Array<int,3> *mask,
                       int const neighborFindingMethod,
                       Array<Real,3> &mass,
                       vector<double> *objectSize,
                       double const minObjectSize)
{
    // get the distinct objects
    vector<Int> objectVolume;
    computeCompactObjects( mask, neighborFindingMethod, &(objectVolume) );
    
    // compute the mass of each object
    objectsMass( *mask, mass, objectSize, objectVolume.size() );
    
    //relabel the objects according to size and discard the small ones
    relabelSignificantObjects( mask, objectSize, minObjectSize, false );
    return objectSize->size();
}



/* Returns the clean response map for the given input parameters: for the given maximum response map with values >= 'threshold' and keeping only the objects with sizes >= 'minSize'.
*/
void cleanResponse(Array<Real,3> &response,
                   Real const threshold,
                   int const neighborFindingMethod,
                   Int const minObjectSize,
                   Array<shortInt,3> *mask,
                   bool const VERBOSE)
{
    if (VERBOSE) cout << "Computing the clean response using optimal threshold =" << threshold << " and discarding objects with volumes smaller than " << minObjectSize << " ... " << flush;
    
    // get the distinct objects
    Array<Int,1> size = response.axisSize<Int>();
    Array<int,3> tempMask( size.ptrData(), size.totalSize() );
    vector<Int> objectVolume;
    computeCompactObjects( response, threshold, neighborFindingMethod, &tempMask, &(objectVolume) );
    
    
    for (Int i=0; i<tempMask.totalSize(); ++i)
        if ( tempMask[i]>=0 and objectVolume[ tempMask[i] ]>=minObjectSize )
            (*mask)[i] = 1;
        else
            (*mask)[i] = 0;
    
    if (VERBOSE)
    {
        int noValidObjects = 0;
        for (size_t i=0; i<objectVolume.size(); ++i)
            if ( objectVolume[i]>=minObjectSize ) ++noValidObjects;
        cout << "Done.\n\tFound " << noValidObjects << " objects with volumes larger than " << minObjectSize << "\n";
    }
}
void cleanResponse(Array<Real,3> &response,
                   Real const threshold,
                   int const neighborFindingMethod,
                   Array<Real,3> &mass,
                   double const minObjectSize,
                   Array<shortInt,3> *mask,
                   bool const VERBOSE)
{
    if (VERBOSE) cout << "Computing the clean response using optimal threshold =" << threshold << " and discarding objects with masses smaller than " << minObjectSize << " ... " << flush;
    
    // get the distinct objects
    Array<Int,1> size = response.axisSize<Int>();
    Array<int,3> tempMask( size.ptrData(), size.totalSize() );
    vector<Int> objectVolume;
    computeCompactObjects( response, threshold, neighborFindingMethod, &tempMask, &(objectVolume) );
    
    // compute the mass of each object
    vector<double> objectMass;
    objectsMass( tempMask, mass, &objectMass, objectVolume.size() );
    
    
    for (Int i=0; i<tempMask.totalSize(); ++i)
        if ( tempMask[i]>=0 and objectMass[ tempMask[i] ]>=minObjectSize )
            (*mask)[i] = 1;
        else
            (*mask)[i] = 0;
    
    if (VERBOSE)
    {
        int noValidObjects = 0;
        for (size_t i=0; i<objectMass.size(); ++i)
            if ( objectMass[i]>=minObjectSize ) ++noValidObjects;
        cout << "Done.\n\tFound " << noValidObjects << " objects with masses larger than " << minObjectSize << "\n";
    }
}

