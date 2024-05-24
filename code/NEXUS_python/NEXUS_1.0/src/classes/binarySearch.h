

#ifndef BINARYSEARCH_HEADER
#define BINARYSEARCH_HEADER


#include <new>
#include <cmath>
#include <vector>

#include <misc.h>
#include <box.h>

/*
 * This class takes as input an array of positions and a separation 'radius'. It than computes the optimal grid to cover the data such that each grid cell has size about 'radius' and one can search for neighbors just by simply looping over all the particles in the given cells and its 26 cell neighbors.
*/

template <typename T>
struct BSR
{
    size_t id;  // stores the id of the point
    T      dis; // stores the square distance of the point
    
    BSR(size_t idx, T disx)
    {
        id = idx;
        dis = disx;
    }
    BSR()
    {
        id = size_t(-1);
        dis = T(-1.);
    }
};




template <typename T, size_t N>
struct BinarySearch
{
    Box<T,N> box;       // keeps the coordinates of the box: xMin, xMax, yMin, yMax,...
    size_t grid[N];     // keeps the grid dimensions along each axis
    T dx[N];            // the grid spacing along each dimension
    T radius;           // the radius of the search - used to compute the grid and the grid spacings
    T R2;               // the square of 'radius'
    
    T *data;            // pointer to the input data and positions
    size_t size;        // the number of elements in the data array
    size_t *indices;    // pointer that keeps track of the particle indices - ordered according to the spatial position
    size_t *offsets;    // pointer to keep track of the offsets of the particles in the 'indices' array for the different grid cells
    int    *noCount;    // pointer to keep track of how many particles are in each cell
    bool dataOn;        // true if the user supplied data
    
    
    
    BinarySearch(T *boxCoords, T const r, bool const periodic=false);   //constructor of the class - if 'periodic=true' than the box coordinates become boxCoords+-radius
    BinarySearch(Box<T,N> &boxCoords, T const r, bool const periodic=false);   
    ~BinarySearch();                                          // class destructor - release the memory allocated in the class
    
    void buildTree(T *inputData, size_t const noElements);    // inserts the input data and builds a search tree using the data
    
    // functions that search for a given particle
    void pointNeighbors_radius(T *pos, std::vector< BSR<T> > *result);                // returns in a vector all the points that are within the distance 'radius' form the target point
    void pointNeighbors_radius(T *pos, T const newR2, std::vector< BSR<T> > *result); // returns in a vector all the points that are within the distance sqrt('newR2') form the target point
    void pointNeighbors_closest(T *pos, int const noPoints, std::vector< BSR<T> > *result); // returns the closest 'noPoints' points within distance 'radius' (if less points were found, it returns size_t(-1) as id and dis=4*R2 )
    
    void massCenter(T *pos, T *result);                       // returns the mass center of the point cloud around that point
    
    
    // functions that search for the full array - faster
    void massCenter(T *result, size_t const sizeResult);      // computes the mass center of the point cloud associated to each input point
    void noNeighbors(int *result, size_t const sizeResult);   // the number of neighbors within radius for each point
    
    
    void inertiaTensor(T *result, size_t const sizeResult );  // computes the reduced inertia tensor of the neighbors for each particle
    
    
    //functions used by the class
    Box<int,N> gridIndicesBox;                          // stores (0,grid[0]-1,0,grid[1]-1,..) -used to check if the point is inside the grid
    int        noMaxNeighbors;                          // maximum number of neighbor cells
    size_t     gridSize;                                // the total number of cells in the grid
    
    inline size_t getPointGridIndex(T *pos);            // returns the grind index corresponding to cell in which the point lies
    inline T distance(T *pos1, T*pos2);                 // returns the square of the distance between the two positions
    // returns the indices for all the cell's neighbors (the return value is the number of neighbors) -> the next 4 functions
    inline int neighborCells(size_t cellIndex, size_t *neighbors);  
    inline int neighborCells(int n0, size_t *neighbors);
    inline int neighborCells(int n0, int n1, size_t *neighbors);
    inline int neighborCells(int n0, int n1, int n2, size_t *neighbors);
    void construct(T *boxCoords, T const r, bool const periodic);   //implements the actions of the constructor
};




/* General constructor for the class. */
template <typename T, size_t N> BinarySearch<T,N>::BinarySearch(T *boxCoords, T const r, bool const periodic)
{
    this->construct( boxCoords, r, periodic );
}
template <typename T, size_t N> BinarySearch<T,N>::BinarySearch(Box<T,N> &boxCoords, T const r, bool const periodic)
{
    this->construct( boxCoords.coords, r, periodic );
}

template <typename T, size_t N> BinarySearch<T,N>::~BinarySearch()
{
    if (dataOn)
    {
        delete[] indices;
        delete[] offsets;
        delete[] noCount;
    }
}


/* Insert the data and build the tree. */
template <typename T, size_t N> void BinarySearch<T,N>::buildTree(T *inputData, size_t const noElements)
{
    dataOn = true;
    size   = noElements;
    data   = inputData;
    
    // allocate memory and check if allocation was succesful
    size_t *tempIndex;
    indices = new (std::nothrow) size_t[size];
    tempIndex = new (std::nothrow) size_t[size];
    offsets = new (std::nothrow) size_t[gridSize];
    noCount = new (std::nothrow) int[gridSize];
    memoryAllocationCheck( indices, "Error while allocating memory for the 'indices' variables inside the function 'BinarySearch::buildTree'." );
    memoryAllocationCheck( tempIndex, "Error while allocating memory for the 'tempIndex' variables inside the function 'BinarySearch::buildTree'." );
    memoryAllocationCheck( offsets, "Error while allocating memory for the 'offsets' variables inside the function 'BinarySearch::buildTree'." );
    memoryAllocationCheck( noCount, "Error while allocating memory for the 'noCount' variables inside the function 'BinarySearch::buildTree'." );
    for (int i=0; i<gridSize; ++i)
        noCount[i] = 0;
    
    
    // find the grid cell associated to each particle
    for (size_t i=0; i<size; ++i)
    {
        tempIndex[i] = getPointGridIndex( &(inputData[i*N]) );
        if ( tempIndex[i]>=gridSize ) continue;
        ++noCount[ tempIndex[i] ];
    }
    
    
    // compute the offsets in the 'indices' array of the particles in each grid cell
    offsets[0] = 0;
    for (size_t i=1; i<gridSize; ++i)
        offsets[i] = offsets[i-1] + noCount[i-1];
    
    
    // write in the 'indices' array the point indices corresponding to each cell, in increasing order of the grid cell index
    for (size_t i=0; i<gridSize; ++i)
        noCount[i] = 0;
    for (size_t i=0; i<size; ++i)
        if ( tempIndex[i]<gridSize )
        {
            size_t ind = tempIndex[i];
            indices[ offsets[ind]+noCount[ind] ] = i;
            ++noCount[ind];
        }
    
    delete[] tempIndex;
}




//! Functions that search for the neighbors of a single point within a given distance

// returns all the neighbors within distance R from the point 'pos'
template <typename T, size_t N> void BinarySearch<T,N>::pointNeighbors_radius(T *pos, std::vector< BSR<T> > *result)
{
    // get the point position
    size_t index = getPointGridIndex( pos );
    if ( index>=gridSize ) return;
    
    // get the neighbors associated to the cell in which the point lies
    size_t neighbors[noMaxNeighbors];
    int noCellNeighbors = this->neighborCells( index, neighbors );
    
    // loop over all the cells and find the neighbors within the distance
    result->clear();
    for (int i=0; i<noCellNeighbors; ++i)
    {
        size_t cellIndex = neighbors[i];
        for (size_t j=offsets[cellIndex]; j<offsets[cellIndex]+noCount[cellIndex]; ++j )
        {
            T *pos1 = &( data[ indices[j]*N ] );
            T dis = this->distance( pos, pos1 );
            if ( dis>R2 )
                continue;
            else
                result->push_back( BSR<T>(indices[j],dis) );
        }
    }
}

/* Returns all the neighbors within distance sqrt( newR2 ) from the point 'pos'.
NOTE: You need to have newR2 < R2.
*/
template <typename T, size_t N> void BinarySearch<T,N>::pointNeighbors_radius(T *pos, T const newR2, std::vector< BSR<T> > *result)
{
    if ( (newR2-R2)/R2 > 1.e-4 )
        throwError( "You are trying to find the neighbors within a give distance, but this distance is larger than the grid size used in the binary tree." );
    // get the point position
    size_t index = getPointGridIndex( pos );
    if ( index>=gridSize ) return;
    
    // get the neighbors associated to the cell in which the point lies
    size_t neighbors[noMaxNeighbors];
    int noCellNeighbors = this->neighborCells( index, neighbors );
    
    // loop over all the cells and find the neighbors within the distance
    result->clear();
    for (int i=0; i<noCellNeighbors; ++i)
    {
        size_t cellIndex = neighbors[i];
        for (size_t j=offsets[cellIndex]; j<offsets[cellIndex]+noCount[cellIndex]; ++j )
        {
            T *pos1 = &( data[ indices[j]*N ] );
            T dis = this->distance( pos, pos1 );
            if ( dis>newR2 )
                continue;
            else
                result->push_back( BSR<T>(indices[j],dis) );
        }
    }
}




//! Functions that return the n-closest neighbors of a point up to distance R
// Returns the closest 'noPoints' from 'pos'
template <typename T, size_t N> void BinarySearch<T,N>::pointNeighbors_closest(T *pos, int const noPoints, std::vector< BSR<T> > *result)
{
    // get the point position
    size_t index = getPointGridIndex( pos );
    if ( index>=gridSize ) return;
    
    // get the neighbors associated to the cell in which the point lies
    size_t neighbors[noMaxNeighbors];
    int noCellNeighbors = this->neighborCells( index, neighbors );
    
    // loop over all the cells and find the closest 'noPoints' neighbors within the distance R2
    result->assign( noPoints, BSR<T>(size_t(-1),4*R2) );   // assign invalid neighbors
    for (int i=0; i<noCellNeighbors; ++i)
    {
        size_t cellIndex = neighbors[i];
        for (size_t j=offsets[cellIndex]; j<offsets[cellIndex]+noCount[cellIndex]; ++j )
        {
            T *pos1 = &( data[ indices[j]*N ] );
            T dis = this->distance( pos, pos1 );
            if ( dis>(*result)[noPoints-1].dis )
                continue;
            
            // add this new point
            int k;
            BSR<T> temp[2];
            for ( k=0; k<noPoints; ++k)    //insert the new point
                if ( dis<(*result)[k].dis )
                {
                    temp[k%2] = (*result)[k];
                    (*result)[k] = BSR<T>(indices[j],dis);
                    break;
                }
            for ( k+=1; k<noPoints; ++k)    //resuffle the other points
            {
                temp[k%2] = (*result)[k];
                (*result)[k] = temp[(k-1)%2];
            }
        }
    }
}







/* Returns the grid index for the given point. */
template <typename T, size_t N> size_t BinarySearch<T,N>::getPointGridIndex(T *pos)
{
    // get the grid indices of the point
    int ind[N];
    for (size_t i=0; i<N; ++i)
        ind[i] = int( std::floor( (pos[i]-box[2*i])/dx[i] ) );
    
    // check if point is inside the grid
    if ( not gridIndicesBox.isPointInBox(ind) )
        return gridSize*10+1;
    
    // if point inside the grid, return grid index
    size_t index = ind[0];
    for (size_t i=1; i<N; ++i)
        index = index*grid[i] + ind[i];
    return index;
}

/* Returns the square distance between two points */
template <typename T, size_t N> T BinarySearch<T,N>::distance(T *pos1, T *pos2)
{
    T dis = T(0.);
    // get the grid indices of the point
    for (size_t i=0; i<N; ++i)
    {
        T temp = pos1[i] - pos2[i];
        dis += temp*temp;
    }
    return dis;
}

/* Returns all the cell neighbors for the given cell */
template <typename T, size_t N> int BinarySearch<T,N>::neighborCells(size_t cellIndex, size_t *neighbors)
{
    if (N==3)
    {
        int n0 = cellIndex / (grid[1]*grid[2]);
        size_t temp = cellIndex - n0*grid[1]*grid[2];
        int n1 = temp / grid[2];
        int n2 = temp % grid[2];
        return neighborCells( n0, n1, n2, neighbors );
    }
    else if (N==2)
    {
        int n0 = cellIndex / grid[1];
        int n1 = cellIndex % grid[1];
        return neighborCells( n0, n1, neighbors );
    }
    else if (N==1)
        return neighborCells( int(cellIndex), neighbors );
}
template <typename T, size_t N> int BinarySearch<T,N>::neighborCells(int n0, size_t *neighbors)
{
    int count = 0;
    for (int i0=n0-1; i0<=n0+1; ++i0)
        if ( i0>=gridIndicesBox[0] and i0<=gridIndicesBox[1] )
            neighbors[count++] = i0;
    return count;
}
template <typename T, size_t N> int BinarySearch<T,N>::neighborCells(int n0, int n1, size_t *neighbors)
{
    int count = 0;
    for (int i0=n0-1; i0<=n0+1; ++i0)
        for (int i1=n1-1; i1<=n1+1; ++i1)
            if ( i0>=gridIndicesBox[0] and i0<=gridIndicesBox[1] and i1>=gridIndicesBox[2] and i1<=gridIndicesBox[3] )
                neighbors[count++] = i0*grid[1] + i1;
    return count;
}
template <typename T, size_t N> int BinarySearch<T,N>::neighborCells(int n0, int n1, int n2, size_t *neighbors)
{
    int count = 0;
    for (int i0=n0-1; i0<=n0+1; ++i0)
        for (int i1=n1-1; i1<=n1+1; ++i1)
            for (int i2=n2-1; i2<=n2+1; ++i2)
                if ( i0>=gridIndicesBox[0] and i0<=gridIndicesBox[1] and i1>=gridIndicesBox[2] and i1<=gridIndicesBox[3] and i2>=gridIndicesBox[4] and i2<=gridIndicesBox[5] )
                    neighbors[count++] = (i0*grid[1] + i1)*grid[2] + i2;
    return count;
}


// general constructor for the class
template <typename T, size_t N> void BinarySearch<T,N>::construct(T *boxCoords, T const r, bool const periodic)
{
    // copy the data
    dataOn = false;
    radius = r;
    R2 = radius*radius;
    box.assign( boxCoords, 2*N );
    if ( periodic ) box.addPadding( radius );
    
    //compute the grid dimensions
    for (size_t i=0; i<N; ++i)
    {
        grid[i] = int( std::floor( (box[2*i+1]-box[2*i])/radius ) );
        dx[i] = (box[2*i+1]-box[2*i]) / grid[i];    //always have dx[i] >= radius
        gridIndicesBox[2*i+0] = 0;
        gridIndicesBox[2*i+1] = int(grid[i]-1);
    }
    
    if ( N>=4 ) throwError( "The class 'BinarySearch' is not fully implemented for 4 or more dimensions." );
    noMaxNeighbors = N==1 ? 3 : (N==2 ? 9 : 27);
    gridSize = N==1 ? grid[0] : (N==2 ? grid[0]*grid[1] : grid[0]*grid[1]*grid[2]);
}



#endif
