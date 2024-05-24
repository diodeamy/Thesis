/*

This class is used to get a variable with the behavior of a k-space momentum for a FFTW k-space array.

*/


#ifndef K_SPACE_HEADER
#define K_SPACE_HEADER

#include <iostream>
#include <miscellaneous.h>
#include <boost/assert.hpp>
#include <cmath>



template <typename T, int N>
struct K_space
{
    K_space(int const sizeX, T const boxSize);   // constructors
    K_space(int const sizeX, int const sizeY, T const boxSizeX, T const boxSizeY);
    K_space(int const sizeX, int const sizeY, int const sizeZ, T const boxSizeX, T const boxSizeY, T const boxSizeZ);
    K_space(int const size1, int const size2, int const size3, int const size4, T const boxSize1, T const boxSize2, T const boxSize3, T const boxSize4);
    K_space(int const *gridSize, int const noDim, T const *boxSize);
    K_space(size_t const *gridSize, size_t const noDim, T const *boxSize);
    
    
    int totalGridSize() const;                 // returns total number of points of the grid
    int getGridSize(int const i) const;        // returns grid size along 'i'-th direction
    int* getGridSize() const;                  // returns pointer to grid size
    int dimension() const;                     // returns the dimension 'N' of the space
    T boxSize(int const i) const;              // returns length of simulation box along axis i
    
    T n(int const n1) const;                   // return 'n=sqrt(sum_i n_i^2)' with n_i in the interval -grid/2+1 to grid/2 - this is the momentum
    T n(int const n1, int const n2) const;
    T n(int const n1, int const n2, int const n3) const;
    T n(int const n1, int const n2, int const n3, int const n4) const;
    T n(int const *n1, int const noN) const;   // returns again the momentum, but the input is an array of k-grid points of dimensions noN==N (checks for this condition)
    
    T nSquare(int const n1) const;             // return 'n=sum_i n_i^2' with n_i in the interval -grid/2+1 to grid/2 - this is the momentum square
    T nSquare(int const n1, int const n2) const;
    T nSquare(int const n1, int const n2, int const n3) const;
    T nSquare(int const n1, int const n2, int const n3, int const n4) const;
    T nSquare(int const *n1, int const noN) const;// returns again the momentum square, but the input is an array of k-grid points of dimensions noN==N (checks for this condition)
    
    int nx(int const i) const;                 // returns 'i' translated to -grid/2+1 to grid/2 for direction 'x'
    int ny(int const i) const;
    int nz(int const i) const;
    int nn(int const i, int const direction) const;// returns 'i' translated to -grid/2+1 to grid/2 for direction 'direction'
    
    
    T k(int const n1) const;                   // return 'k=sqrt(sum_i n_i^2)*scaling' with n_i in the interval -grid/2+1 to grid/2 - this is the momentum
    T k(int const n1, int const n2) const;
    T k(int const n1, int const n2, int const n3) const;
    T k(int const n1, int const n2, int const n3, int const n4) const;
    T k(int const *n1, int const noN) const;   // returns again the momentum, but the input is an array of k-grid points of dimensions noN==N (checks for this condition)
    
    T kSquare(int const n1) const;             // return 'n=(sum_i n_i^2)*scaling^2' with n_i in the interval -grid/2+1 to grid/2 - this is the momentum square
    T kSquare(int const n1, int const n2) const;
    T kSquare(int const n1, int const n2, int const n3) const;
    T kSquare(int const n1, int const n2, int const n3, int const n4) const;
    T kSquare(int const *n1, int const noN) const;// returns again the momentum square, but the input is an array of k-grid points of dimensions noN==N (checks for this condition)
    
    T kx(int const i) const;                   // returns 'i*scaling' translated to -grid/2+1 to grid/2 for direction 'x'
    T ky(int const i) const;
    T kz(int const i) const;
    T kn(int const i, int const direction) const;// returns 'i*scaling' translated to -grid/2+1 to grid/2 for direction 'direction'
    
    
    
    private:
    int grid[N];
    T boxLength[N];
    T dx[N];
    T scaling[N];
    void setGrid(int const *gridSize);         //sets the values for 'grid' from pointer 'gridSize'
    void setScale(T const *boxSize);           //sets 'boxLength' and 'dx'
};

/* Class constructors. */
template <typename T, int N> K_space<T,N>::K_space(int const sizeX, T const boxSize)
{
    if ( N!=1 ) throwError( "You called the 'K_space' class constructor for a 1 dimensional space but your 'K_space' class has another dimension." );
    
    int size[] = {sizeX};
    this->setGrid( size );
    T _boxSize[] = {boxSize};
    this->setScale( _boxSize );
}
template <typename T, int N> K_space<T,N>::K_space(int const sizeX, int const sizeY, T const boxSizeX, T const boxSizeY )
{
    if ( N!=2 ) throwError( "You called the 'K_space' class constructor for a 2 dimensional space but your 'K_space' class has another dimension." );
    
    int size[] = {sizeX,sizeY};
    this->setGrid( size ); 
    T _boxSize[] = {boxSizeX, boxSizeY};
    this->setScale( _boxSize );
}
template <typename T, int N> K_space<T,N>::K_space(int const sizeX, int const sizeY, int const sizeZ, T const boxSizeX,  T const boxSizeY,  T const boxSizeZ)
{
    if ( N!=3 ) throwError( "You called the 'K_space' class constructor for a 3 dimensional space but your 'K_space' class has another dimension." );
    
    int size[] = {sizeX,sizeY,sizeZ};
    this->setGrid( size ); 
    T _boxSize[] = {boxSizeX, boxSizeY, boxSizeZ};
    this->setScale( _boxSize );
}
template <typename T, int N> K_space<T,N>::K_space(int const size1, int const size2, int const size3, int const size4, T const boxSize1,  T const boxSize2,  T const boxSize3, T const boxSize4)
{
    if ( N!=4 ) throwError( "You called the 'K_space' class constructor for a 4 dimensional space but your 'K_space' class has another dimension." );
    
    int size[] = {size1,size2,size3,size4};
    this->setGrid( size );
    T _boxSize[] = {boxSize1, boxSize2, boxSize3, boxSize4 };
    this->setScale( _boxSize );
}
template <typename T, int N> K_space<T,N>::K_space(int const *gridSize, int const noDim, T const *boxSize)
{
    if ( N!=noDim ) throwError( "You called the 'K_space' class constructor with a different dimension than the one assigned to the 'K_space' class element." );
    
    this->setGrid( gridSize );
    this->setScale( boxSize );
}
template <typename T, int N> K_space<T,N>::K_space(size_t const *gridSize, size_t const noDim, T const *boxSize)
{
    if ( N!=noDim ) throwError( "You called the 'K_space' class constructor with a different dimension than the one assigned to the 'K_space' class element." );
    
    for (int i=0; i<N; ++i)
        grid[i] = gridSize[i];
    this->setScale( boxSize );
}



//! Function to get access to the class variables
/* Returns total number of points of the grid. */
template <typename T, int N> int K_space<T,N>::totalGridSize() const
{
    int temp = 1;
    for (int i=0; i<N; ++i)
        temp *= grid[i];
    return temp;
}
/* Returns the grid size along 'i'-th direction. */
template <typename T, int N> int K_space<T,N>::getGridSize(int const i) const
{
    BOOST_ASSERT( i<N );
    return grid[i];
}
/* Returns a pointer to the grid size array. */
template <typename T, int N> int* K_space<T,N>::getGridSize() const
{
    return grid;
}
/* Returns the dimension of the space. */
template <typename T, int N> int K_space<T,N>::dimension() const
{
    return N;
}
/* Returns the box length along axis i. */
template <typename T, int N> T K_space<T,N>::boxSize(int const i) const
{
    BOOST_ASSERT( i<N );
    return boxLength[i];
}






//! Functions that compute 'n', 'n'^2 and 'n'_i
/* Functions that compute 'n' = sqrt( sum_i n_i^2 ) */
template <typename T, int N> inline T K_space<T,N>::n(int const n1) const
{
    return sqrt( this->nSquare(n1) );
}
template <typename T, int N> inline T K_space<T,N>::n(int const n1, int const n2) const
{
    return sqrt( this->nSquare(n1,n2) );
}
template <typename T, int N> inline T K_space<T,N>::n(int const n1, int const n2, int const n3) const
{
    return sqrt( this->nSquare(n1,n2,n3) );
}
template <typename T, int N> inline T K_space<T,N>::n(int const n1, int const n2, int const n3, int const n4) const
{
    return sqrt( this->nSquare(n1,n2,n3,n4) );
}
template <typename T, int N> inline T K_space<T,N>::n(int const *n1, int const noN) const
{
    return sqrt( this->nSquare(n1,noN) );
}

/* Functions that compute 'n^2' = sum_i n_i^2 */
template <typename T, int N> inline T K_space<T,N>::nSquare(int const n1) const
{
    BOOST_ASSERT( N==1 );
    int temp = nx(n1)*nx(n1);
    return T( temp );
}
template <typename T, int N> inline T K_space<T,N>::nSquare(int const n1, int const n2) const
{
    BOOST_ASSERT( N==2 );
    BOOST_ASSERT( grid[0]==grid[1] );
    int temp = nx(n1)*nx(n1) + ny(n2)*ny(n2);
    return T( temp );
}
template <typename T, int N> inline T K_space<T,N>::nSquare(int const n1, int const n2, int const n3) const
{
    BOOST_ASSERT( N==3 );
    BOOST_ASSERT( grid[0]==grid[1] and grid[0]==grid[2] );
    int temp = nx(n1)*nx(n1) + ny(n2)*ny(n2) + nz(n3)*nz(n3);
    return T( temp );
}
template <typename T, int N> inline T K_space<T,N>::nSquare(int const n1, int const n2, int const n3, int const n4) const
{
    BOOST_ASSERT( N==4 );
    BOOST_ASSERT( grid[0]==grid[1] and grid[0]==grid[2] and grid[0]==grid[3] );
    int temp = nx(n1)*nx(n1) + ny(n2)*ny(n2) + nz(n3)*nz(n3) + nn(n4,3)*nx(n4,3);
    return T( temp );
}
template <typename T, int N> inline T K_space<T,N>::nSquare(int const *n1, int const noN) const
{
    BOOST_ASSERT( N==noN );
    for (int i=1; i<N; ++i)
        BOOST_ASSERT( grid[0]==grid[i] );
    
    int temp = 0;
    for (int i=0; i<N; ++i)
        temp += nn(n1[i],i) * nn(n1[i],i);
    return T( temp );
}

/* Functions that compute n_i. */
template <typename T, int N> inline int K_space<T,N>::nx(int const i) const
{
    BOOST_ASSERT( N>=1 );
    if (i>grid[0]/2)
        return i-grid[0];
    else
        return i;
}
template <typename T, int N> inline int K_space<T,N>::ny(int const i) const
{
    BOOST_ASSERT( N>=2 );
    if (i>grid[1]/2)
        return i-grid[1];
    else
        return i;
}
template <typename T, int N> inline int K_space<T,N>::nz(int const i) const
{
    BOOST_ASSERT( N>=3 );
    if (i>grid[2]/2)
        return i-grid[2];
    else
        return i;
}
template <typename T, int N> inline int K_space<T,N>::nn(int const i, int const direction) const
{
    BOOST_ASSERT( N>=direction+1 );
    if (i>grid[direction]/2)
        return i-grid[direction];
    else
        return i;
}





//! Functions that compute 'k', 'k'^2 and 'k'_i
/* Functions that compute 'k' = 'scaling' * 'n' */
template <typename T, int N> inline T K_space<T,N>::k(int const n1) const
{
    return sqrt( this->kSquare(n1) );
}
template <typename T, int N> inline T K_space<T,N>::k(int const n1, int const n2) const
{
    return sqrt( this->kSquare(n1,n2) );
}
template <typename T, int N> inline T K_space<T,N>::k(int const n1, int const n2, int const n3) const
{
    return sqrt( this->kSquare(n1,n2,n3) );
}
template <typename T, int N> inline T K_space<T,N>::k(int const n1, int const n2, int const n3, int const n4) const
{
    return sqrt( this->kSquare(n1,n2,n3,n4) );
}
template <typename T, int N> inline T K_space<T,N>::k(int const *n1, int const noN) const
{
    return sqrt( this->kSquare(n1,noN) );
}

/* Functions that compute 'k'^2 = 'scaling'^2 * 'n'^2 */
template <typename T, int N> inline T K_space<T,N>::kSquare(int const n1) const
{
    BOOST_ASSERT( N==1 );
    T tempX = 2./dx[0] * std::sin( scaling[0]/2.*nx(n1) );
    return tempX*tempX;
}
template <typename T, int N> inline T K_space<T,N>::kSquare(int const n1, int const n2) const
{
    BOOST_ASSERT( N==2 );
    T tempX = 2./dx[0] * std::sin( scaling[0]/2.*nx(n1) );
    T tempY = 2./dx[1] * std::sin( scaling[1]/2.*ny(n2) );
    return tempX*tempX + tempY*tempY;
}
template <typename T, int N> inline T K_space<T,N>::kSquare(int const n1, int const n2, int const n3) const
{
    BOOST_ASSERT( N==3 );
    T tempX = 2./dx[0] * std::sin( scaling[0]/2.*nx(n1) );
    T tempY = 2./dx[1] * std::sin( scaling[1]/2.*ny(n2) );
    T tempZ = 2./dx[2] * std::sin( scaling[2]/2.*nz(n3) );
    return tempX*tempX + tempY*tempY + tempZ*tempZ;
}
template <typename T, int N> inline T K_space<T,N>::kSquare(int const n1, int const n2, int const n3, int const n4) const
{
    BOOST_ASSERT( N==4 );
    T tempX = 2./dx[0] * std::sin( scaling[0]/2.*nx(n1) );
    T tempY = 2./dx[1] * std::sin( scaling[1]/2.*ny(n2) );
    T tempZ = 2./dx[2] * std::sin( scaling[2]/2.*nz(n3) );
    T temp4 = 2./dx[3] * std::sin( scaling[3]/2.*nn(n4,3) );
    return tempX*tempX + tempY*tempY + tempZ*tempZ + temp4*temp4;
}
template <typename T, int N> inline T K_space<T,N>::kSquare(int const *n1, int const noN) const
{
    BOOST_ASSERT( N==noN );
    T temp = 0.;
    for (int i=0; i<noN; ++i)
    {
        T rest = 2./dx[i] * std::sin( scaling[i]/2.*nn(n1[i],i) );
        temp += rest * rest;
    }
    return temp;
}

/* Functions that compute 'k'_i = 'scaling' * 'n'_i */
template <typename T, int N> inline T K_space<T,N>::kx(int const i) const
{
    return std::sin(scaling[0]*nx(i)) / dx[0];
}
template <typename T, int N> inline T K_space<T,N>::ky(int const i) const
{
    return std::sin(scaling[1]*ny(i)) / dx[1];
}
template <typename T, int N> inline T K_space<T,N>::kz(int const i) const
{
    return std::sin(scaling[2]*nz(i)) / dx[2];
}
template <typename T, int N> inline T K_space<T,N>::kn(int const i, int const direction) const
{
    return std::sin(scaling[direction]*nn(i,direction)) / dx[direction];;
}



//! Some private functions
/* Set the values of 'grid' variable. */
template <typename T, int N> inline void K_space<T,N>::setGrid(int const *gridSize)
{
    for (int i=0; i<N; ++i)
        grid[i] = gridSize[i];
}
/* Sets the values of 'boxLength' and 'scaling'. */
template <typename T, int N> inline void K_space<T,N>::setScale(T const *boxSize)
{
    for (int i=0; i<N; ++i)
    {
        boxLength[i] = boxSize[i];
        dx[i] = boxLength[i] / grid[i];
        scaling[i] = 2. * PI / grid[i];
    }
}



#endif
