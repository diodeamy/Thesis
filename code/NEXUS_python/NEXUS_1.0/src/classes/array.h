/*

This is a general as possible class for storing multidimensional arrays and performing operations with the arrays. An important factor when designing the class was computational speed (generality may be droped to increase speed for very common cases).

There are multiple ways of accessing the array elements. For example for a 3D array, several possibilities are:
      A.at(global index); A.at( i,j,k );
      A[global index];  A(global index);   A(i,j,k);
      A[global index] - doesn't do a range check

There are also many operations implemented by default between the class and similar class as well as for T elements (e.g. sum the entries of two clases, multiply the entries of one class by a constant).

NOTE:
     -1) IT DOES NOT PERFORM MATRIX OPERATIONS (with the exception of addition and substraction).
     0) The class performs index and error checks. The most computationally demanding ones can be dissabled using the compiler directive BOOST_DISABLE_ASSERTS (since it uses BOOST_ASSERT). Use BOOST_ENABLE_ASSERT_HANDLER to get additional information about asserts.
     1) There is a very small overhead for storing the arrays, '2*N+1' intergers which can become important for many small arrays.
     2) The class was not tested with 4D arrays or higher.

*/


#ifndef ARRAY_HEADER
#define ARRAY_HEADER

#include <iostream>
#include <new>
#include <miscellaneous.h>
#include <boost/assert.hpp>




template <typename T, size_t N>
struct Array
{
    T *data;        // pointer to the data
    
    Array(size_t const *dimSize, int const noDimensions);
    Array(int const *dimSize, int const noDimensions);
    Array(Array const &other);
    Array(size_t const dimSize1);
    Array(size_t const dimSize1, size_t const dimSize2 );
    Array(size_t const dimSize1, size_t const dimSize2, size_t const dimSize3 );
    Array(size_t const dimSize1, size_t const dimSize2, size_t const dimSize3, size_t const dimSize4 );
    // keeps track of outside data - it does not delete or free the data at exit
    Array(T *otherData, size_t const *dimSize, int const noDimensions);
    Array(T *otherData, int const *dimSize, int const noDimensions);
    Array(T *otherData, size_t const dimSize1);
    Array(T *otherData, size_t const dimSize1, size_t const dimSize2 );
    Array(T *otherData, size_t const dimSize1, size_t const dimSize2, size_t const dimSize3 );
    Array(T *otherData, size_t const dimSize1, size_t const dimSize2, size_t const dimSize3, size_t const dimSize4 );
    ~Array();
    
    size_t dimensions() const;              // returns the number of dimensions of the array
    size_t size() const;                    // get total size of the array
    size_t totalSize() const;               // get total size of the array
    template<typename SIZE_T> void axisSize(SIZE_T *result, SIZE_T const noElements) const; // returns the array size along each dimension in the 'result' posize_ter
    size_t axisSize(size_t const i) const;     // returns the array size along the given dimension
    template<typename SIZE_T> Array<SIZE_T,1> axisSize() const; // returns the size along each dimension
    template<typename SIZE_T> void getSize(SIZE_T *result, SIZE_T const noElements) const; // returns the array size along each dimension in the 'result' posize_ter
    size_t getSize(size_t const i) const;   // returns the array size along the given dimension
    template<typename SIZE_T> Array<SIZE_T,1> getSize() const;  // returns the size along each dimension
    void print();                           // prsize_ts the array values to the standard output
    void freeMemory();                      // free memory occupied by array elements
    bool similar(Array const &newArray) const;  // check if the two arrays have the same sizes along each dimension
    
    void assign(T const value);             // assign this constant value to all the entries
    void assign(T const *values, size_t const noElements);  // copy the values pointed by the pointer to this class data (if both have the same number of elements)
    void assign(Array const &toCopy); 
    T* ptrData() const;                     // returns pointer to first entry in the 'data'
    T at(size_t const i) const;             // returns the 'i'-th entry of 'data' after performing a check that the index 'i' is within the allowed value of entry indexes (DOES NOT RETURN A REFERENCE)
    T at(size_t const i, size_t const j) const; // returns the (i,j) entry of a 2D array after performing a check that the indeces are within the allowed values of entry indexes (DOES NOT RETURN A REFERENCE)
    T at(size_t const i, size_t const j, size_t const k) const;
    T at(size_t const i, size_t const j, size_t const k, size_t const l) const;
    
    Array& operator =(Array const &right);  // define the copy constructor
    Array& operator =(T const right);       // assign a constant value
    
    bool operator ==(Array const &right) const; // return true if the two arrays are the same ( checks each entry individually )
    bool operator !=(Array const &right) const; // return true if the two arrays are not the same ( checks each entry individually )
    Array  operator +(Array const &right);    // adds the elements of the two arrays ( each entry individually )
    Array& operator +=(Array const &right);   // adds the elements of another array to the present one
    Array  operator +(T const right);         // adds to each array entry a value of type 'T'
    Array& operator +=(T const right);
    Array  operator -(Array const &right);    // substracts the elements of the two arrays ( each entry individually )
    Array& operator -=(Array const &right);
    Array  operator -(T const right);         // substracts from each array entry a value of type 'T'
    Array& operator -=(T const right);
    Array  operator *(Array const &right);    // computes the product of the two arrays (NOT MATRIX PRODUCT) - entry by entry
    Array& operator *=(Array const &right);   // multiplies the entries with the entries of a second array (NOT MATRIX PRODUCT) - entry by entry
    Array  operator *(T const right);         // multiplies all entries by the same element of type 'T'
    Array& operator *=(T const right);
    Array  operator /(Array const &right);    // computes the ratio of the two arrays (NOT MATRIX DIVISION) - entry by entry
    Array& operator /=(Array const &right);   // divides the entries by those of a second array (NOT MATRIX DIVISION) - entry by entry
    Array  operator /(T const right);         // divides all entries by the same element of type 'T'
    Array& operator /=(T const right); 
    
    T& operator [](size_t const i);              // returns a reference to the 'i'-th entry of 'data'
    T& operator ()(size_t const i);              // returns a reference to the 'i'-th entry of 'data'
    T& operator ()(size_t const i, size_t const j); // returns a reference to the (i,j) entry of a 2D array
    T& operator ()(size_t const i, size_t const j, size_t const k); // returns a reference to the (i,j,k) entry of a 3D array
    T& operator ()(size_t const i, size_t const j, size_t const k, size_t const l); // returns a reference to the (i,j,k,l) entry of a 4D array
    
//! define 'iterator'
    
    private:
        size_t dim[N];          // the size of each dimension of N
        size_t indexOffset[N];  // keeps track of 'n[i]' such that for example the position in a 3D array is given by i*n[0] + j*n[1] + k*n[2] where (i,j,k) are the indices of a 3D array
        size_t totSize;         // total size of the array
        bool free;              // true if memory was freed using the 'free()' function
        bool ownData;           // true if data was allocatted by this class
        
        template<typename SIZE_T> void setDimensionSize(SIZE_T const *dimSize, SIZE_T const noDimensions );  // size_ternally sets the dimension sizes
        void computeTotalSize();                                        // computes total number of elements in the array
        void computeIndexOffset();                                      // computes total number of elements in the array
        T dataValue(size_t const i) const;                              // outputs read-only values of 'data[i]'
};


/* General constructor for the class. */
template <typename T, size_t N> Array<T,N>::Array(size_t const *dimSize, int const noElements)
{
    this->setDimensionSize<size_t>( dimSize, noElements );
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(int const *dimSize, int const noElements)
{
    this->setDimensionSize<int>( dimSize, noElements );
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(Array<T,N> const &other)
{
    for (size_t i=0; i<N; ++i)
        dim[i] = other.getSize(i);
    this->computeTotalSize();
    this->computeIndexOffset();
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
    this->assign( other );
}
template <typename T, size_t N> Array<T,N>::Array(size_t const dimSize1)
{
    size_t dimSize[1] = { dimSize1 };
    this->setDimensionSize<size_t>( dimSize, size_t(1) );
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(size_t const dimSize1, size_t const dimSize2)
{
    size_t dimSize[2] = { dimSize1, dimSize2 };
    this->setDimensionSize<size_t>( dimSize, size_t(2) );
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(size_t const dimSize1, size_t const dimSize2, size_t const dimSize3 )
{
    size_t dimSize[3] = { dimSize1, dimSize2, dimSize3 };
    this->setDimensionSize<size_t>( dimSize, size_t(3) );
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(size_t const dimSize1, size_t const dimSize2, size_t const dimSize3, size_t const dimSize4 )
{
    size_t dimSize[4] = { dimSize1, dimSize2, dimSize3, dimSize4 };
    this->setDimensionSize<size_t>( dimSize, size_t(4) );
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );  //check if allocation was succesful
    free = false;
    ownData = true;
}

template <typename T, size_t N> Array<T,N>::Array(T *otherData, size_t const *dimSize, int const noElements)
{
    this->setDimensionSize<size_t>( dimSize, noElements );
    
    data = otherData;
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(T *otherData, int const *dimSize, int const noElements)
{
    this->setDimensionSize<int>( dimSize, noElements );
    
    data = otherData;
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(T *otherData, size_t const dimSize1)
{
    size_t dimSize[1] = { dimSize1 };
    this->setDimensionSize<size_t>( dimSize, size_t(1) );
    
    data = otherData;
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(T *otherData, size_t const dimSize1, size_t const dimSize2)
{
    size_t dimSize[2] = { dimSize1, dimSize2 };
    this->setDimensionSize<size_t>( dimSize, size_t(2) );
    
    data = otherData;
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(T *otherData, size_t const dimSize1, size_t const dimSize2, size_t const dimSize3 )
{
    size_t dimSize[3] = { dimSize1, dimSize2, dimSize3 };
    this->setDimensionSize<size_t>( dimSize, size_t(3) );
    
    data = otherData;
    free = false;
    ownData = true;
}
template <typename T, size_t N> Array<T,N>::Array(T *otherData, size_t const dimSize1, size_t const dimSize2, size_t const dimSize3, size_t const dimSize4 )
{
    size_t dimSize[4] = { dimSize1, dimSize2, dimSize3, dimSize4 };
    this->setDimensionSize<size_t>( dimSize, size_t(4) );
    
    data = otherData;
    free = false;
    ownData = true;
}

/* Destructor for the class. */
template <typename T, size_t N> Array<T,N>::~Array()
{
    if (not free and ownData) delete[] data;
}




//! Functions that return array characteristics
/* Return number of array dimensions. */
template <typename T, size_t N> size_t Array<T,N>::dimensions() const
{
    return N;
}

/* Return number of array dimensions. */
template <typename T, size_t N> size_t Array<T,N>::totalSize() const
{
    return totSize;
}
template <typename T, size_t N> size_t Array<T,N>::size() const
{
    return totSize;
}

/* Return the size of each array dimension. */
template <typename T, size_t N>
template<typename SIZE_T>
void Array<T,N>::getSize(SIZE_T *result, SIZE_T const noElements) const
{
    if ( N!=size_t(noElements) ) throwError( "Cannot return the size of class 'Array' along each dimension if you don't supply a posize_ter with the same size as the dimension of class 'Array'." );
    for (size_t i=0; i<N; ++i)
        result[i] = SIZE_T( dim[i] );
}
template <typename T, size_t N>
template<typename SIZE_T>
void Array<T,N>::axisSize(SIZE_T *result, SIZE_T const noElements) const
{
    this->getSize( result, noElements );
}

/* Return the size of a given array dimension. */
template <typename T, size_t N> size_t Array<T,N>::getSize(size_t const i) const
{
    if ( i>=N or i<0 ) throwError( "You are asking for the size of the class 'Array' along a dimension larger than the number of dimensions of class 'Array'." );
    return dim[i];
}
template <typename T, size_t N> size_t Array<T,N>::axisSize(size_t const i) const
{
    return this->getSize( i );
}

/* Return the size of each array dimension. */
template <typename T, size_t N>
template<typename SIZE_T>
Array<SIZE_T,1> Array<T,N>::getSize() const
{
    Array<SIZE_T,1> temp(N);
    for (size_t i=0; i<N; ++i)
        temp.data[i] = SIZE_T( this->dim[i] );
    return temp;
}
template <typename T, size_t N>
template<typename SIZE_T>
Array<SIZE_T,1> Array<T,N>::axisSize() const
{
    return this->getSize<SIZE_T>();
}

/* Free the memory posize_ted by 'data'. */
template <typename T, size_t N> void Array<T,N>::freeMemory()
{
    if ( not free and ownData ) delete[] data;
    else throwWarning( "You are trying to free an 'Array<T,N>' object that has already been freed or it was allocated by another class." );
    free = true;
}

/* Check if two arrays have the same size along each dimension. */
template <typename T, size_t N> bool Array<T,N>::similar(Array<T,N> const &newArray) const
{
    for (size_t i=0; i<N; ++i)
        if ( dim[i]!=newArray.dim[i] )
            return false;
    return true;
}



//! Functions that manipulate the values inside 'data' and also let the user acces 'data'
/* Initialize data to a constant value. */
template <typename T, size_t N> void Array<T,N>::assign(T const value)
{
    for (size_t i=0; i<totSize; ++i)
        data[i] = value;
}

/* Initialize data to an array of values. */
template <typename T, size_t N> void Array<T,N>::assign(T const *values, size_t const noElements)
{
    if ( totSize!=noElements ) throwError( "Cannot assign the values to the data in class 'Array' since the array containing those values has a different size from the size of class 'Array'." );
    for (size_t i=0; i<totSize; ++i)
        data[i] = values[i];
}

/* Initialize data to an array of values. */
template <typename T, size_t N> void Array<T,N>::assign(Array<T,N> const &toCopy)
{
    if ( totSize!=toCopy.totalSize() ) throwError( "Cannot assign the values to the data in class 'Array' since the array containing those values has a different size from the size of class 'Array'." );
    for (size_t i=0; i<totSize; ++i)
        data[i] = toCopy.dataValue(i);
}

/* Returns posize_ter to 1st element of 'data'. */
template <typename T, size_t N> T* Array<T,N>::ptrData() const
{
    return data;
}

/* Returns reference to the 'i'-th element of 'data'. */
template <typename T, size_t N> inline T Array<T,N>::at(size_t const i) const
{
    BOOST_ASSERT( i>=0 or i<totSize );// throwError( "The argument of function 'Array::at' is outside the range of values contained in class 'Array'." );
    return data[i];
}

/* Returns reference to the (i,j) element of a 2D array. */
template <typename T, size_t N> inline T Array<T,N>::at(size_t const i, size_t const j ) const
{
    BOOST_ASSERT( N==2 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    return data[ i*indexOffset[0] + j ];
}
template <typename T, size_t N> inline T Array<T,N>::at(size_t const i, size_t const j, size_t const k ) const
{
    BOOST_ASSERT( N==3 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    BOOST_ASSERT( k>=0 or k<dim[2] );
    return data[ i*indexOffset[0] + j*indexOffset[1] + k ];
}
template <typename T, size_t N> inline T Array<T,N>::at(size_t const i, size_t const j, size_t const k, size_t const l ) const
{
    BOOST_ASSERT( N==4 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    BOOST_ASSERT( k>=0 or k<dim[2] );
    BOOST_ASSERT( l>=0 or l<dim[3] );
    return data[ i*indexOffset[0] + j*indexOffset[1] + k*indexOffset[1] + l ];
}





//! Define operators acting on the arrays (on the arrays themselves or on all the entries of the arrays)
/* copy two classes. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator =(T const right)
{
    if ( free )
        throwError( "Error while assigning a constant value in the 'Array' class. You must first allocate memory for the array before assigning values." );
    this->assign( right );
    return *this;
}
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator =(Array<T,N> const &right)
{
    if ( this==&right )
        return *this;
    
    delete[] data;
    for (size_t i=0; i<N; ++i)
        dim[i] = right.axisSize( i );
    this->computeTotalSize();
    this->computeIndexOffset();
    
    data = new (std::nothrow) T[totSize];
    memoryAllocationCheck( data, "Error while allocating memory for the 'Array' class." );
    free = false;
    
    this->assign( right );
    return *this;
}
/* Check if two arrays have exactly the same entries. */
template <typename T, size_t N> bool Array<T,N>::operator ==(Array<T,N> const &right) const
{
    if ( not this->similar(right) ) throwError( "You are comparing two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        if ( data[i]!=right.data[i] )
            return false;
    return true;
}

/* Check if two arrays have different entries. */
template <typename T, size_t N> bool Array<T,N>::operator !=(Array<T,N> const &right) const
{
    if ( not this->similar(right) ) throwError( "You are comparing two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        if ( data[i]!=right.data[i] )
            return true;
    return false;
}

/* Add the entries of the two arrays. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator +(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to add the entries of two class 'Array' objects of different sizes." );
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] + right.data[i];
    return temp;
}

/* Add the entries of another array to the present array. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator +=(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to add the entries of two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        data[i] += right.data[i];
    return (*this);
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator +(T const right)
{
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] + right;
    return temp;
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator +=(T const right)
{
    for (size_t i=0; i<totSize; ++i)
        data[i] += right;
    return (*this);
}

/* Substract the entries of the two arrays. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator -(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to substract the entries of two class 'Array' objects of different sizes." );
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] - right.data[i];
    return temp;
}

/* Substract the entries of another array from the present array. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator -=(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to substract the entries of two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        data[i] -= right.data[i];
    return (*this);
}

/* Substracts from each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator -(T const right)
{
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] - right;
    return temp;
}

/* Substracts from each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator -=(T const right)
{
    for (size_t i=0; i<totSize; ++i)
        data[i] -= right;
    return (*this);
}

/* Multiplies the entries of the two arrays. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator *(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to multiply the entries of two class 'Array' objects of different sizes." );
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] * right.data[i];
    return temp;
}

/* Multiplies the entries of another array to the present array. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator *=(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to multiply the entries of two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        data[i] *= right.data[i];
    return (*this);
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator *(T const right)
{
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] * right;
    return temp;
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator *=(T const right)
{
    for (size_t i=0; i<totSize; ++i)
        data[i] *= right;
    return (*this);
}

/* Add the entries of the two arrays. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator /(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to divide the entries of two class 'Array' objects of different sizes." );
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] / right.data[i];
    return temp;
}

/* Add the entries of another array to the present array. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator /=(Array<T,N> const &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to divide the entries of two class 'Array' objects of different sizes." );
    for (size_t i=0; i<totSize; ++i)
        data[i] /= right.data[i];
    return (*this);
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N> Array<T,N>::operator /(T const right)
{
    Array<T,N> temp(*this);
    for (size_t i=0; i<totSize; ++i)
        temp.data[i] = this->data[i] / right;
    return temp;
}

/* Add to each entry an element of type 'T'. */
template <typename T, size_t N> Array<T,N>& Array<T,N>::operator /=(T const right)
{
    for (size_t i=0; i<totSize; ++i)
        data[i] /= right;
    return (*this);
}





//! Address access operators
/* Returns reference to the 'i'-th element of 'data'. */
template <typename T, size_t N> inline T& Array<T,N>::operator [](size_t const i)
{
    return data[i];
}

/* Returns reference to the 'i'-th element of 'data'. */
template <typename T, size_t N> inline T& Array<T,N>::operator ()(size_t const i)
{
    BOOST_ASSERT( i>=0 or i<totSize );// throwError( "The argument of function 'Array::at' is outside the range of values contained in class 'Array'." );
    return data[i];
}

/* Returns reference to the (i,j) element of a 2D array. */
template <typename T, size_t N> inline T& Array<T,N>::operator ()(size_t const i, size_t const j )
{
    BOOST_ASSERT( N==2 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    return data[ i*indexOffset[0] + j ];
}

/* Returns reference to the (i,j,k) element of a 3D array. */
template <typename T, size_t N> inline T& Array<T,N>::operator ()(size_t const i, size_t const j, size_t const k )
{
    BOOST_ASSERT( N==3 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    BOOST_ASSERT( k>=0 or k<dim[2] );
    return data[ i*indexOffset[0] + j*indexOffset[1] + k ];
}

/* Returns reference to the (i,j,k,l) element of a 4D array. */
template <typename T, size_t N> inline T& Array<T,N>::operator ()(size_t const i, size_t const j, size_t const k, size_t const l )
{
    BOOST_ASSERT( N==4 );
    BOOST_ASSERT( i>=0 or i<dim[0] );
    BOOST_ASSERT( j>=0 or j<dim[1] );
    BOOST_ASSERT( k>=0 or k<dim[2] );
    BOOST_ASSERT( l>=0 or l<dim[3] );
    return data[ i*indexOffset[0] + j*indexOffset[1] + k*indexOffset[1] + l ];
}




//! Additional functions
/* Returns reference to the (i,j,k) element of a 3D array. */
template <typename T, size_t N> void Array<T,N>::print()
{
    switch (N)
    {
        case 1:
            cout << "[";
            for (size_t i1=0; i1<dim[0]; ++i1)
                cout << " " << (*this)(i1);
            cout << " ]\n";
            break;
        case 2:
            cout << "[";
            for (size_t i1=0; i1<dim[0]; ++i1)
            {
                cout << " [";
                for (size_t i2=0; i2<dim[1]; ++i2)
                    cout << " " << (*this)(i1,i2);
                cout << " ]\n";
            }
            cout << " ]\n";
            break;
        case 3:
            cout << "[";
            for (size_t i1=0; i1<dim[0]; ++i1)
            {
                cout << " [";
                for (size_t i2=0; i2<dim[1]; ++i2)
                {
                    cout << " [";
                    for (size_t i3=0; i3<dim[2]; ++i3)
                        cout << " " << (*this)(i1,i2,i3);
                    cout << " ]\n";
                }
                cout << " ]\n";
            }
            cout << " ]\n";
            break;
        case 4:
            cout << "[";
            for (size_t i1=0; i1<dim[0]; ++i1)
            {
                cout << " [";
                for (size_t i2=0; i2<dim[1]; ++i2)
                {
                    cout << " [";
                    for (size_t i3=0; i3<dim[2]; ++i3)
                    {
                        cout << " [";
                        for (size_t i4=0; i4<dim[3]; ++i4)
                            cout << " " << (*this)(i1,i2,i3,i4);
                        cout << " ]\n";
                    }
                    cout << " ]\n";
                }
                cout << " ]\n";
            }
            cout << " ]\n";
            break;
        default:
            throwError( "The function 'Array::prsize_t()' was not defined for dimensions outside the size_terval 1 - 4." );
    }
}



//! Private functions
/* Sets the dimension sizes size_ternally. */
template <typename T, size_t N>
template <typename SIZE_T>
void Array<T,N>::setDimensionSize(SIZE_T const *dimSize, SIZE_T const noElements)
{
    if ( N!=noElements ) throwError( "When creating the class 'Array' since the dimension of the array is different from the number of elements specified to the class constructor." );
    for (SIZE_T i=0; i<noElements; ++i )
        if ( dimSize[i]<=0 ) throwError( "When creating the class 'Array' since the dimension ", i, " has <=0 value. All array dimensions must be positive!" );
    
    for (size_t i=0; i<N; ++i)
        dim[i] = dimSize[i];
    
    this->computeTotalSize();
    this->computeIndexOffset();
}

/* computes total number of elements in the array */
template <typename T, size_t N> void Array<T,N>::computeTotalSize()
{
    totSize = 1;
    for (size_t i=0; i<N; ++i)
        totSize *= dim[i];
}

/* computes total number of elements in the array */
template <typename T, size_t N> void Array<T,N>::computeIndexOffset()
{
    if ( N==1 )
        indexOffset[0] = 1;
    else
    {
        indexOffset[N-1] = 1;
        for (int i=int(N-2); i>=0; --i)
            indexOffset[i] = indexOffset[i+1] * dim[i+1];
    }
}

/* Returns the value of of i-th data element. */
template <typename T, size_t N> T Array<T,N>::dataValue(size_t const i) const
{
    return data[i];
}



#endif
