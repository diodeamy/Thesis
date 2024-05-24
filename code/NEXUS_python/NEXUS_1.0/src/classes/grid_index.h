
/*

This class implements a grid index along a dimension for a periodic grid, i.e. the grid index can ONLY take values from 0 to n, where n is the size of the grid along that direction. 

NOTE:
     1) To implement the concept the class saves twice as much memory than it would need to keep track of the index only.
     2) This computation can be applied to other problems, like the position of particles in a periodic box (even though ++ and -- operators do not make sense).

*/

#ifndef GRID_INDEX_HEADER
#define GRID_INDEX_HEADER


#include <misc.h>








template <typename T>
struct Grid_index
{
    Grid_index(T const grid_size);                // constructors
    Grid_index(Grid_index const & other);
    Grid_index(T const grid_size, T const initial_value);
    
    T get() const;                            // return the value of 'index'
    T operator ()() const;                    // return the value of 'index'
    T operator ()(T const value);             // set the 'index' to 'value' and return 'index'
    T grid() const;                           // return the value of 'gridSize'
    bool similar(Grid_index &newIndex) const; // check if two Grid_index objects are similar (i.e. have the same 'gridSize')
    
    Grid_index& operator ++();                // increase the value of 'index' by one
    Grid_index operator ++(int);
    Grid_index& operator --();                // decrease the value of 'index' by one
    Grid_index operator --(int);
    
    Grid_index& operator =(T right);          // assign this value to the 'index' variable
    Grid_index& operator =(Grid_index right); // copy the right value to the left one
    
    bool operator ==(Grid_index &right) const;// define equality check (same 'index' and 'gridSize')
    bool operator !=(Grid_index &right) const;// define inequality check (different 'index' or 'gridSize')
    
    Grid_index operator +(Grid_index &right);  // defines addition of two Grid_index objects (add 'index' value if the two are similar)
    Grid_index& operator +=(Grid_index &right);// defines self addition
    Grid_index operator +(T right);            // addition to a 'T' value
    Grid_index& operator +=(T right);
    Grid_index operator -(Grid_index &right);  // defines substraction of two Grid_index objects (substract 'index' value if the two are similar)
    Grid_index& operator -=(Grid_index &right);
    Grid_index operator -(T right);
    Grid_index& operator -=(T right);
    Grid_index operator *(T right);            // defines multiplication to 'T' value
    Grid_index& operator *=(T right); 
    Grid_index operator /(T right);            // defines multiplication to 'T' value
    Grid_index& operator /=(T right);
    Grid_index operator %(T right);            // apply the % operator on the 'index' variable
    Grid_index& operator %=(T right);
    
    
    private:
        T gridSize;
        T index;
        void translate();                      // translate the value of index to be between 0 to gridSize
};







/* Constructors: */
template <typename T> Grid_index<T>::Grid_index(T const grid_size)
{
    gridSize = grid_size;
    index = T(0);
}
template <typename T> Grid_index<T>::Grid_index(Grid_index<T> const &other)
{
    gridSize = other.grid();
    index = other.get();
}
template <typename T> Grid_index<T>::Grid_index(T const grid_size, T const initial_value)
{
    gridSize = grid_size;
    index = initial_value;
    this->translate(); // translate 'index' to be between 0 to 'gridSize'
}

/* Output the value of 'index'. */
template <typename T> T Grid_index<T>::get() const
{
    return index;
}
template <typename T> T Grid_index<T>::operator ()() const
{
    return index;
}
template <typename T> T Grid_index<T>::operator ()(T const value)
{
    index = value;
    this->translate(); // translate 'index' to be between 0 to 'gridSize'
    return index;
}
/* Output the value of 'gridSize'. */
template <typename T> T Grid_index<T>::grid() const
{
    return gridSize;
}
/* Check if both objects have the same 'gridSize'. */
template <typename T> bool Grid_index<T>::similar(Grid_index <T>&other) const
{
    if ( gridSize==other.grid() )
        return true;
    return false;
}



//! Define the increment and decrement operators
template <typename T> Grid_index<T>& Grid_index<T>::operator ++()
{
    ++index;
    this->translate();
    return *this;
}
template <> Grid_index<int>& Grid_index<int>::operator ++()
{
    ++index;
    if ( index==gridSize )
        index = 0;
    return *this;
}
template <typename T> Grid_index<T> Grid_index<T>::operator ++(int)
{
    Grid_index<T> temp( *this );
    ++(*this);
    return temp;
}
template <typename T> Grid_index<T>& Grid_index<T>::operator --()
{
    --index;
    this->translate();
    return *this;
}
template <> Grid_index<int>& Grid_index<int>::operator --()
{
    --index;
    if ( index==-1 )
        index = gridSize-1;
    return *this;
}
template <typename T> Grid_index<T> Grid_index<T>::operator --(int)
{
    Grid_index<T> temp( *this );
    --(*this);
    return temp;
}




//! Redefine the = operator
/* Initialize 'index'. */
template <typename T> Grid_index<T>& Grid_index<T>::operator =(T right)
{
    index = right;
    this->translate(); // translate 'index' to be between 0 to 'gridSize'
    return (*this);
}

/* Transfer 'index' value between objects. */
template <typename T> Grid_index<T>& Grid_index<T>::operator =(Grid_index<T> right)
{
    this->gridSize = right.grid();
    this->index = right.get();
    return (*this);
}



//! Boolean equality/inequality operators
/* Check for equality of two objects. */
template <typename T> bool Grid_index<T>::operator ==(Grid_index<T> &right) const
{
    if ( this->similar(right) and index==right.index )
       return true;
    return false;
}

/* Check for inequality of two objects. */
template <typename T> bool Grid_index<T>::operator !=(Grid_index<T> &right) const
{
    if ( this->similar(right) and index==right.index )
       return false;
    return true;
}



//! Mathematical operators
/* Define addition of two objects. */
template <typename T> Grid_index<T> Grid_index<T>::operator +(Grid_index<T> &right)
{
    Grid_index<T> temp(*this);
    temp += right;
    return temp;
}

/* Define self addition of two objects. */
template <typename T> Grid_index<T>& Grid_index<T>::operator +=(Grid_index<T> &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to add two 'Grid_index' objects of different grid sizes." );
    this->index += right.get();
    this->translate();
    return (*this);
}

/* Define addition to 'T' value. */
template <typename T> Grid_index<T> Grid_index<T>::operator +(T right)
{
    Grid_index<T> temp(*this);
    temp += right;
    return temp;
}

/* Define self addition to 'T' value. */
template <typename T> Grid_index<T>& Grid_index<T>::operator +=(T right)
{
    index += right;
    this->translate();
    return (*this);
}

/* Define substract of two objects. */
template <typename T> Grid_index<T> Grid_index<T>::operator -(Grid_index<T> &right)
{
    Grid_index<T> temp(*this);
    temp -= right;
    return temp;
}

/* Define self substract of two objects. */
template <typename T> Grid_index<T>& Grid_index<T>::operator -=(Grid_index<T> &right)
{
    if ( not this->similar(right) ) throwError( "You are trying to substract two 'Grid_index' objects of different grid sizes." );
    this->index -= right.get();
    this->translate();
    return (*this);
}

/* Define substract to 'T' value. */
template <typename T> Grid_index<T> Grid_index<T>::operator -(T right)
{
    Grid_index<T> temp(*this);
    temp -= right;
    return temp;
}

/* Define self substract to 'T' value. */
template <typename T> Grid_index<T>& Grid_index<T>::operator -=(T right)
{
    index -= right;
    this->translate();
    return (*this);
}

/* Define multiplication to 'T' value. */
template <typename T> Grid_index<T> Grid_index<T>::operator *(T right)
{
    Grid_index<T> temp(*this);
    temp *= right;
    return temp;
}
template <typename T> Grid_index<T>& Grid_index<T>::operator *=(T right)
{
    index *= right;
    this->translate();
    return (*this);
}

/* Define division to 'T' value. */
template <typename T> Grid_index<T> Grid_index<T>::operator /(T right)
{
    Grid_index<T> temp(*this);
    temp /= right;
    return temp;
}
template <typename T> Grid_index<T>& Grid_index<T>::operator /=(T right)
{
    index /= right;
    this->translate();
    return (*this);
}

/* Define the modulus operator. */
template <typename T> Grid_index<T> Grid_index<T>::operator %(T right)
{
    Grid_index<T> temp(*this);
    temp %= right;
    return temp;
}
template <typename T> Grid_index<T>& Grid_index<T>::operator %=(T right)
{
    index = index % right;
    this->translate();
    return *this;
}



//! Private functions
/* Get the value of index in the range 0 to gridsize. */
template <typename T> void Grid_index<T>::translate()
{
    while ( index<T(0) )          // translate 'index' to be between 0 to 'gridSize'
        index += gridSize;
    while ( index>=gridSize )
        index -= gridSize;
}
template <> void Grid_index<int>::translate()
{
    if ( index<0 or index>=gridSize )    // translate 'index' to be between 0 to 'gridSize'
        index = index % gridSize;
    if ( index<0 )
        index += gridSize;
}





#endif
