/*

This defines a physical vector class with the properties and operations specific to N-dimensional vectors.

NOTE:
     1) The functions a(i) and a.at(i) do a range check for out of range elements, a[i] does not do so.
     2) The range check can be disable using the compiler directive BOOST_DISABLE_ASSERTS.
*/


#ifndef VECTOR_HEADER
#define VECTOR_HEADER


#include <cmath>
#include <boost/assert.hpp>

#include <defines.h>


/* A class that has implemented the physical properties of a N-dimensional vector. */
template<typename T, int N>
struct Vector
{
    T comp[N];
    
    Vector();                 // initialize all components to 0
    Vector(T const value);
    Vector(Vector const &other);
    
    void assign(T const value);         // assign to all entries a cosntant value
    void assign(T *ptr, int const noElements);// assign from an array of values
    void print() const;                 // output the elements of the vector
    T* ptrData();                       // return a pointer to the vector components
    
    Vector& operator =(Vector const &other); // define vector assigment
    
    Vector operator +(Vector const &right);  // adds two vectors and saves the result in the 3rd
    Vector& operator +=(Vector const &right);// adds another vector to the present one
    Vector operator -(Vector const &right);
    Vector& operator -=(Vector const &right);
    
    T& operator [](int const i);        // adds a pointer to the 'i'-th component
    T& operator ()(int const i);        // adds a pointer to the 'i'-th component - does range check
    T& at(int const i);                 // adds a pointer to the 'i'-th component - does range check
    
    T operator *(Vector const &right) const;// computes scalar product
    Vector operator *(T const right);       // multiplication by scalar
    Vector& operator *=(T const right);
    Vector operator /(T const right);       // division by scalar
    Vector& operator /=(T const right);
    
    T cosinusAngle(Vector const &y) const;  // computes the cos(angle between two vectors)
    T magnitude() const;                    // get magnitude of vector
    
    void normalize();                       // normalize the components to get magnitude of vector 1
    void normalize(Vector<T,N> vec[], int const size);   // normalize the magnitude of an array of vectors
    
    private:
        T value(int const i) const;     // returns the entries of the Vector
};





// Constructors
template<typename T, int N> Vector<T,N>::Vector()
{
    for (int i=0; i<N; ++i)
        comp[i] = T(0);
}
template<typename T, int N> Vector<T,N>::Vector(T const value)
{
   this->assign( value );
}
template<typename T, int N> Vector<T,N>::Vector(Vector<T,N> const &other)
{
    for (int i=0; i<N; ++i)
        comp[i] = other.value(i);
}



//! Assign values to the vector entries
template<typename T, int N> void Vector<T,N>::assign(T const value)
{
    for (int i=0; i<N; ++i)
        comp[i] = value;
}
template<typename T, int N> void Vector<T,N>::assign(T *ptr, int const noElements)
{
    BOOST_ASSERT( N==noElements );
    for (int i=0; i<N; ++i)
        comp[i] = ptr[i];
}
template<typename T, int N> T* Vector<T,N>::ptrData()
{
    return &(comp[0]);
}



//! Define the equality operator
template<typename T, int N> Vector<T,N>& Vector<T,N>::operator =(Vector<T,N> const &other)
{
    for (int i=0; i<N; ++i)
        comp[i] = other.value(i);
    return (*this);
}



//! mathematical operations
/* Define addition of 2 vector.*/
template<typename T, int N> inline Vector<T,N> Vector<T,N>::operator +(Vector<T,N> const &right)
{
    Vector<T,N> temp;
    for (int i=0; i<N; ++i)
        temp[i] = comp[i] + right.value(i);
    return temp;
}

/* Define addition to oneself.*/
template<typename T, int N> inline Vector<T,N>& Vector<T,N>::operator +=(Vector<T,N> const &right)
{
    for (int i=0; i<N; ++i)
        comp[i] += right.value(i);
    return *this;
}

/* Define difference of 2 vector.*/
template<typename T, int N> inline Vector<T,N> Vector<T,N>::operator -(Vector<T,N> const &right)
{
    Vector<T,N> temp;
    for (int i=0; i<N; ++i)
        temp[i] = comp[i] - right.value(i);
    return temp;
}

/* Define difference from oneself.*/
template<typename T, int N> inline Vector<T,N>& Vector<T,N>::operator -=(Vector<T,N> const &right)
{
    for (int i=0; i<N; ++i)
        this->comp[i] -= right.value(i);
    return *this;
}


//! Define element access operators
/* Define [] operator = element access.*/
template<typename T, int N> inline T& Vector<T,N>::operator [](int const i)
{
    return this->comp[i];
}
template<typename T, int N> inline T& Vector<T,N>::operator ()(int const i)
{
    BOOST_ASSERT( i>=0 and i<N );
    return this->comp[i];
}
template<typename T, int N> inline T& Vector<T,N>::at(int const i)
{
   return (*this)(i);
}



//! Define vector operations and operations with scalars
/* Define the scalar product.*/
template<typename T, int N> inline T Vector<T,N>::operator *(Vector<T,N> const &right) const
{
    T temp = T(0.);
    for (int i=0; i<N; ++i)
        temp += comp[i]*right.value(i);
    return temp;
}

/* Define multiplication with a scalar.*/
template<typename T, int N> inline Vector<T,N> Vector<T,N>::operator *(T const right)
{
    Vector<T,N> temp;
    for (int i=0; i<N; ++i)
        temp[i] = comp[i]*right;
    return temp;
}
template<typename T, int N> inline Vector<T,N>& Vector<T,N>::operator *=(T const right)
{
    for (int i=0; i<N; ++i)
        comp[i] *= right;
    return (*this);
}

/* Define division with a scalar.*/
template<typename T, int N> inline Vector<T,N> Vector<T,N>::operator /(T const right)
{
    Vector<T,N> temp;
    for (int i=0; i<N; ++i)
        temp[i] = comp[i]/right;
    return temp;
}
template<typename T, int N> inline Vector<T,N>& Vector<T,N>::operator /=(T const right)
{
    for (int i=0; i<N; ++i)
        comp[i] /= right;
    return (*this);
}



/* Define the cosinus(angle between two vectors).*/
template<typename T, int N> inline T Vector<T,N>::cosinusAngle(Vector<T,N> const &other) const
{
    T temp = this->magnitude() * other.magnitude();
    T temp1 = (*this) * other;
    return (temp1 / temp);
}



/* Compute the magnitude of the vector. */
template<typename T, int N> inline T Vector<T,N>::magnitude() const
{
    T temp = T(0);
    for (int i=0; i<N; ++i)
        temp += comp[i]*comp[i];
    return sqrt( temp );
}


/* Normalize the vector to magnitude 1. */
template<typename T, int N> inline void Vector<T,N>::normalize()
{
    T temp = this->magnitude();
    for (int i=0; i<N; ++i)
        comp[i] /= temp;
}


/* Normalize an array of vectors to magnitude 1. */
template<typename T, int N> void Vector<T,N>::normalize(Vector<T,N> vec[], int const size)
{
    for (int i=0; i<size; ++i)
        vec[i].normalize();
}

/* Output the contents of the vector */
template<typename T, int N> void Vector<T,N>::print() const
{
    cout << "\n";
    for (int i=0; i<N; ++i)
        cout << comp[i] <<"   ";
    cout << "\n";
}



//! Private functions
/* Return the 'i'-th component of the vector. */
template<typename T, int N> inline T Vector<T,N>::value(int const i) const
{
    return comp[i];
}




















// computes the magnitude of a vector
template<typename T> inline T magnitude(T vec[])
{
    T temp = T(0.);
    for (int j=0; j<NO_DIM; ++j)
        temp += vec[j]*vec[j];
    return std::sqrt( temp );
}
template<typename T> inline T scalarProduct(T vec1[], T vec2[])
{
    T temp = T(0.);
    for (int j=0; j<NO_DIM; ++j)
        temp += vec1[j]*vec2[j];
    return temp;
}
// normalizes a vector. Returns false if the operation failed.
template<typename T> inline bool normalize(T vec[])
{
    T mag = magnitude(vec);
    for (int j=0; j<NO_DIM; ++j)
        vec[j] /= mag;
    if ( mag<1.e-5 ) return false;
    return true;
}

// computes the cosinus between two vectors
template<typename T> inline T cosAngle(T vec1[], T vec2[])
{
    return scalarProduct(vec1,vec2) /( magnitude(vec1)*magnitude(vec2) );
}

// distance between two points
template<typename T> inline T distance(T pos1[], T pos2[])
{
    T res[NO_DIM], dis = 0.;
    for (int j=0; j<NO_DIM; ++j)
    {
        res[j] += pos1[j] - pos2[j];
        dis += res[j]*res[j];
    }
    return std::sqrt(dis);
}
template<typename T> inline T distance(T pos1[], T pos2[], T result[])
{
    T dis = 0.;
    for (int j=0; j<NO_DIM; ++j)
    {
        result[j] += pos1[j] - pos2[j];
        dis += result[j]*result[j];
    }
    return std::sqrt(dis);
}


// this functions gives the direction along which to displace an object from pos1 to pos2 along a given direction 'direction'.
template<typename T>
inline T distanceAlongDirection(T pos1[],
                                T pos2[],
                                T direction[],
                                T moveDirection[] )
{
    T mag = magnitude(direction);    // the magnitude of the vector
    T temp[NO_DIM];  //the displacement vector
    for (int i=0; i<NO_DIM; ++i)
        temp[i] = pos2[i] - pos1[i];
    T factor = scalarProduct(temp,direction)/(mag*mag);
    T distance = 0.;
    for (int i=0; i<NO_DIM; ++i)
    {
        moveDirection[i] = factor * direction[i];      // the distance along the given direction
        distance += moveDirection[i]*moveDirection[i];
    }
    return std::sqrt(distance);
}
// this functions gives the direction along which to displace an object from pos1 to pos2 along a direction perpendicular to 'perpDirection'.
template<typename T>
inline T distanceAlongPerpendicularDirection(T pos1[],
                                             T pos2[],
                                             T perpDirection[],
                                             T moveDirection[] )
{
    T mag = magnitude(perpDirection);    // the magnitude of the vector
    T temp[NO_DIM];  //the displacement vector
    for (int i=0; i<NO_DIM; ++i)
        temp[i] = pos2[i] - pos1[i];
    T factor = scalarProduct(temp,perpDirection)/(mag*mag);
    T distance = 0.;
    for (int i=0; i<NO_DIM; ++i)
    {
        moveDirection[i] = temp[i] - factor*perpDirection[i];  // substract the perpendicular direction to get the distance perpendicular to 'perpDirection'
        distance += moveDirection[i]*moveDirection[i];
    }
    return std::sqrt(distance);
}



// checks if a particle is in the box
template<typename T>
inline bool particleInBox(T *pos,
                          T *box)
{
    for (int i=0; i<NO_DIM; ++i)
        if (pos[i]<box[2*i] or pos[i]>box[2*i+1])
            return false;
    return true;
}



/* This function computes the N-th eigenvector given N-1 eigenvectors. Works for 2D and 3D cases. */
template<typename T>
void lastEigenvector(T *e,
                     T *result)
{
#if NO_DIM==2
    result[0] = -e[1];
    result[1] = e[0];
#elif NO_DIM==3
    result[0] = e[2]*e[4] - e[1]*e[5];
    result[1] = e[0]*e[5] - e[2]*e[3];
    result[2] = e[1]*e[3] - e[0]*e[4];
    
    normalize( result );   // normalize the vector to magnitude 1
#endif
}







#endif
