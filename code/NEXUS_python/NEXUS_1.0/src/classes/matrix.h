


#ifndef MATRIX_HEADER
#define MATRIX_HEADER

#include <iostream>
#include <cmath>
#include <boost/assert.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <miscellaneous.h>
#include <defines.h>

#define SWAP_VARIABLES(lvX,lvY,VARIABLE_TYPE) {VARIABLE_TYPE temporary=lvX; lvX=lvY; lvY=temporary;}


namespace MATRIX
{

template <typename T, size_t N>
struct Square
{
    static const size_t size; // the number of elements in the matrix
    static const T tolerance; // constant that characterizes the numerical inaccuracy
    static const T largeTolerance; // constant that characterizes the numerical accuracy when switching from float to double
    T elem[N*(N+1)/2];  // don't need to store all the elements
    
    // functions
    Square();                        // constructor
    Square(T *elements);             // constructor
    void data(T *elements);          // change the entries of the matrix
    
    // overloading of operators
    typedef Square<T,N> SQ;
    
    bool operator ==(SQ &right);  // return true if the two matrices are the same ( checks each entry individually )
    bool operator !=(SQ &right);  // return true if the two matrices are not the same ( checks each entry individually )
    SQ  operator +(SQ &right);    // adds the elements of the two matrices ( each entry individually )
    SQ& operator +=(SQ &right);   // adds the elements of another matrix to the present one
    SQ  operator -(SQ &right);    // substracts the elements of the two matrices ( each entry individually )
    SQ& operator -=(SQ &right);
    SQ  operator *(SQ  &right);   // computes the product of the two matrices
    SQ operator *=(SQ &right);
    SQ  operator *(T const right);      // multiplies all entries by the same element of type 'T'
    SQ& operator *=(T const right);
    SQ  operator /(T const right);      // divides all entries by the same element of type 'T'
    SQ& operator /=(T const right); 
    
    inline T& operator [](size_t const i);              // returns a reference to the 'i'-th entry of 'data'
    inline T& operator ()(size_t const i, size_t const j); // returns a reference to the (i,j) entry of the matrix
    
    // matrix operations
    inline T det();                            // returns the determinant
    inline T trace();                          // returns the trace
    int rank();                                // returns the rank of the matrix
    void solveLinearSystem(T *right, T *solution);  // returns the rank of the matrix
    template <typename COMPLEX> void eigenvalues(COMPLEX *eigenvalues, COMPLEX *eigenvectors = NULL, bool allEigenvectors=true);   // computes the matrix eigenvalues and eigenvectors (if an eigenvectors pointer is supplied)
    inline T scale();       // returns the sum of the absolute values of the matrix entries (used to normalize the entries to reduce errors due to too large/small numbers)
    
    private:
        int rank2();
        int rank3();
        int rankN();
        void solveLinearSystem2(T *right, T *solution);
        void solveLinearSystem3(T *right, T *solution);
        void solveLinearSystemN(T *right, T *solution);
        template <typename COMPLEX> void eigenvaluesN(COMPLEX *eigenvalues, COMPLEX *eigenvectors, bool allEigenvectors);
};
template <typename T, size_t N> const size_t Square<T,N>::size = N*N;
template <typename T, size_t N> const T Square<T,N>::tolerance = sizeof(T)==4 ? 1.e-5 : 1.e-14 ;
template <typename T, size_t N> const T Square<T,N>::largeTolerance = sizeof(T)==4 ? 1.e-4 : 1.e-13 ;




template <typename T, size_t N>
        struct Symmetric
{
    static const size_t size; // the number of elements in the matrix
    static const T tolerance; // constant that characterizes the numerical accuracy
    static const T largeTolerance; // constant that characterizes the numerical accuracy when switching from float to double
    T elem[N*(N+1)/2];  // don't need to store all the elements
    
    
    // functions
    Symmetric();                        // constructor
    Symmetric(T *elements);             // constructor
    void data(T *elements);             // change the entries of the matrix
    
    // overloading of operators
    typedef Symmetric<T,N> SM;
    
    bool operator ==(SM &right);  // return true if the two matrices are the same ( checks each entry individually )
    bool operator !=(SM &right);  // return true if the two matrices are not the same ( checks each entry individually )
    SM  operator +(SM &right);    // adds the elements of the two matrices ( each entry individually )
    SM& operator +=(SM &right);   // adds the elements of another matrix to the present one
    SM  operator -(SM &right);    // substracts the elements of the two matrices ( each entry individually )
    SM& operator -=(SM &right);
    SM  operator *(SM  &right);    // computes the product of the two matrices
    SM operator *=(SM &right);
    SM  operator *(T const right);      // multiplies all entries by the same element of type 'T'
    SM& operator *=(T const right);
    SM  operator /(T const right);      // divides all entries by the same element of type 'T'
    SM& operator /=(T const right); 
    
    inline T& operator [](size_t const i);              // returns a reference to the 'i'-th entry of 'data'
    inline T& operator ()(size_t const i, size_t const j); // returns a reference to the (i,j) entry of the matrix
    
    
    // matrix operations
    inline T det();                            // returns the determinant
    inline T trace();                          // returns the trace
    int rank();                                // returns the rank of the matrix
    void solveLinearSystem(T *right, T *solution);  // returns the rank of the matrix
    void eigenvalues(T *eigenvalues, T *eigenvectors = NULL, bool allEigenvectors=true);   // computes the matrix eigenvalues and eigenvectors (if an eigenvectors pointer is supplied)
    inline T scale();       // returns the sum of the absolute values of the matrix entries (used to normalize the entries to reduce errors due to too large/small numbers)
    
    private:
        void eigenvalues2(T *eigenvalues, T *eigenvectors, bool allEigenvectors);
        void eigenvalues3(T *eigenvalues, T *eigenvectors, bool allEigenvectors);
        void eigenvaluesN(T *eigenvalues, T *eigenvectors, bool allEigenvectors);
};
template <typename T, size_t N> const size_t Symmetric<T,N>::size = N*(N+1)/2;
template <typename T, size_t N> const T Symmetric<T,N>::tolerance = sizeof(T)==4 ? 1.e-5 : 1.e-14 ;
template <typename T, size_t N> const T Symmetric<T,N>::largeTolerance = sizeof(T)==4 ? 1.e-4 : 1.e-13 ;






//! Functions for the 'template <typename T, size_t N> Square<T,N>' class
// constructors
template <typename T, size_t N> Square<T,N>::Square()
{
    for (size_t i=0; i<size; ++i)
        elem[i] = T(0.);
}
template <typename T, size_t N> Square<T,N>::Square(T *elements)
{
    for (size_t i=0; i<size; ++i)
        elem[i] = elements[i];
}
template <typename T, size_t N> void Square<T,N>::data(T *elements)
{
    for (size_t i=0; i<size; ++i)
        elem[i] = elements[i];
}


// overloading of operators
template <typename T, size_t N> bool Square<T,N>::operator ==(Square<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        if ( elem[i]!=right[i] )
            return false;
    return true;
}
template <typename T, size_t N> bool Square<T,N>::operator !=(Square<T,N> &right)
{
    return not (*this)==right;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator +(Square<T,N> &right)
{
    Square<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] + right[i];
    return temp;
}
template <typename T, size_t N> Square<T,N>& Square<T,N>::operator +=(Square<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] += right[i];
    return *this;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator -(Square<T,N> &right)
{
    Square<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] - right[i];
    return temp;
}
template <typename T, size_t N> Square<T,N>& Square<T,N>::operator -=(Square<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] -= right[i];
    return *this;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator *(Square<T,N> &right)
{
    Square<T,N> temp;
    for (size_t i=0; i<N; ++i)
        for (size_t j=i; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                temp(i,j) += (*this)(i,k) * right(k,j);
    return temp;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator *=(Square<T,N> &right)
{
    return (*this)*right;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator *(T const right)
{
    Square<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] * right;
    return temp;
}
template <typename T, size_t N> Square<T,N>& Square<T,N>::operator *=(T const right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] *= right;
    return *this;
}
template <typename T, size_t N> Square<T,N> Square<T,N>::operator /(T const right)
{
    Square<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] / right;
    return temp;
}
template <typename T, size_t N> Square<T,N>& Square<T,N>::operator /=(T const right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] /= right;
    return *this;
}

template <typename T, size_t N> T& Square<T,N>::operator [](size_t const i)
{
    BOOST_ASSERT( i<size );
    return elem[i];
}
template <typename T, size_t N> T& Square<T,N>::operator ()(size_t const i, size_t const j)
{
    BOOST_ASSERT( i<N and j<N );
    return elem[i*N+j];
}


// Matrix operations
/* Computes the determinant of the matrix. */
template <typename T, size_t N> T Square<T,N>::det()
{
    if (N==2)
        return elem[0]*elem[3] - elem[1]*elem[2];
    else if (N==3)
        return elem[0]*elem[4]*elem[8] + elem[1]*elem[5]*elem[6] + elem[2]*elem[3]*elem[7] - elem[2]*elem[4]*elem[6] - elem[0]*elem[5]*elem[7] - elem[1]*elem[3]*elem[8];
    else
        throwError( "The 'MATRIX::Square::det()' function is implemented only for 2x2 and 3x3 matrices." );
}

/* Computes the trace of the matrix. */
template <typename T, size_t N> T Square<T,N>::trace()
{
    if (N==2)
        return elem[0]+elem[3];
    else if (N==3)
        return elem[0]+elem[4]+elem[8];
    else
    {
        T temp = T(0.);
        for (size_t i=0; i<N; ++i)
            temp += (*this)(i,i);
        return temp;
    }
}

/* Computes the sum of the absolute values of the matrix entries. */
template <typename T, size_t N> T Square<T,N>::scale()
{
    T temp = T(0.);
    for (size_t i=0; i<size; ++i)
        temp += std::fabs( elem[i] );
    return temp;
}


/* Returns matrix properties as rank or solves a system of linear equations. */
template <typename T, size_t N> int Square<T,N>::rank()
{
    // scale the matrix elements
    T const scale = this->scale();  // returns the sum of the absolute values of the elements
    Square<T,N> tempMat = (*this) / scale;
    
    if (N==2)
        return tempMat.rank2();
    else if (N==3)
        return tempMat.rank3();
    else
        return tempMat.rankN();
}
template <typename T, size_t N> int Square<T,N>::rank2()
{
    int rank = 0;
    if ( std::fabs(elem[2])>std::fabs(elem[0]) )
    {
        SWAP_VARIABLES( elem[0], elem[2], T );
        SWAP_VARIABLES( elem[1], elem[3], T );
    }
    if ( std::fabs(elem[0])>tolerance )
    {
        ++rank;
        elem[3] -= elem[1] * elem[2]/elem[0];
    }
    
    if ( std::fabs(elem[0])>tolerance ) ++rank;
    return rank;
}
template <typename T, size_t N> int Square<T,N>::rank3()
{
    int rank = 0;
    // reduce the 1st column
    if ( std::fabs(elem[3])>std::fabs(elem[0]) )
        for (int i=0; i<N; ++i) SWAP_VARIABLES( elem[0+i], elem[3+i], T );
    if ( std::fabs(elem[6])>std::fabs(elem[0]) )
        for (int i=0; i<N; ++i) SWAP_VARIABLES( elem[0+i], elem[6+i], T );
    
    if ( std::fabs(elem[0])>tolerance )
    {
        ++rank;
        elem[4] -= elem[1] * elem[3]/elem[0];
        elem[5] -= elem[2] * elem[3]/elem[0];
        elem[7] -= elem[1] * elem[6]/elem[0];
        elem[8] -= elem[2] * elem[6]/elem[0];
    }
    
    // reduce the 2nd column
    if ( std::fabs(elem[7])>std::fabs(elem[4]) )
        for (int i=1; i<N; ++i) SWAP_VARIABLES( elem[3+i], elem[6+i], T );
    if ( std::fabs(elem[4])>tolerance )
    {
        ++rank;
        elem[8] -= elem[5] * elem[7]/elem[4];
    }
    
    // reduce the 3rd column
    if ( std::fabs(elem[8])>tolerance )
        ++rank;
    
    return rank;
}
template <typename T, size_t N> int Square<T,N>::rankN()
{
    throwError( "There is no method implemented to return the rank of N>3 matrices." );
}


template <typename T, size_t N> void Square<T,N>::solveLinearSystem(T *right, T *solution)
{
    T const scale = this->scale();  // returns the sum of the absolute values of the elements
    Square<T,N> tempMat = (*this) / scale;
    
    if (N==2)
        tempMat.solveLinearSystem2(right,solution);
    else if (N==3)
        tempMat.solveLinearSystem3(right,solution);
    else
        tempMat.solveLinearSystemN(right,solution);
    
    for (int i=0; i<N; ++i)
        solution[i] *= scale;
}
template <typename T, size_t N> void Square<T,N>::solveLinearSystem2(T *right, T *solution)
{
    int order[N];
    T result[N];
    for (int i=0; i<N; ++i)
        order[i] = i;
    
    // reduce the first column
    if ( std::fabs(elem[2])>std::fabs(elem[0]) )
    {
        for (int i=0; i<N; ++i) SWAP_VARIABLES( elem[0+i], elem[N+i], T );
        SWAP_VARIABLES( right[0], right[1], T );
        order[0] = 1; order[1] = 0;
    }
    if ( std::fabs(elem[0])>tolerance )
    {
        elem[3] -= elem[1] * elem[2]/elem[0];
        right[1] -= right[0] * elem[2]/elem[0];
    }
    else throwError( "Asking to solve a linear system for a 2x2 matrix of rank < 2. The solution is not unique." );
    
    // reduce the 2nd column
    if ( not std::fabs(elem[3])>tolerance )
        throwError( "Asking to solve a linear system for a 2x2 matrix of rank < 2. The solution is not unique." );
    
    // get the solution
    result[1] = right[1]/elem[3];
    result[0] = (right[0] - elem[1]*result[1]) /elem[0];
    for (int i=0; i<N; ++i)
        solution[ order[i] ] = result[i];
}
template <typename T, size_t N> void Square<T,N>::solveLinearSystem3(T *right, T *solution)
{
    int order[N];
    T result[N];
    for (int i=0; i<N; ++i)
        order[i] = i;
    
    // reduce the first column
    if ( std::fabs(elem[3])>std::fabs(elem[0]) )
    {
        for (int i=0; i<N; ++i) SWAP_VARIABLES( elem[0+i], elem[3+i], T );
        SWAP_VARIABLES( right[0], right[1], T );
        SWAP_VARIABLES( order[0], order[1], int );
    }
    if ( std::fabs(elem[6])>std::fabs(elem[0]) )
    {
        for (int i=0; i<N; ++i) SWAP_VARIABLES( elem[0+i], elem[6+i], T );
        SWAP_VARIABLES( right[0], right[2], T );
        SWAP_VARIABLES( order[0], order[2], int );
    }
    if ( std::fabs(elem[0])>tolerance )
    {
        elem[4] -= elem[1] * elem[3]/elem[0];
        elem[5] -= elem[2] * elem[3]/elem[0];
        right[1] -= right[0] * elem[3]/elem[0];
        elem[7] -= elem[1] * elem[6]/elem[0];
        elem[8] -= elem[2] * elem[6]/elem[0];
        right[2] -= right[0] * elem[6]/elem[0];
    }
    else throwError( "Asking to solve a linear system for a 3x3 matrix of rank < 3. The solution is not unique." );
    
    // reduce the 2nd column
    if ( std::fabs(elem[7])>std::fabs(elem[4]) )
    {
        for (int i=1; i<N; ++i) SWAP_VARIABLES( elem[3+i], elem[6+i], T );
        SWAP_VARIABLES( right[1], right[2], T );
        SWAP_VARIABLES( order[1], order[2], int );
    }   
    if ( std::fabs(elem[4])>tolerance )
    {
        elem[8] -= elem[2] * elem[7]/elem[4];
        right[2] -= right[0] * elem[7]/elem[4];
    }
    else throwError( "Asking to solve a linear system for a 3x3 matrix of rank < 3. The solution is not unique." );
    
    // reduce the 3rd column
    if ( not std::fabs(elem[8])>tolerance )
        throwError( "Asking to solve a linear system for a 3x3 matrix of rank < 3. The solution is not unique." );
    
    // get the solution
    result[2] = right[2]/elem[8];
    result[1] = (right[1] - elem[5]*result[2]) /elem[4];
    result[0] = (right[0] - elem[1]*result[1] - elem[2]*result[2]) / elem[0];
    for (int i=0; i<N; ++i)
        solution[ order[i] ] = result[i];
}
template <typename T, size_t N> void Square<T,N>::solveLinearSystemN(T *right, T *solution)
{
    double data[N*N], bVec[N];
    for (int i=0; i<N*N; ++i)
        data[i] = elem[i];
    for (int i=0; i<N; ++i)
        bVec[i] = right[i];
    
    gsl_matrix_view m = gsl_matrix_view_array( data, N, N );
    gsl_vector_view b = gsl_vector_view_array( bVec, N );
    gsl_vector *x = gsl_vector_alloc( N );
    int s;
    gsl_permutation * p = gsl_permutation_alloc( N );
    gsl_linalg_LU_decomp( &m.matrix, p, &s );
    gsl_linalg_LU_solve( &m.matrix, p, &b.vector, x );
    
    for (int i=0; i<N; ++i)
        solution[i] = T( gsl_vector_get( x, i ) );
    
    gsl_permutation_free (p);
    gsl_vector_free (x);
}

/* Computes the eigenvalues and eigenvectors of the matrix. The eigenvalues are sorted ascending. It return only the first n-1 eigenvectors, the last eigenvector can be easily found from the n-1 eigenvectors since must be orthogonal to those. */
template <typename T, size_t N>
template <typename COMPLEX>
void Square<T,N>::eigenvalues(COMPLEX *eigenvalues, COMPLEX *eigenvectors, bool allEigenvectors)
{
    // scale the matrix elements
    T const scale = this->scale();  // returns the sum of the absolute values of the elements
    Square<T,N> tempMat = (*this) / scale;
    
    tempMat.eigenvaluesN( eigenvalues, eigenvectors, allEigenvectors);
    for (size_t i=0; i<N; ++i)
        eigenvalues[i] *= scale;
}


template <typename T, size_t N>
template <typename COMPLEX>
void Square<T,N>::eigenvaluesN(COMPLEX *eigenvalues, COMPLEX *eigenvectors, bool allEigenvectors)
{
    // alocate memory and computational variables
    gsl_vector_complex *eval = gsl_vector_complex_alloc(N);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(N,N);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(N);
    
    
    double data[N*N];
    for (int i1=0; i1<N; ++i1)
        for (int i2=0; i2<N; ++i2)
            data[i1*N+i2] = double( (*this)(i1,i2) );
    
    // compute the eigenvalues and eigenvectors
    gsl_matrix_view m  = gsl_matrix_view_array( data, N, N );
    gsl_eigen_nonsymmv( &m.matrix, eval, evec, w );
    gsl_eigen_nonsymmv_sort( eval, evec, GSL_EIGEN_SORT_VAL_ASC );
    
    // save the results to the input arrays
    for (int i1=0; i1<N; ++i1)
        eigenvalues[i1] = COMPLEX( gsl_vector_complex_get( eval, i1 ) );
    int eigenvecSize = allEigenvectors ? N : (N-1);
    if ( eigenvectors!=NULL )
        for (int i1=0; i1<eigenvecSize; ++i1)
            for (int i2=0; i2<N; ++i2)
                eigenvectors[N*i1+i2] = COMPLEX( gsl_matrix_complex_get( evec, i2, i1 ) );
    
    // free memory
    gsl_eigen_nonsymmv_free (w);
    gsl_vector_complex_free (eval);
    gsl_matrix_complex_free (evec);
}










//! Functions for the 'template <typename T, size_t N> Symmetric<T,N>' class
// constructors
template <typename T, size_t N> Symmetric<T,N>::Symmetric()
{
    if ( N>3 or N<=1 ) throwError( "The class 'MATRIX::Symmetric' has been implemented only for 2x2 and 3x3 symmetric matrices." );
    for (size_t i=0; i<size; ++i)
        elem[i] = T(0.);
}
template <typename T, size_t N> Symmetric<T,N>::Symmetric(T *elements)
{
    if ( N>3 or N<=1 ) throwError( "The class 'MATRIX::Symmetric' has been implemented only for 2x2 and 3x3 symmetric matrices." );
    for (size_t i=0; i<size; ++i)
        elem[i] = elements[i];
}
template <typename T, size_t N> void Symmetric<T,N>::data(T *elements)
{
    for (size_t i=0; i<size; ++i)
        elem[i] = elements[i];
}


// overloading of operators
template <typename T, size_t N> bool Symmetric<T,N>::operator ==(Symmetric<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        if ( elem[i]!=right[i] )
            return false;
    return true;
}
template <typename T, size_t N> bool Symmetric<T,N>::operator !=(Symmetric<T,N> &right)
{
    return not (*this)==right;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator +(Symmetric<T,N> &right)
{
    Symmetric<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] + right[i];
    return temp;
}
template <typename T, size_t N> Symmetric<T,N>& Symmetric<T,N>::operator +=(Symmetric<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] += right[i];
    return *this;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator -(Symmetric<T,N> &right)
{
    Symmetric<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] - right[i];
    return temp;
}
template <typename T, size_t N> Symmetric<T,N>& Symmetric<T,N>::operator -=(Symmetric<T,N> &right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] -= right[i];
    return *this;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator *(Symmetric<T,N> &right)
{
    Symmetric<T,N> temp;
    for (size_t i=0; i<N; ++i)
        for (size_t j=i; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                temp(i,j) += (*this)(i,k) * right(k,j);
    return temp;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator *=(Symmetric<T,N> &right)
{
    return (*this)*right;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator *(T const right)
{
    Symmetric<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] * right;
    return temp;
}
template <typename T, size_t N> Symmetric<T,N>& Symmetric<T,N>::operator *=(T const right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] *= right;
    return *this;
}
template <typename T, size_t N> Symmetric<T,N> Symmetric<T,N>::operator /(T const right)
{
    Symmetric<T,N> temp;
    for (size_t i=0; i<size; ++i)
        temp[i] = elem[i] / right;
    return temp;
}
template <typename T, size_t N> Symmetric<T,N>& Symmetric<T,N>::operator /=(T const right)
{
    for (size_t i=0; i<size; ++i)
        elem[i] /= right;
    return *this;
}

template <typename T, size_t N> T& Symmetric<T,N>::operator [](size_t const i)
{
    BOOST_ASSERT( i<size );
    return elem[i];
}
template <typename T, size_t N> T& Symmetric<T,N>::operator ()(size_t const i, size_t const j)
{
    BOOST_ASSERT( i<N and j<N );
    if (N==2)
        return (i<j) ? elem[i*N+j] : elem[j*N+i];
    if (N==3)
    {
        size_t k = (i<j) ? i : j;
        size_t l = (i>j) ? i : j;
        switch (k)
        {
            case 0: return elem[l]; break;
            case 1: return elem[3+l-1]; break;
            case 2: return elem[5];
        }
    }
}


// Matrix operations
/* Computes the determinant of the matrix. */
template <typename T, size_t N> T Symmetric<T,N>::det()
{
    if (N==2)
        return elem[0]*elem[2] - elem[1]*elem[1];
    else if (N==3)
        return elem[0]*elem[3]*elem[5] + T(2.)*elem[1]*elem[2]*elem[4] - elem[0]*elem[4]*elem[4] - elem[3]*elem[2]*elem[2] - elem[5]*elem[1]*elem[1];
    else
        throwError( "The 'MATRIX::Symetric::det()' function is implemented only for 2x2 and 3x3 matrices." );
}

/* Computes the trace of the matrix. */
template <typename T, size_t N> T Symmetric<T,N>::trace()
{
    if (N==2)
        return elem[0]+elem[2];
    else if (N==3)
        return elem[0]+elem[3]+elem[5];
    else
    {
        T temp = T(0.);
        for (size_t i=0; i<N; ++i)
            temp += (*this)(i,i);
        return temp;
    }
}

template <typename T, size_t N> int Symmetric<T,N>::rank()
{
    T tempData[N*N];
    for (int i=0; i<N; ++i)
        for (int j=0; j<N; ++j)
            tempData[i*N+j] = (*this)(i,j);
    Square<T,N> tempM( tempData );
    return tempM.rank();
}

template <typename T, size_t N> void Symmetric<T,N>::solveLinearSystem(T *right, T *solution)
{
    T tempData[N*N];
    for (int i=0; i<N; ++i)
        for (int j=0; j<N; ++j)
            tempData[i*N+j] = (*this)(i,j);
    Square<T,N> tempM( tempData );
    tempM.solveLinearSystem( right, solution );
}

/* Computes the sum of the absolute values of the matrix entries. */
template <typename T, size_t N> T Symmetric<T,N>::scale()
{
    T temp = T(0.);
    for (size_t i=0; i<size; ++i)
        temp += std::fabs( elem[i] );
    return temp;
}


/* Computes the eigenvalues and eigenvectors of the matrix. The eigenvalues are sorted ascending. It return only the first n-1 eigenvectors, the last eigenvector can be easily found from the n-1 eigenvectors since must be orthogonal to those. */
template <typename T, size_t N> void Symmetric<T,N>::eigenvalues(T *eigenvalues, T *eigenvectors, bool allEigenvectors)
{
    // scale the matrix elements
    T const scale = this->scale();  // returns the sum of the absolute values of the elements
    Symmetric<T,N> tempMat = (*this) / scale;
    
    if (N==2)
        tempMat.eigenvalues2( eigenvalues, eigenvectors, allEigenvectors);
    else if (N==3)
        tempMat.eigenvalues3( eigenvalues, eigenvectors, allEigenvectors);
    else
        tempMat.eigenvaluesN( eigenvalues, eigenvectors, allEigenvectors);
        
    for (size_t i=0; i<N; ++i)
        eigenvalues[i] *= scale;
}
template <typename T, size_t N> void Symmetric<T,N>::eigenvalues2(T *eigenvalues, T *eigenvectors, bool allEigenvectors)
{
    // find the eigenvalues
    T const temp1 = elem[0] + elem[2];
    T const temp2 = elem[0] - elem[2];
    T temp3 = temp2*temp2 + T(4.)*elem[1]*elem[1];
    temp3 = sqrt( temp3 );
    eigenvalues[0] = (temp1-temp3) / T(2.);
    eigenvalues[1] = (temp1+temp3) / T(2.);
    
    
    // check that the precision of the variables does not influence the values of the matrix eigenvalues
    Symmetric<T,N> tempM(elem);
    bool precisionOK = true;
    for (size_t i=0; i<N; ++i)
    {
        tempM[0] = elem[0] - eigenvalues[i]; tempM[2] = elem[2] - eigenvalues[i];
        if ( std::fabs(tempM.det())>largeTolerance ) precisionOK = false;
    }
    if ( not precisionOK )
    {
        this->eigenvaluesN( eigenvalues, eigenvectors, allEigenvectors );
        return;
    }
    
    
    //find the eigenvectors
    if ( eigenvectors==NULL ) return;
    T const temp4 = eigenvalues[0] - elem[0];
    temp3 = elem[1]*elem[1] + temp4*temp4;
    temp3 = std::sqrt( temp3 );
    if ( temp3<tolerance )  //both elements are zero (since their square is 0)
    {
        eigenvectors[0] = T(1.);
        eigenvectors[1] = T(0.);
    }
    else
    {
        eigenvectors[0] = -elem[1] / temp3;
        eigenvectors[1] = temp4 / temp3;
    }
    
    // the second eigenvector is orthogonal on the first
    if ( not allEigenvectors ) return;
    eigenvectors[2] = -eigenvectors[1];
    eigenvectors[3] = eigenvectors[0];
}
/* See: http://en.wikipedia.org/wiki/Eigenvalue_algorithm#Eigenvalues_of_2.C3.972_matrices */
template <typename T, size_t N> void Symmetric<T,N>::eigenvalues3(T *eigenvalues, T *eigenvectors, bool allEigenvectors)
{
    // find the eigenvalues
    T const m = this->trace() / T(3.);
    Symmetric<T,N> temp = *this;
    temp[0] -= m; temp[3] -= m; temp[5] -= m;
    T const q = temp.det() / T(2.);
    
    T p = temp[0]*temp[0] + temp[3]*temp[3] + temp[5]*temp[5] + T(2.)*(temp[1]*temp[1] + temp[2]*temp[2] + temp[4]*temp[4]);
    p /= T(6.);
    T p1 = std::sqrt(p);
    T p2 = p * p1;
    
    T phi = T(0.);
    
    if ( std::fabs(q)<=std::fabs(p2) and p>tolerance )
        phi = std::acos( q/p2 ) / T(3.);
    else if ( q<p2 ) phi = PI / T(3.);
    if ( phi<0. ) phi += PI / T(3.);
    T const cosPhi = std::cos( phi );
    T const sinPhi = std::sin( phi ) * sqrt( T(3.) );

    eigenvalues[2] = m + T(2.) * p1 * cosPhi;
    eigenvalues[1] = m - p1 * ( cosPhi + sinPhi );
    eigenvalues[0] = m - p1 * ( cosPhi - sinPhi );
    
    
    // // check that the precision of the variables does not influence the values of the matrix eigenvalues
    // Symmetric<T,N> tempM(elem);
    // bool precisionOK = true;
    // for (int i=0; i<N; ++i)
    // {
        // tempM[0] = elem[0] - eigenvalues[i]; tempM[3] = elem[3] - eigenvalues[i]; tempM[5] = elem[5] - eigenvalues[i];
        // if ( std::fabs(tempM.det())>largeTolerance ) precisionOK = false;
    // }
    // if ( not precisionOK )
    // {
        // std::cout << "Wrong precision. Calling gsl routine.";
        // this->eigenvaluesN( eigenvalues, eigenvectors, allEigenvectors );
        // return;
    // }
    
    
    // order the eigenvalues ascending
    if ( eigenvalues[0]>eigenvalues[1] )
        SWAP_VARIABLES( eigenvalues[0], eigenvalues[1], T );
    if ( eigenvalues[1]>eigenvalues[2] )
        SWAP_VARIABLES( eigenvalues[1], eigenvalues[2], T );
    if ( eigenvalues[0]>eigenvalues[1] )
        SWAP_VARIABLES( eigenvalues[0], eigenvalues[1], T );
    
    
    //find the eigenvectors
    if ( eigenvectors==NULL ) return;
    int eigenvecSize = allEigenvectors ? N : (N-1);
    Symmetric<T,N> tempM(elem);
    for (int i=0; i<eigenvecSize; ++i)
    {
        tempM[0] = elem[0] - eigenvalues[i]; tempM[3] = elem[3] - eigenvalues[i]; tempM[5] = elem[5] - eigenvalues[i];
        
        eigenvectors[i*3+0] = tempM[2]*tempM[3] - tempM[1]*tempM[4];
        eigenvectors[i*3+1] = tempM[0]*tempM[4] - tempM[1]*tempM[2];
        eigenvectors[i*3+2] = tempM[1]*tempM[1] - tempM[0]*tempM[3];
        T temp3 = eigenvectors[i*3]*eigenvectors[i*3] + eigenvectors[i*3+1]*eigenvectors[i*3+1] + eigenvectors[i*3+2]*eigenvectors[i*3+2];
        temp3 = std::sqrt( temp3 );
        if ( temp3 > largeTolerance )
        {
            eigenvectors[i*3] /= temp3; eigenvectors[i*3+1] /= temp3; eigenvectors[i*3+2] /= temp3;
            continue;
        }
        
        eigenvectors[i*3+0] = tempM[1]*tempM[5] - tempM[2]*tempM[4];
        eigenvectors[i*3+1] = tempM[2]*tempM[2] - tempM[0]*tempM[5];
        eigenvectors[i*3+2] = tempM[0]*tempM[4] - tempM[1]*tempM[2];
        temp3 = eigenvectors[i*3]*eigenvectors[i*3] + eigenvectors[i*3+1]*eigenvectors[i*3+1] + eigenvectors[i*3+2]*eigenvectors[i*3+2];
        temp3 = std::sqrt( temp3 );
        if ( temp3 > largeTolerance )
        {
            eigenvectors[i*3] /= temp3; eigenvectors[i*3+1] /= temp3; eigenvectors[i*3+2] /= temp3;
            continue;
        }
        
        eigenvectors[i*3+0] = tempM[3]*tempM[5] - tempM[4]*tempM[4];
        eigenvectors[i*3+1] = tempM[2]*tempM[4] - tempM[1]*tempM[5];
        eigenvectors[i*3+2] = tempM[1]*tempM[2] - tempM[2]*tempM[3];
        temp3 = eigenvectors[i*3]*eigenvectors[i*3] + eigenvectors[i*3+1]*eigenvectors[i*3+1] + eigenvectors[i*3+2]*eigenvectors[i*3+2];
        temp3 = std::sqrt( temp3 );
        if ( temp3 > largeTolerance )
        {
            eigenvectors[i*3] /= temp3; eigenvectors[i*3+1] /= temp3; eigenvectors[i*3+2] /= temp3;
            continue;
        }
        
        //special case -> send it to the slower function
        this->eigenvaluesN( eigenvalues, eigenvectors, allEigenvectors);
        break;
    }
}
template <typename T, size_t N> void Symmetric<T,N>::eigenvaluesN(T *eigenvalues, T *eigenvectors, bool allEigenvectors)
{
//   cout << "\t Called gsl routine.\n";
   // alocate memory and computational variables
    gsl_vector *eval = gsl_vector_alloc (N);
    gsl_matrix *evec = gsl_matrix_alloc (N,N);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (N);

    double data[N*N];
    for (size_t i1=0; i1<N; ++i1)
        for (size_t i2=0; i2<N; ++i2)
            data[i1*N+i2] = double( (*this)(i1,i2) );
        
    // compute the eigenvalues and eigenvectors
    gsl_matrix_view m  = gsl_matrix_view_array( data, N, N );
    gsl_eigen_symmv( &m.matrix, eval, evec, w );
    gsl_eigen_symmv_sort( eval, evec, GSL_EIGEN_SORT_VAL_ASC );

    // save the results to the input arrays
    for (size_t i1=0; i1<N; ++i1)
        eigenvalues[i1] = T( gsl_vector_get( eval, i1 ) );
    int eigenvecSize = allEigenvectors ? N : (N-1);
    if ( eigenvectors!=NULL )
        for (int i1=0; i1<eigenvecSize; ++i1)
            for (size_t i2=0; i2<N; ++i2)
                eigenvectors[N*i1+i2] = T( gsl_matrix_get( evec, i2, i1 ) );
   
    // free memory
    gsl_eigen_symmv_free (w);
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
}



} //end namespace MATRIX
#endif


