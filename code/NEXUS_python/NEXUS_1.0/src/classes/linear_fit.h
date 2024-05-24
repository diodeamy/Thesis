
#ifndef FIT_HEADER
#define FIT_HEADER

#include <cmath>
#include <boost/assert.hpp>
#include <gsl/gsl_multifit.h>

#include <miscellaneous.h>
#include <arrayProperties.h>
#include <array.h>


namespace LEAST_SQUARE_METHOD
{
/* Returns the sign of the entries of an array. */
template <typename T, typename T_INT> int signEntries(T *x, T_INT const size);
static int const POSITIVE_ENTRIES = 1;
static int const NEGATIVE_ENTRIES = -1;
static int const MIXED_SIGN_ENTRIES = 0;



/*! This class fits a linear function to a set of 'x,y' data with the form:
        y = c0 + c1 * x
NOTE: Fitting result taken from: http://mathworld.wolfram.com/LeastSquaresFitting.html
*/
template< typename T >
class FitLinearFunction
{
    public:
    T c[2];                        // the coefficients of the fit
    T residuals;                   // gives the sum of the square residuals
    
    void dataToFit(T *x, T *y, int const size);    // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<2 and n>=0 ); return c[n];}
    void fitParameters(T * param1, T * param2) const // this function retunrs both 'c0' and 'c1'
    { *param1 = c[0]; *param2 = c[1];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t  minDataSize() const { return 3;}       // returns the minimum number of data needed for the fit
    size_t  noFitParamaters() const { return 2;}   // returns the number of fit parameters
    T yValue(T const x) const {return c[0]+c[1]*x;}// returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const {return c[1];}  // returns the dy/dx value corresponding to 'x' from the fit function
};

/* This function reads the data and does the actual fitting. */
template< typename T > void FitLinearFunction<T>::dataToFit(T *x, T *y, int const size)
{
    if ( size<=2 )
        throwError( "The data array to be fitted in class 'LinearFunctionFit' must have more than 2 elements." );
    
    T sum_x = T(0);
    T sum_x2 = T(0);
    T sum_y = T(0);
    T sum_y2 = T(0);
    T sum_xy = T(0);
    
    for (int i=0; i<size; ++i)
    {
        sum_x += x[i];
        sum_x2 += x[i] * x[i];
        sum_y += y[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    
    T ss_xx = sum_x2 - sum_x*sum_x / size;
    T ss_xy = sum_xy - sum_x*sum_y / size;
    
    c[1] = ss_xy / ss_xx;            // the slope
    c[0] = sum_y/size - c[1] * sum_x/size;   // the constant term
    residuals = sum_y2 + size*c[0]*c[0] + c[1]*c[1]*sum_x2 - 2.*c[0]*sum_y - 2.*c[1]*sum_xy + 2.*c[0]*c[1]*sum_x; //the sum of the residuals square 
}




/*! This class fits a linear function to a set of 'x,y' data and computes extended statistical properties for the fit. The fit is given by:
        y = c[0] + c[1] * x
NOTE: Fitting result taken from: http://mathworld.wolfram.com/LeastSquaresFitting.html
*/
template< typename T >
class FitLinearFunction_extendedProperties
{
    public:
    T c[2];                        // the coefficients of the fit
    T residuals;                   // gives the sum of the square residuals
    T sigma, err_c[2], r2;         // errors and goodness of the fit
    
    FitLinearFunction_extendedProperties();
    void dataToFit(T *x, T *y, int const size);    // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<2 and n>=0 ); return c[n];}
    void fitParameters(T * param1, T * param2) const // this function retunrs both 'c0' and 'c1'
    { *param1 = c[0]; *param2 = c[1];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t  minDataSize() const { return 3;}       // returns the minimum number of data needed for the fit
    size_t  noFitParamaters() const { return 2;}   // returns the number of fit parameters
    T yValue(T const x) const {return c[0]+c[1]*x;}// returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const {return c[1];}  // returns the dy/dx value corresponding to 'x' from the fit function

    T standardDeviation() const {return sigma;}    // return standard deviation
    T error_0() const {return err_c[0];}           // returns error of 'c0'
    T error_1() const {return err_c[1];}           // returns error of 'c1'
    T R2() const {return r2;}                      // returns 'R2' = goodness of fit
};

template< typename T > FitLinearFunction_extendedProperties<T>::FitLinearFunction_extendedProperties()
{
    c[0] = T(0); c[1] = T(0);
    residuals = T(0);
    sigma = T(0); err_c[0] = T(0); err_c[1] = T(0); r2 = T(0);
}

/* This function reads the data and does the actual fitting. */
template< typename T > void FitLinearFunction_extendedProperties<T>::dataToFit(T *x, T *y, int const size)
{
    if ( size<=2 )
        throwError( "The data array to be fitted in class 'LinearFunctionFit' must have more than 2 elements." );
    
    T sum_x = T(0);
    T sum_x2 = T(0);
    T sum_y = T(0);
    T sum_y2 = T(0);
    T sum_xy = T(0);
    
    for (int i=0; i<size; ++i)
    {
        sum_x += x[i];
        sum_x2 += x[i] * x[i];
        sum_y += y[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    
    T ss_xx = sum_x2 - sum_x*sum_x / size;
    T ss_yy = sum_y2 - sum_y*sum_y / size;
    T ss_xy = sum_xy - sum_x*sum_y / size;
    
    c[1] = ss_xy / ss_xx;            // the slope
    c[0] = sum_y/size - c[1] * sum_x/size;   // the constant term
    residuals = sum_y2 + size*c[0]*c[0] + c[1]*c[1]*sum_x2 - 2.*c[0]*sum_y - 2.*c[1]*sum_xy + 2.*c[0]*c[1]*sum_x; //the sum of the residuals square 
    
    r2 = ss_xy*ss_xy / (ss_xx*ss_yy); // the correlation coefficient - goodness of fit
    sigma = sqrt( (ss_yy-c[1]*ss_xy) / (size-2) ); // standard deviation
    err_c[0] = sigma * sqrt( 1./size + sum_x*sum_x/(ss_xx*size*size) ); // error in determination constant coefficient
    err_c[0] = sigma / sqrt( ss_xx ); //error in slope determination
}





/*! This class fits a polynomial function to a set of 'x,y' data with the form:
        y = sum_{i=0}^{N} c[i] * x^i
*/
template< typename T, size_t const N >
class FitPolynomialFunction
{
    public:
    T c[N+1];                        // the coefficients of the fit
    T residuals;                     // gives the sum of the square residuals
    
    FitPolynomialFunction();
    void dataToFit(T *x, T *y, size_t const size); // this function does the actual fit
    
    T param_n(size_t const n) const;               // this function returns 'cn' wih n=0, 1 or 2
    void fitParameters(T *param, size_t const size) const; // this function retunrs both 'c0' and 'c1'
    T residual() const {return residuals;}         // returns the fit residuals
    size_t  minDataSize() const { return N+1;}     // returns the minimum number of data needed for the fit
    size_t  noFitParamaters() const { return N+1;} // returns the number of fit parameters
    T yValue(T const x) const;                     // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const;                           // returns the dy/dx value corresponding to 'x' from the fit function
};

template< typename T, size_t const N > FitPolynomialFunction<T,N>::FitPolynomialFunction()
{
    for (size_t i=0; i<=N; ++i)
        c[i] = T(0);
    residuals = T(0);
}

/* This function reads the data and does the actual fitting. */
template< typename T, size_t const N > void FitPolynomialFunction<T,N>::dataToFit(T *x, T *y, size_t const size)
{
    if ( size<N+1 )
        throwError( "The data array to be fitted in class 'FitPolynomialFunction' must have at least " , N+1, " elements." );
    
    gsl_matrix *X, *COV;
    gsl_vector *Y, *C;

    X = gsl_matrix_alloc( size, N+1 );
    Y = gsl_vector_alloc( size );
    
    C = gsl_vector_alloc( N+1 );
    COV = gsl_matrix_alloc( N+1, N+1 );
    
    for (size_t i=0; i<size; ++i)
    {
        double temp = 1.;
        for (size_t j=0; j<=N; ++j)
        {
            gsl_matrix_set( X, i, j, temp );
            temp *= x[i];
        }
        gsl_vector_set (Y, i, y[i]);
    }
    
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc( size, N+1 );
    double chisq;
    gsl_multifit_linear( X, Y, C, COV, &chisq, work );
    gsl_multifit_linear_free( work );
    
    for (size_t j=0; j<=N; ++j)
        c[j] = gsl_vector_get( C, j );
    
    gsl_matrix_free (X);
    gsl_vector_free (Y);
    gsl_vector_free (C);
    gsl_matrix_free (COV);
    
    // compute the residuals
    residuals = 0.;
    for (size_t i=0; i<size; ++i)
    {
        T res = y[i];
        T temp = 1.;
        for (size_t j=0; j<=N; ++j)
        {
            res -= c[i] * temp;
            temp *= x[i];
        }
        residuals += res * res;
    }
}

template< typename T, size_t const N > T FitPolynomialFunction<T,N>::param_n(size_t const n) const
{
    BOOST_ASSERT( n<=N and n>=0 );
    return c[n];
}

template< typename T, size_t const N > void FitPolynomialFunction<T,N>::fitParameters(T *param, size_t const size) const
{
    BOOST_ASSERT( size==N+1 );
    for (size_t i=0; i<N; ++i)
        param[i] = c[i];
}

/* Return 'y' value for a given 'x' using fit coefficients. */
template< typename T, size_t const N > T FitPolynomialFunction<T,N>::yValue(T const x) const
{
    T res = 0.;
    T temp = 1.;
    for (size_t j=0; j<=N; ++j)
    {
        res += c[j] * temp;
        temp *= x;
    }
    return res;
}
/* Return 'dy/dx' value for a given 'x' using fit coefficients. */
template< typename T, size_t const N > T FitPolynomialFunction<T,N>::yDerivative(T const x) const
{
    T res = 0.;
    T temp = 1.;
    for (size_t j=1; j<=N; ++j)
    {
        res += c[j]*j * temp;
        temp *= x;
    }
    return res;
}




/*! This class fits a parabola function to a set of 'x,y' data with the form:
        y = c[0] + c[1] * x + c[2] * x*x
 */
template< typename T>
class FitParabolaFunction: public FitPolynomialFunction<T,2>
{
    public:
    T param_0() const {return this->c[0];}            // this function returns 'c0'
    T param_1() const {return this->c[1];}            // this function returns 'c1'
    T param_2() const {return this->c[2];}            // this function returns 'c2'
    T xAtMaximum() const;                             // returns the maximum of the fit function, if any
    T xAtMinimum() const;                             // returns the maximum of the fit function, if any
};


template< typename T > T FitParabolaFunction<T>::xAtMaximum() const
{
    if ( this->c[2]>0. )
    {
        throwWarning( "You are trying to get the maximum of the quadratic function, but it does not have a maximum since the second derivative is positive." );
        return 0.;
    }
    return -0.5 * this->c[1] / this->c[2];
}
template< typename T > T FitParabolaFunction<T>::xAtMinimum() const
{
    if ( this->c[2]<0. )
    {
        throwWarning( "You are trying to get the minimum of the quadratic function, but it does not have a minimum since the second derivative is negative." );
        return 0.;
    }
    return -0.5 * this->c[1] / this->c[2];
}





/*! This function fits an std::exponential to a data set:
        y = c[0] * std::exp( c[1] * x )
NOTE: For the difference between equal and non-equal weight see: http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
*/
template< typename T >
class FitExponentialFunction
{
    public:
    T c[2];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size, bool equalWeight = false); // this function does the actual fit; if 'equalWeight', all data get the same wait
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<2 and n>=0 ); return c[n];}
    void fitParameters(T * param1, T * param2) const // this function retunrs both 'c0' and 'c1'
    { *param1 = c[0]; *param2 = c[1];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t minDataSize() const { return 3;}        // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return 2;}    // returns the number of fit parameters
    T yValue(T const x) const {return c[0]*std::exp(x*c[1]);} // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const                 // returns the 'dy/dx' value corresponding to 'x' from the fit function
    { return c[0]*c[1] * std::exp( x*c[1] );}
    
};

template< typename T > void FitExponentialFunction<T>::dataToFit(T *x, T *y,
                                                                 size_t const size,
                                                                 bool equalWeight)
{
    //check that all the data have the same sign - cannot take log(0)
    if ( MIXED_SIGN_ENTRIES==signEntries(y,size) )
        throwError( "Cannot fit the data with an exponential function since the 'y' values do not have all the same sign. Error in function 'FitExponentialFunction'." );
    T signY = T(1.);
    if ( NEGATIVE_ENTRIES==signEntries(y,size) )
        signY = T(-1.);
    
    T lnY[size];
    for (size_t i=0; i<size; ++i)
        lnY[i] = std::log( signY*y[i] );
    
    if ( not equalWeight )
    {
        FitLinearFunction<T> linearFit;
        linearFit.dataToFit( x, lnY, size );
        
        c[0] = std::exp( linearFit.param_0() );
        c[1] = linearFit.param_1();
    }
    else if ( equalWeight )
    {
        T sum_xy=T(0.), sum_x2y=T(0.), sum_xylnY=T(0.), sum_y=T(0.), sum_ylnY=T(0.);
        for (size_t i=0; i<size; ++i)
        {
            sum_xy    += x[i] * y[i];
            sum_x2y   += x[i]*x[i] * y[i];
            sum_xylnY += x[i] * y[i] * lnY[i];
            sum_y     += y[i];
            sum_ylnY  += y[i] * lnY[i];
        }
        
        c[0] = (sum_x2y*sum_ylnY - sum_xy*sum_xylnY) / (sum_y*sum_x2y - sum_xy*sum_xy);
        c[0] = std::exp( c[0] );
        c[1] = (sum_y*sum_xylnY - sum_xy*sum_ylnY) / (sum_y*sum_x2y - sum_xy*sum_xy);
    }
    
    c[0] *= signY;
    residuals = T(0.);
    for (size_t i=0; i<size; ++i)
    {
        T res = y[i] - c[0] * std::exp( c[1]*x[i] );
        residuals += res * res;
    }
}





/*! This function fits a power law to a data set:
        y = c[0] * x^c[1]
*/
template< typename T >
class FitPowerLawFunction
{
    public:
    T c[2];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size); // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<2 and n>=0 ); return c[n];}
    void fitParameters(T * param1, T * param2) const // this function retunrs both 'c0' and 'c1'
    { *param1 = c[0]; *param2 = c[1];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t minDataSize() const { return 3;}           // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return 2;}       // returns the number of fit parameters
    T yValue(T const x) const {return c[0]*std::pow(x,c[1]);} // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const                 // returns the 'dy/dx' value corresponding to 'x' from the fit function
    { return c[0]*c[1] * std::pow( x, T(c[1]-1.) );}
};

template< typename T > void FitPowerLawFunction<T>::dataToFit(T *x, T *y,
                                                              size_t const size)
{
    //check that all the data have the same sign - cannot take log(0)
    if ( MIXED_SIGN_ENTRIES==signEntries(y,size) )
        throwError( "Cannot fit the data with a power law function since the 'y' values do not have all the same sign. Error in function 'FitPowerLawFunction'." );
    if ( MIXED_SIGN_ENTRIES==signEntries(x,size) or NEGATIVE_ENTRIES==signEntries(x,size) )
        throwError( "Cannot fit the data with a power law function since the 'x' values have negative values too. Error in function 'FitPowerLawFunction'." );
    T signY = T(1.);
    if ( NEGATIVE_ENTRIES==signEntries(y,size) )
        signY = T(-1.);
    
    T lnX[size];
    T lnY[size];
    for (size_t i=0; i<size; ++i)
    {
        lnX[i] = std::log( x[i] );
        lnY[i] = std::log( signY*y[i] );
    }
    
    FitLinearFunction<T> linearFit;
    linearFit.dataToFit( lnX, lnY, size );
        
    c[0] = signY * std::exp( linearFit.param_0() );
    c[1] = linearFit.param_1();
    
    residuals = T(0.);
    for (size_t i=0; i<size; ++i)
    {
        T res = y[i] - c[0] * std::pow( x[i], c[1] );
        residuals += res * res;
    }
}






/*! This function fits a logarithm function to a data set:
        y = c[0] + c[1] * ln(x)
 */
template< typename T >
class FitLogarithmFunction
{
    public:
    T c[2];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size); // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<2 and n>=0 ); return c[n];}
    void fitParameters(T * param1, T * param2) const // this function retunrs both 'c0' and 'c1'
    { *param1 = c[0]; *param2 = c[1];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t minDataSize() const { return 3;}           // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return 2;}       // returns the number of fit parameters
    T yValue(T const x) const {return c[0]+c[1]*std::log(x);} // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const { return c[1]/x;} // returns the 'dy/dx' value corresponding to 'x' from the fit function
};

template< typename T > void FitLogarithmFunction<T>::dataToFit(T *x, T *y,
                                                               size_t const size)
{
    //check that all the data have the same sign - cannot take log(0)
    if ( MIXED_SIGN_ENTRIES==signEntries(x,size) or NEGATIVE_ENTRIES==signEntries(x,size) )
        throwError( "Cannot fit the data with a log(x) function since the 'x' values have negative values too. Error in function 'FitLogarithmFunction'." );
    
    T lnX[size];
    for (size_t i=0; i<size; ++i)
        lnX[i] = std::log( x[i] );
    
    FitLinearFunction<T> linearFit;
    linearFit.dataToFit( lnX, y, size );
        
    c[0] = linearFit.param_0();
    c[1] = linearFit.param_1();
    
    residuals = T(0.);
    for (size_t i=0; i<size; ++i)
    {
        T res = y[i] - c[0] - c[1] * lnX[i];
        residuals += res * res;
    }
}






/*! This function fits a gaussian function to a data set:
        y = c[0] * exp(c[1] * (x-c[2])**2)
 */
template< typename T >
class FitGaussianFunction
{
    public:
    T c[3];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size); // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_x0() const {return c[2];}              // this function returns 'x0'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<3 and n>=0 ); return c[n];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t minDataSize() const { return 4;}        // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return 3;}    // returns the number of fit parameters
    T yValue(T const x) const {return c[0] * std::exp(c[1]*(x-c[2])*(x-c[2]));} // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const                 // returns the 'dy/dx' value corresponding to 'x' from the fit function
    {return 2.*c[0]*c[1]*(x-c[2]) * std::exp(c[1]*(x-c[2])*(x-c[2]));}
};

template< typename T > void FitGaussianFunction<T>::dataToFit(T *x, T *y,
                                                              size_t const size)
{
    //check that all the data have the same sign - cannot take log(0)
    if ( MIXED_SIGN_ENTRIES==signEntries(y,size) )
        throwError( "Cannot fit the data with a gaussian function since the 'y' values have both positive and negative values (cannot take log(negative value)). Error in function 'FitGaussianFunction'." );
    T signY = T(1.);
    if ( NEGATIVE_ENTRIES==signEntries(y,size) )
        signY = T(-1.);
    
    T lnY[size];
    for (size_t i=0; i<size; ++i)
        lnY[i] = std::log( signY*y[i] );
    
    FitParabolaFunction<T> quadraticFit;
    quadraticFit.dataToFit( x, lnY, size );
    
    c[1] = quadraticFit.param_n(2);
    c[2] = quadraticFit.param_n(1) / (-2.*c[1]);
    c[0] = signY * std::exp( quadraticFit.param_n(0) - c[1] * c[2]*c[2] );
    
    residuals = T(0.);
    for (size_t i=0; i<size; ++i)
    {
        T res = y[i] - this->yValue( x[i] );
        residuals += res * res;
    }
}






/*! This function fits a lognormal function to a data set:
        y = c[0] * exp(c[1] * (log10(x)-c[2])**2)
 */
template< typename T >
class FitLogNormalFunction
{
    public:
    T c[3];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size); // this function does the actual fit
    
    T param_0() const {return c[0];}               // this function returns 'c0'
    T param_1() const {return c[1];}               // this function returns 'c1'
    T param_x0() const {return c[2];}              // this function returns 'x0'
    T param_n(size_t const n) const                // return fit parameter 'n'
    { BOOST_ASSERT( n<3 and n>=0 ); return c[n];}
    T residual() const {return residuals;}         // returns the fit residuals
    size_t minDataSize() const { return 4;}        // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return 3;}    // returns the number of fit parameters
    T yValue(T const x) const                      // returns the y value corresponding to 'x' from the fit function
    {return c[0] * std::exp(c[1]*(std::log10(x)-c[2])*(std::log10(x)-c[2]));}
    T yDerivative(T const x) const                 // returns the 'dy/dx' value corresponding to 'x' from the fit function
    {return 2.*c[0]*c[1]*(std::log10(x)-c[2])/x * std::exp(c[1]*(std::log10(x)-c[2])*(std::log10(x)-c[2]));}
};

template< typename T > void FitLogNormalFunction<T>::dataToFit(T *x, T *y,
                                                               size_t const size)
{
    //check that all the data have the same sign - cannot take log(0)
    if ( MIXED_SIGN_ENTRIES==signEntries(y,size) )
        throwError( "Cannot fit the data with a lognormal distribution since the 'y' values have both positive and negative values (cannot take log(negative value)). Error in function 'FitGaussianFunction'." );
    if ( MIXED_SIGN_ENTRIES==signEntries(x,size) or NEGATIVE_ENTRIES==signEntries(x,size) )
        throwError( "Cannot fit the data with a lognormal distribution since the 'x' values have negative values too. Error in function 'FitLogNormalFunction'." );
    
    T logX[size];
    for (size_t i=0; i<size; ++i)
        logX[i] = std::log10( x[i] );
    
    FitGaussianFunction<T> gaussianFit;
    gaussianFit.dataToFit( logX, y, size );
    
    for (size_t i=0; i<3; ++i)
        c[i] = gaussianFit.param_n(i);
    
    residuals = gaussianFit.residual();
}




/*! This function fits a power law to a data set:
        y = c[0] * x^c[1]
*/
template< typename T, typename Fit>
class FitFunctionWithOffset
{
    public:
    Fit f;
    T offsetValue;
    T c[2];
    T residuals;
    
    void dataToFit(T *x, T *y, size_t const size);         // this function does the actual fit
    
    T param_n(size_t const n) const { return f.param_n(n);}// return fit parameter 'n'
    T residual() const {return f.residual();}              // returns the fit residuals
    size_t minDataSize() const { return f.minDataSize();}  // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return f.noFitParamaters();}// returns the number of fit parameters
    T yValue(T const x) const {return offsetValue+f.yValue(x);}// returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const { return f.yDerivative(x);} // returns the 'dy/dx' value corresponding to 'x' from the fit function
    void setOffset(T const offset) { offsetValue = offset;}    // sets the offset value
};
template< typename T, typename Fit>
void FitFunctionWithOffset<T,Fit>::dataToFit(T *x, T *y,
                                             size_t const size)
{
    T tempY[size];
    for (size_t i=0; i<size; ++i)
        tempY[i] = y[i] - offsetValue;
    
    f.dataToFit( x, tempY, size );
}








/*! Use two functions to fit the data on different data intervals. Decide the switch point between the two fit function by minimizing the residuals.
*/
template< typename T, typename Fit_1, typename Fit_2 >
class FitMultipleFunctions
{
    public:
    Fit_1 f1;
    Fit_2 f2;
    T residuals;
    T xBoundary;
    int boundaryIndex;
    T dataBoundary;
    
    void dataToFit(T *x, T *y, size_t const size);             // this function does the actual fit
    
    T param_f1_n(size_t const n) const {return f1.param_n(n);}// this function returns fit parameter 'n' of function 1
    T param_f2_n(size_t const n) const {return f2.param_n(n);}// this function returns fit parameter 'n' of function 2
    T param_n(size_t const n) const;                          // return fit parameter 'n'
    T residual() const {return residuals;}                    // returns the fit residuals
    size_t minDataSize() const { return f1.minDataSize()+f2.minDataSize()-1;} // returns the minimum number of data needed for the fit
    size_t noFitParamaters() const { return f1.noFitParamaters()+f2.noFitParamaters();} // returns the number of fit parameters
    T yValue(T const x) const;                                // returns the y value corresponding to 'x' from the fit function
    T yDerivative(T const x) const;                           // returns the dy/dx value corresponding to 'x' from the fit function
    T boundary() const {return xBoundary;}                    // returns the intersection point of the 2 functions
};

/* Compute the best fit. */
template< typename T, typename Fit_1, typename Fit_2 >
void FitMultipleFunctions<T,Fit_1,Fit_2>::dataToFit(T *x, T *y, size_t const size)
{
    size_t const n1 = f1.minDataSize();
    size_t const n2 = f2.minDataSize();
    size_t const noFits = size+1 - (n1+n2);
    if ( size<=n1+n2)
        throwError( "The data array to be fitted in class 'FitMultipleFunctions' must have at least " , n1+n2, " elements." );
    T tempResiduals[noFits];
    
    for (size_t i=n1; i<size+1-n2; ++i)
    {
        f1.dataToFit( x, y, i );
        f2.dataToFit( &(x[i-1]), &(y[i-1]), size+1-i );
        
        tempResiduals[i-n1] = f1.residual() + f2.residual();
    }
    
    ArrayProperties<T,size_t> arrayProp( tempResiduals, noFits );
    boundaryIndex = arrayProp.indexMinimum();
    dataBoundary = x[boundaryIndex];
    
    // do a fit again to find the parameters of best fit
    f1.dataToFit( x, y, boundaryIndex );
    f2.dataToFit( &(x[boundaryIndex-1]), &(y[boundaryIndex-1]), size+1-boundaryIndex );
    residuals = f1.residual() + f2.residual();
    
    
    // find the value of 'xBoundary'
    T const factor = 5.;
    T X = dataBoundary;
    T dX = ( x[boundaryIndex+1] - x[boundaryIndex-1] ) / factor;
    int direction = 0;
    while ( dX/X>1.e-4 )
    {
        T tempValue = f1.yValue(X) - f2.yValue(X);
        T tempDerivative = f1.yDerivative(X) - f2.yDerivative(X);
        if ( (tempValue>0. and tempDerivative>0.) or (tempValue<0. and tempDerivative<0.) )
        {
            if ( direction==1 )  // if previous direction was positive, we must be close to the intersection point of the two functions - so decrease iteration step
                dX /= factor;
            X -= dX;
            direction = -1;
        }
        else if ( (tempValue>0. and tempDerivative<0.) or (tempValue<0. and tempDerivative>0.) )
        {
             if ( direction==-1 ) 
                 dX /= factor;
            X += dX;
            direction = 1;
        }
        else
            break;
    }
    xBoundary = X;
}

/* Returns the fit parameter 'n' of the fit. */
template< typename T, typename Fit_1, typename Fit_2 >
T FitMultipleFunctions<T,Fit_1,Fit_2>::param_n(size_t const n) const
{
    BOOST_ASSERT( n>=0 and n<this->noFitParamaters() );
    if ( n<f1.noFitParamaters() )
        return f1.param_n( n );
    else
        return f2.param_n( n-f1.noFitParamaters() );
}

/* Return the value of 'y' for a given 'x' from the fit parameters. */
template< typename T, typename Fit_1, typename Fit_2 >
T FitMultipleFunctions<T,Fit_1,Fit_2>::yValue(T const x) const
{
    if ( x<=xBoundary )
        return f1.yValue( x );
    else
        return f2.yValue( x );
}

/* Return the value of 'dy/dx' for a given 'x' from the fit parameters. */
template< typename T, typename Fit_1, typename Fit_2 >
T FitMultipleFunctions<T,Fit_1,Fit_2>::yDerivative(T const x) const
{
    if ( x<=xBoundary )
        return f1.yDerivative( x );
    else
        return f2.yDerivative( x );
}



/* Checks the sign of all the array elements and returns:
          POSITIVE_ENTRIES - if all entries are positive
          NEGATIVE_ENTRIES - if all entries are negative
          MIXED_SIGN_ENTRIES - if there are both + and - entries
*/
template <typename T, typename T_INT> int signEntries(T *x, T_INT const size)
{
    bool positive = false;
    if ( x[0]>T(0.) ) positive = true;
    if ( positive )
    {
         for (T_INT i=1; i<size; ++i)
             if ( x[i]<T(0.) )
                 return MIXED_SIGN_ENTRIES;
         return POSITIVE_ENTRIES;
    }
    else
    {
         for (T_INT i=1; i<size; ++i)
             if ( x[i]>T(0.) )
                 return MIXED_SIGN_ENTRIES;
         return NEGATIVE_ENTRIES;
    }
}




}   //! end namespace
#endif


///*! This function fits a power law to a data set:
//        y = c[0] + c[1] * x^c[2]
// */
//template< typename T, typename Func >
//class FitFunctionWithOffset
//{
//    public:
//        Func<T> f;
//        T c[3];
//        T residuals;
//    
//        void dataToFit(T *x, T *y, size_t const size, checkSign = true); // this function does the actual fit
//    
//        T param_n(size_t const n) const { return f.param_n(c);}// return fit parameter 'n'
//        T residual() const {return residuals;}                 // returns the fit residuals
//        size_t minDataSize() const { return f.minDataSize();}  // returns the minimum number of data needed for the fit
//        size_t noFitParamaters() const { return f.noFitParamaters();}  // returns the number of fit parameters
//        T yValue(T const x) const {return f.yValue(x);}        // returns the y value corresponding to 'x' from the fit function
//        T yDerivative(T const x) const { return f.yDerivative(x);} // returns the 'dy/dx' value corresponding to 'x' from the fit function
//};
//
//template< typename T, typename Func > 
//void FitFunctionWithOffset<T,Func>::dataToFit(T *x, T *y,
//                                              size_t const size,
//                                              bool checkSign)
//{
//    //check that all the data have the same sign - cannot take log(0)
//    if ( checkSign )
//    {
//        if ( MIXED_SIGN_ENTRIES==signEntries(y,size) )
//            throwError( "Cannot fit the data since the 'y' values do not have all the same sign. Error in function 'FitFunctionWithOffset'." );
//        if ( MIXED_SIGN_ENTRIES==signEntries(x,size) or NEGATIVE_ENTRIES==signEntries(x,size) )
//            throwError( "Cannot fit the data since the 'x' values have negative values too. Error in function 'FitFunctionWithOffset'." );
//    }
//    T signY = T(1.);
//    if ( NEGATIVE_ENTRIES==signEntries(y,size) )
//        signY = T(-1.);
//    Array<T,1> tempY(size);  //store only positive Y values
//    for (size_t i=0; i<size; ++i)
//        tempY[i] = signY * y[i];
//    
//    
//    // do the actual fit
//    size_t const No = 10;
//    T tempResiduals[No];
//    T A, dA;
//    Array<T,1> tempY2(size);
//    // keep track of the min and max values of y
//    ArrayProperties<T,size_t> prop( tempY.ptrData(), size );
//    T minY = prop.valueMinimum() * (1. - 1.e-6);     //use small offset factor to don't get y-A==0 !
//    T maxY = prop.valueMaximum() * (1. + 1.e-6);
//    
//    
//    // if 'checkSign'=true, the y values must have one single sign - must check separately for A>maxY and A<minY
//    if ( checkSign )
//    {
//        // first check for A>maxY
//        A = 1.5*maxY;
//        dA = maxY / No;
//        T minA = minY;
//        bool boundary = false;   // true if the iteration moves to the boundary, used to increase step length
//        T upperMinimum, lowerMinimum; // keep track of the minimum residual for A>maxY and A<minY respectively
//        while ( dA/A>1.e-4 )
//        {
//             tempY2 = tempY - minA;
//             for (size_t i=0; i<No; ++i) // loop over different A values
//             {
//                 tempY2 -= dA;
//                 tempResiduals[i] = f.dataToFit( x, tempY2.ptrData(), size ); // get the residual for each fit
//             }
//             
//             prop.newData( tempResiduals, No );
//             
//             if ( prop.indexMinimum==No-1 )
//             {
//                 A += No/2 * dA;
//                 if ( boundary ) dA *= 2;
//                 boundary = true;
//             }
//             else if ( prop.indexMinimum==0 )
//                 break;       // stop if the minimum is for A<maxY
//             else
//             {
//                 A = minA + prop.indexMinimum*dA;
//                 dA /= T( No/2 );
//                 boundary = false;
//             }
//        }
//    }
//    
//    T lnX[size];
//    T lnY[size];
//    for (size_t i=0; i<size; ++i)
//    {
//        if ( x[i]<=0. or y[i]<=0. ) throwError( "Found values of 0. in the data set to be fitted in function 'FitPowerLawOffsetFunction'. The fit cannot continue since one cannot take logarithm from a non-pozitive number." );
//        lnX[i] = std::log( x[i] );
//        lnY[i] = std::log( y[i] );
//    }
//    
//    
//    
//    
//    while ( (da/a>1.e-4 and a>1.e-7 ) or a< )
//    FitLinearFunction<T> linearFit;
//    linearFit.dataToFit( lnX, lnY, size );
//        
//    c[0] = std::exp( linearFit.param_0() );
//    c[1] = linearFit.param_1();
//    
//    residuals = T(0.);
//    for (size_t i=0; i<size; ++i)
//    {
//        T res = y[i] - c[0] * std::pow( x[i], c[1] );
//        residuals += res * res;
//    }
//}

