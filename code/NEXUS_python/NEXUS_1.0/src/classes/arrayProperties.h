
#ifndef ARRAYPROPERTIES_HEADER
#define ARRAYPROPERTIES_HEADER

#include <cmath>
#include <algorithm>
#include <boost/assert.hpp>

#include <miscellaneous.h>


static int const ORDER_ASCENDING = 1;
static int const ORDER_DESCENDING = -1;
static int const LINEAR_BIN = 0;
static int const LOGARITHMIC_BIN = 1;



template <typename T, typename T_INT>
class ArrayProperties
{
    T *ptr;
    T_INT size;
    
    //! variables to keep track of min and max values
    T minimumValue, maximumValue;       // min and max values of array entries
    T absoluteMaximumValue, absoluteMinimumValue;       // min and max absolute values of array entries
    T_INT minimumIndex, maximumIndex;   // the indices of the min and max values
    T_INT absoluteMinimumIndex, absoluteMaximumIndex;   // the indices of the min and max absolute values
    bool minMaxComputation;             // true if the min and max values were computed - computed only on request
    
    //! variables to keep track of the sum of entries
    double sumValue;                    // the sum of all array entries
    double sum2Value;                   // the sum of the square of the array entries
    bool sumComputation;                // true if sum computation was performed
    
    //! variables to keep track of the median
    T medianValue;                      // median value of the array values
    T medianApprox;                     // approximative median value
    T medianApproxError;                // maximal deviation for the approximative median value
    bool medianComputed;                // true if median value was computed
    bool medianApproxComputed;          // true if approximative median value was computed
    
    //! variables to keep track of the histogram
    T_INT histogramElements;            // number of histogram bins
    std::vector<T_INT> valuesHistogram; // the histogram values corresponding to each bin
    std::vector<double> valuesHistogramWeighted; // the values of the weighted histogram
    std::vector<T> binsHistogram;       // the edges of each histogram bin (there are n+1 values)
    bool equalBins;                     // true if the bin size are equal
    bool logarithmicBins;               // true if the bins are logarithmically spaced
    bool histogramComputed;             // true if the histogram was computed
    bool histogramWeightedComputed;     // true if a weighted histogram was computed
    
    //! additional options
    bool sorted;                        // true if array was sorted
    bool uniqueOnly;                    // true if all array elements are unique
    
    
    //! private functions
    void defaultClassVariables();       // initializes all class variables to default values
    void computeMinMax();               // computes the min and max values of entries
    void computeSum();                  // computes sum of the elements
    void computeMedian(int const binsType);  // computes the true median of the data set
    void computeMedianApprox(size_t const noBins,
                             int const binsType); // computes an approximation for the median
    void setHistogramScale(T const lowerLimit, T const upperLimit); //sets the bins of the histogram
    void computeHistogram();            // computes the histogram
    void computeHistogramWeighted(T *weight);  // computes a weighted histogram
    T getHistogramBinBoundary(size_t const i); // returns the lower histogram bin boundary corresponding to bin i
    void sortArray();                   // sorts the array entries
    void uniqueElements();              // keeps only unique elements
    
    public:
        ArrayProperties();
        ArrayProperties(T *pointerData, T_INT const noElements);
        
        
        void newData(T *pointerData, T_INT const noElements);   // insert new array as argument for the class
        
        //! functions for min and max values
        T valueMaximum();
        T valueMinimum();
        T_INT indexMaximum();
        T_INT indexMinimum();
        
        //! functions for sum values
        double sum();                    // computes the sum of all array entries
        double squareSum();              // computes the sum of all array entries square
        std::vector<double> cumulativeSum(int const sumDirection=ORDER_ASCENDING); //return the cumulative sum of the entries
        void cumulativeSum(double *returnValues, T_INT const noElements, int const sumDirection=ORDER_ASCENDING); //return the cumulative sum of the entries
        
        //! functions for average and median
        T mean();   // returns the arithmetic mean
        T median(int const binsType = LINEAR_BIN);  // returns the mean of the data set - it uses bining to find the median rapidly
        T medianApproximative(size_t const noBins = 100,
                              int const binsType = LINEAR_BIN); // returns an approximative mean, the more bins, the more accurate; It uses bining to find an approximation for the median - the higher the bin size, the better the approximation.
        T medianApproximativeError(); // returns the maximum error in determining the true median
        
        //! functions for histogram computation
        void histogram(T lowerLimit, T upperLimit, T_INT noElements, int const binsType = LINEAR_BIN); // computes the histogram using equal spaced bins (or logarithmically spaced)
        void histogramWeighted(T lowerLimit, T upperLimit, T_INT noElements, T *weight, int const binsType = LINEAR_BIN); // computes a weighted histogram (instead of increasing the count by 1, the count increases by weight[i])
        void histogram(T const *pointerBins, T_INT noBins, int const binsType = LINEAR_BIN);           // computes the histogram using the bin limits specified in the pointer. NOTE: the pointer 'pointerBins' must have length noBins+1.
        T_INT noHistogramBins();                                // returns number of histogram bins
        std::vector<T_INT> histogramValues();                   // returns a vector with histogram values in each bin
        std::vector<double> histogramWeightedValues();          // returns a vector with the values of the weighted histogram
        void histogramValues(T_INT *returnValues, T_INT const size); // returns the histogram values in the pointer supplied
        std::vector<T> histogramBins();                         // returns a vector with the central position of the bins
        void histogramBins(T *returnValues, T_INT const size);  // returns the central position of the bins in the pointer supplied
        
        //! additional functions
        void sort();                // sorts the array entries in ascending order
        T_INT unique();             // keeps only unique array entries and returns the number of unique entries
};



//! Functions constructors
template <typename T, typename T_INT> ArrayProperties<T,T_INT>::ArrayProperties()
{
    this->defaultClassVariables();
    ptr = 0;
    size = 0;
}
template <typename T, typename T_INT> ArrayProperties<T,T_INT>::ArrayProperties(T *pointerData, T_INT const noElements)
{
    this->defaultClassVariables();
    ptr = pointerData;
    size = noElements;
}



//! Public functions
/* Change the array whose properties are computed in the class. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::newData(T *pointerData, T_INT const noElements)
{
    this->defaultClassVariables();
    ptr = pointerData;
    size = noElements;
}


//! Functions related to the min and max entries of the array
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::valueMaximum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not minMaxComputation )
        this->computeMinMax();
    return maximumValue;
}
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::valueMinimum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not minMaxComputation )
        this->computeMinMax();
    return minimumValue;
}
template <typename T, typename T_INT> T_INT ArrayProperties<T,T_INT>::indexMaximum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not minMaxComputation )
        this->computeMinMax();
    return maximumIndex;
}
template <typename T, typename T_INT> T_INT ArrayProperties<T,T_INT>::indexMinimum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not minMaxComputation )
        this->computeMinMax();
    return minimumIndex;
}


//! Functions related to the element sum of the matrix
template <typename T, typename T_INT> double ArrayProperties<T,T_INT>::sum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not sumComputation )
        this->computeSum();
    return sumValue;
}
template <typename T, typename T_INT> double ArrayProperties<T,T_INT>::squareSum()
{
    BOOST_ASSERT( ptr!=0 );
    if ( not sumComputation )
        this->computeSum();
    return sum2Value;
}
/* Returns the cumulative sum of the entries. */
template <typename T, typename T_INT> std::vector<double> ArrayProperties<T,T_INT>::cumulativeSum(int const sumDirection)
{
    std::vector<double> temp;
    temp.reserve( size );
    if ( sumDirection==ORDER_ASCENDING )
    {
        temp.push_back( ptr[0] );
        for (T_INT i=1; i<size; ++i)
            temp.push_back( ptr[i]+temp.back() );
    }
    else if ( sumDirection==ORDER_DESCENDING )
    {
        temp.assign( size, T(0) );
        temp.at(size-1) = ptr[size-1];
        for (T_INT i=size-1; i--; )
        {
            temp[i] = ptr[i] + temp[i+1];
        }
    }
    else throwError( "Unknow value for the 'sumDirection' argument of function 'ArrayProperties::cumulativeSum'. This parameter gives the direction of the cumulative sum and can take the values 'ORDER_ASCENDING' or 'ORDER_DESCENDING'." );
    return temp;
}
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::cumulativeSum(double *returnValues,
                                                                                   T_INT const noElements,
                                                                                   int const sumDirection)
{
    if( size!=noElements ) throwError( "The value of the second argument of the function 'ArrayProperties::cumulativeSum' needs to be the same as the size of the array whose cumulative sum is computed." );
    if ( sumDirection==ORDER_ASCENDING )
    {
        returnValues[0] = ptr[0];
        for (T_INT i=1; i<size; ++i)
            returnValues[i] = ptr[i] + returnValues[i-1];
    }
    else if ( sumDirection==ORDER_DESCENDING )
    {
        returnValues[size-1] = ptr[size-1];
        for (T_INT i=size-1;  i--; )
            returnValues[i] = ptr[i] + returnValues[i+1];
    }
    else throwError( "Unknow value for the 'sumDirection' argument of function 'ArrayProperties::cumulativeSum'. This parameter gives the direction of the cumulative sum and can take the values 'ORDER_ASCENDING' or 'ORDER_DESCENDING'." );
}



//! Functions related to computing the mean and median of the values
/* Returns the arithmetic mean of the values. */
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::mean()
{
    return sum() / size;
}
/* Returns the median of a set of values. */
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::median(int const binsType)
{
    if ( not medianComputed )
        this->computeMedian( binsType );
    return medianValue;
}
/* Returns and approximative value for the median. */
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::medianApproximative(size_t const noBins,
                                                                                      int const binsType)
{
    if ( not medianApproxComputed )
        computeMedianApprox( noBins, binsType );
    return medianApprox;
}
/* Returns the maximal error for the approximative value of the median. */
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::medianApproximativeError()
{
    if ( not medianApproxComputed )
        throwError ( "Must call the function 'medianApproximative()' before asking for the errors in the determination of the approximative meadian." );
    return medianApproxError;
}




//! Functions related to histogram computation
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::histogram(T lowerLimit,
                                                                               T upperLimit,
                                                                               T_INT noElements,
                                                                               int const binsType)
{
    BOOST_ASSERT( ptr!=0 );
    equalBins = true;
    if ( binsType==LOGARITHMIC_BIN )
        logarithmicBins = true;
    else if ( binsType!=LINEAR_BIN ) throwError( "Unknow value for argument 'binsType' in function 'ArrayProperties::histogram'." );
    histogramElements = noElements;
    this->setHistogramScale( lowerLimit, upperLimit );
    this->computeHistogram();
    histogramComputed = true;
    histogramWeightedComputed = false;
}
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::histogramWeighted(T lowerLimit,
                                                                                       T upperLimit,
                                                                                       T_INT noElements,
                                                                                       T *weight,
                                                                                       int const binsType)
{
    BOOST_ASSERT( ptr!=0 );
    equalBins = true;
    if ( binsType==LOGARITHMIC_BIN )
        logarithmicBins = true;
    else if ( binsType!=LINEAR_BIN ) throwError( "Unknow value for argument 'binsType' in function 'ArrayProperties::histogramWeighted'." );
    histogramElements = noElements;
    this->setHistogramScale( lowerLimit, upperLimit );
    this->computeHistogramWeighted( weight );
    histogramWeightedComputed = true;
    histogramComputed = false;
}
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::histogram(T const *pointerBins,
                                                                               T_INT noElements,
                                                                               int const binsType)
{
    BOOST_ASSERT( ptr!=0 );
    equalBins = false;
    if ( binsType==LOGARITHMIC_BIN )
        logarithmicBins = true;
    else if ( binsType==LINEAR_BIN ) throwError( "Unknow value for argument 'binsType' in function 'ArrayProperties::histogram'." );
    histogramElements = noElements;
    binsHistogram.clear();  //save the bin limits
    binsHistogram.reserve(histogramElements+1);
    for (size_t i=0; i<=histogramElements; ++i)
        binsHistogram.push_back( pointerBins[i] );
    this->computeHistogram();
    histogramComputed = true;
    histogramWeightedComputed = false;
}
template <typename T, typename T_INT> T_INT ArrayProperties<T,T_INT>::noHistogramBins()
{
    if ( not histogramComputed and not histogramWeightedComputed ) throwError( "You must call the function 'ArrayProperties::histogram' (which computes the actual histogram) before calling the function 'ArrayProperties::noHistogramBins' (which returns the number of histogram bins)." );
    return histogramElements;
}
template <typename T, typename T_INT> std::vector<T_INT> ArrayProperties<T,T_INT>::histogramValues()
{
    if ( not histogramComputed) throwError( "You must call the function 'ArrayProperties::histogram' (which computes the actual histogram) before calling the function 'ArrayProperties::histogramValues' (which returns the histogram values)." );
    return valuesHistogram;
}
template <typename T, typename T_INT> std::vector<double> ArrayProperties<T,T_INT>::histogramWeightedValues()
{
    if ( not histogramWeightedComputed ) throwError( "You must call the function 'ArrayProperties::histogramWeighted' (which computes the actual weighted histogram) before calling the function 'ArrayProperties::histogramWeightedValues' (which returns the weighted histogram values)." );
    return valuesHistogramWeighted;
}
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::histogramValues(T_INT *returnValues, T_INT const noElements)
{
    if ( not histogramComputed) throwError( "You must call the function 'ArrayProperties::histogram' (which computes the actual histogram) before calling the function 'ArrayProperties::histogramValues' (which returns the histogram values)." );
    if( valuesHistogram.size()!=noElements ) throwError( "The value of the second argument of the function 'ArrayProperties::histogramValues' needs to be the same as the number of bins of the histogram." );
    for (size_t i=0; i<valuesHistogram.size(); ++i)
        returnValues[i] = valuesHistogram[i];
}
template <typename T, typename T_INT> std::vector<T> ArrayProperties<T,T_INT>::histogramBins()
{
    if ( not histogramComputed and not histogramWeightedComputed ) throwError( "You must call the function 'ArrayProperties::histogram' (which computes the actual histogram) before calling the function 'ArrayProperties::histogramBins' (which returns the center of each histogram bins)." );
    std::vector<T> temp;
    temp.reserve(histogramElements);
    if ( not logarithmicBins )
        for (size_t i=0; i<histogramElements; ++i)
            temp.push_back( (binsHistogram[i]+binsHistogram[i+1])/2. );
    else
        for (size_t i=0; i<histogramElements; ++i)
            temp.push_back( T(std::sqrt( binsHistogram[i]*binsHistogram[i+1] )) );
    return temp;
}
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::histogramBins(T *returnValues, T_INT const noElements)
{
    if ( not histogramComputed and not histogramWeightedComputed ) throwError( "You must call the function 'ArrayProperties::histogram' (which computes the actual histogram) before calling the function 'ArrayProperties::histogramBins' (which returns the center of each histogram bins)." );
    if( valuesHistogram.size()!=noElements ) throwError( "The value of the second argument of the function 'ArrayProperties::histogramBins' needs to be the same as the number of bins of the histogram." );
    if ( not logarithmicBins )
        for (size_t i=0; i<noElements; ++i)
            returnValues[i] = (binsHistogram[i]+binsHistogram[i+1]) / 2. ;
    else
        for (size_t i=0; i<noElements; ++i)
            returnValues[i] = std::sqrt( binsHistogram[i]*binsHistogram[i+1]);
}




//! Private functions
/* Computes the min and max value of the array. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeMinMax()
{
#define ABSOLUTE(x) (((x)<0)?(-(x)):(x))
    minimumValue = maximumValue = ptr[0];
    absoluteMinimumValue = ABSOLUTE( ptr[0] );
    minimumIndex = maximumIndex = absoluteMinimumIndex = T_INT(0);
    for (T_INT i=1; i<size; ++i)
    {
        if ( minimumValue>ptr[i] )
        {
            minimumValue = ptr[i];
            minimumIndex = i;
        }
        if ( maximumValue<ptr[i] )
        {
            maximumValue = ptr[i];
            maximumIndex = i;
        }
        if ( absoluteMinimumValue>ABSOLUTE(ptr[i]) )
        {
            absoluteMinimumValue = ABSOLUTE( ptr[i] );
            absoluteMinimumIndex = i;
        }
    }
    
    //compute the absolute maximum
    if ( ABSOLUTE(maximumValue)>ABSOLUTE(minimumValue) )
    {
        absoluteMaximumValue = ABSOLUTE( maximumValue );
        absoluteMaximumIndex = maximumIndex;
    }
    else
    {
        absoluteMaximumValue = ABSOLUTE( minimumValue );
        absoluteMaximumIndex = minimumIndex;
    }
    
    
    minMaxComputation = true;
}

/* Computes the sum of all the array entries. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeSum()
{
    sumValue = T(0);
    sum2Value = T(0);
    for (T_INT i=0; i<size; ++i)
    {
        sumValue += ptr[i];
        sum2Value += ptr[i]*ptr[i];
    }
    sumComputation = true;
}

/* Computes the median of a set of values. It uses a bining technique to do so. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeMedian(int const binsType)
{
    size_t const noBins = 1000;
    T_INT position = size/2, remainingSize = size;
    T minVal = this->valueMinimum();
    T maxVal = this->valueMaximum();
    if ( minVal==maxVal or size<=1 )
    {
        medianValue = minVal;
        medianComputed = true;
        return;
    }
    
    if ( binsType==LOGARITHMIC_BIN and minVal==T(0.) )
        minVal = 1.e-8 * maxVal;
    else if ( binsType==LOGARITHMIC_BIN and minVal*maxVal<T(0.) )
        throwError ("Cannot compute the median using logarithmic bins since the data has both positive and negative values.");
    
    
    while ( remainingSize>10*noBins )
    {
        // compute histogram between minVal to maxVal
        ArrayProperties<T,T_INT> prop( ptr, size );
        prop.histogram( minVal, maxVal, noBins, binsType );
        vector<T_INT> histValues = prop.histogramValues();
        
        // find histogram bin with the median value
        T_INT res = 0, i=0;
        while ( res<position )
            res += histValues[i++];
        --i;
        
        // find new minimum and maximum values where the median is
        minVal = prop.getHistogramBinBoundary( i );
        maxVal = prop.getHistogramBinBoundary( i+1 );
        
        // update the position of the median in the new histogram interval
        position -= (res - histValues[i]);
        remainingSize = histValues[i];
    }
    
    
    // Now when the are only very little elements -> do sort and find the median
    vector<T> temp;    //temporary variable to store the value for sorting
    temp.reserve( remainingSize );
    for (T_INT i=0; i<size; ++i)
        if ( ptr[i]>=minVal and ptr[i]<maxVal )
            temp.push_back( ptr[i] );
    
    sort( temp.begin(), temp.end() );
    
    medianValue = temp[position];
    medianComputed = true;
}
/* Computes the approximative value for the median of a set of values. It uses a bining technique to do so. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeMedianApprox(size_t const noBins,
                                                                                         int const binsType)
{
    ArrayProperties<T,T_INT> temp( ptr, size );
    T minValue = this->valueMinimum();
    T maxValue = this->valueMaximum();
    if ( minValue==maxValue or size<=1 )
    {
        medianApprox = minValue;
        medianApproxError = 0.;
        medianApproxComputed = true;
        return;
    }
    
    if ( binsType==LOGARITHMIC_BIN and minValue==T(0.) )
        minValue = 1.e-6 * maxValue;
    else if ( binsType==LOGARITHMIC_BIN and minValue*maxValue<=T(0.) )
        throwError ("Cannot compute the median approximation using logarithmic bins since the data has both positive and negative values.");
    temp.histogram( minValue, maxValue, noBins, binsType );
    std::vector<T_INT> tempValues = temp.histogramValues();
    T_INT const half = size / 2;
    T_INT res = 0;
    int i=0; // keep track of histogram bin where the median value is
    while ( res<half )
        res += tempValues[i++];
    --i;
    
    T minVal = temp.getHistogramBinBoundary( i );
    T maxVal = temp.getHistogramBinBoundary( i+1 );
    
    medianApprox = (maxVal + minVal) / 2.;
    medianApproxError = ( maxVal - minVal ) / 2.;
    medianApproxComputed = true;
}



/* Creates the bin values for the histogram computation. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::setHistogramScale(T const lowerLimit,
                                                                                       T const upperLimit)
{
    binsHistogram.clear();
    binsHistogram.reserve(histogramElements+1);
    if ( not logarithmicBins )
    {
        T dx = (upperLimit-lowerLimit) / histogramElements;
        T temp = lowerLimit;
        for (size_t i=0; i<=histogramElements; ++i)
        {
            binsHistogram.push_back( temp );
            temp += dx;
        }
    }
    else
    {
        T dx = std::log10(upperLimit/lowerLimit) / histogramElements;
        dx = pow( double(10.), double(dx) );
        T temp = lowerLimit;
        for (size_t i=0; i<=histogramElements; ++i)
        {
            binsHistogram.push_back( temp );
            temp *= dx;
        }
    }
}
/* This function computes the histogram for the given array. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeHistogram()
{
    valuesHistogram.assign( histogramElements, T_INT(0) );
    if ( equalBins and not logarithmicBins )
    {
        T dx = (binsHistogram.at(histogramElements)-binsHistogram.at(0)) / histogramElements;
        T minX = binsHistogram.at(0);
        for (T_INT i=0; i<size; ++i)
        {
            if ( ptr[i]<binsHistogram[0] or ptr[i]>binsHistogram[histogramElements] )
                continue;
            T_INT temp = T_INT( (ptr[i]-minX) / dx );
            if ( temp<histogramElements ) ++valuesHistogram[temp];
        }
    }
    else if ( equalBins )
    {
        T dx = std::log10(binsHistogram.at(histogramElements)/binsHistogram.at(0)) / histogramElements;
        T minX = binsHistogram.at(0);
        for (T_INT i=0; i<size; ++i)
        {
            if ( ptr[i]<binsHistogram[0] or ptr[i]>binsHistogram[histogramElements] )
                continue;
            T_INT temp = T_INT( std::log10(ptr[i]/minX) / dx );
            if ( temp<histogramElements ) ++valuesHistogram[temp];
        }
    }
    else
    {
        for (T_INT i=0; i<size; ++i)    //binary search the position of the value in the histogram table for random bin sizes
        {
            if ( ptr[i]<binsHistogram[0] or ptr[i]>binsHistogram[histogramElements] )
                continue;
            T_INT dx = histogramElements/4;
            T_INT temp = histogramElements/2;
            while( dx>1 )
            {
                if ( ptr[i]>binsHistogram[temp] ) temp += dx;
                else temp -= dx;
                dx /= 2;
            }
            // now 'temp' is very close to the bin number where lies 'ptr[i]' - so do a small linear search
            if ( ptr[i]>=binsHistogram[temp] )
                for (T_INT j=temp; j<histogramElements; ++j)
                    if ( ptr[i]>=binsHistogram[j] and ptr[i]<binsHistogram[j+1] )
                    {
                         ++valuesHistogram[j];
                         break;
                    }
            else
                for (T_INT j=temp; j>0; --j)
                    if ( ptr[i]>=binsHistogram[j-1] and ptr[i]<binsHistogram[j] )
                    {
                         ++valuesHistogram[j];
                         break;
                    }
        }
    }
}
/* This function computes the histogram for the given array. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::computeHistogramWeighted(T *weight)
{
    valuesHistogramWeighted.assign( histogramElements, T(0.) );
    if ( equalBins and not logarithmicBins )
    {
        T dx = (binsHistogram.at(histogramElements)-binsHistogram.at(0)) / histogramElements;
        T minX = binsHistogram.at(0);
        for (T_INT i=0; i<size; ++i)
        {
            if ( ptr[i]<binsHistogram[0] or ptr[i]>binsHistogram[histogramElements] )
                continue;
            T_INT temp = T_INT( (ptr[i]-minX) / dx );
            if ( temp<histogramElements )
                valuesHistogramWeighted[temp] += weight[i];
        }
    }
    else if ( equalBins )
    {
        T dx = std::log10(binsHistogram.at(histogramElements)/binsHistogram.at(0)) / histogramElements;
        T minX = binsHistogram.at(0);
        for (T_INT i=0; i<size; ++i)
        {
            if ( ptr[i]<binsHistogram[0] or ptr[i]>binsHistogram[histogramElements] )
                continue;
            T_INT temp = T_INT( std::log10(ptr[i]/minX) / dx );
            if ( temp<histogramElements )
                valuesHistogramWeighted[temp] += weight[i];
        }
    }
}

/* Returns the lower histogram bin boundary corresponding to bin i. */
template <typename T, typename T_INT> T ArrayProperties<T,T_INT>::getHistogramBinBoundary(size_t const i)
{
    if ( histogramElements==0 )
        throwError( "Trying to access the bin boundaries for the histogram computation but not histogram computation was initialized in class 'ArrayProperties'." );
    if ( i<0 or i>binsHistogram.size()-1 )
    {
        cout << i << "\t" << binsHistogram.size() << "\t" << size << "\n";
        throwError( "Invalid argument value in function 'ArrayProperties::getHistogramBinBoundary'." );
    }
    return binsHistogram[i];
}



/* Initializes all class variables to default values. */
template <typename T, typename T_INT> void ArrayProperties<T,T_INT>::defaultClassVariables()
{
    minimumValue = maximumValue = T(0);
    minimumIndex = maximumIndex = T_INT(-1);
    minMaxComputation = false;
    
    sumValue = sum2Value = T(0);
    sumComputation = false;
    
    medianComputed = false;
    medianApproxComputed = false;
    
    histogramElements = 0;
    equalBins = false;
    logarithmicBins = false;
    histogramComputed = false;
    histogramWeightedComputed = false;
    
    sorted = false;
    uniqueOnly = false;
}

#endif
