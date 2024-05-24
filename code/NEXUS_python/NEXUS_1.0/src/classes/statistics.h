

#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
typedef boost::mt19937 base_generator_type;


template <typename T> 
class Statistics
{
    public:
    std::vector<T> data;// pointer to the data location
    size_t size;        // number of elements in the data array
    
    T mean;             // the value of the mean
    T meanErrorLower;   // the error in the mean determination (16% of mean values using bootstrap)
    T meanErrorUpper;   // the error in the mean determination (84% of mean values using bootstrap)
    
    T median;           // the value of the median
    T medianErrorLower; // the error in the median (16% of median values using bootstrap)
    T medianErrorUpper; // the error in the median (84% of median values using bootstrap)
    
    T variance;         // the variance of the data
    T varianceErrorLower; // the error in the variance (16% of variance values using bootstrap)
    T varianceErrorUpper; // the error in the variance (84% of variance values using bootstrap)
    
    T lowerPercentile;  // value of 16% of the data percentile
    T upperPercentile;  // value of 84% of the data percentile
    
    int noBootstrapSamples;     // can give the number of of bootstrap samples to be used
    int maxSizeForBootstrap;    // if there are more data points than this number than the program does not compute the error using bootstrap, but only divides the std by sqrt{N}
    T   lowerPercentileValue;   // can give a different value for the lower percentile than 16%
    T   upperPercentileValue;   // can give a different value for the upper percentile than 84%
    
    bool computeMean;           // true if to compute the average
    bool computeMedian;         // true if to compute the median
    bool computeVariance;       // true if to compute the variance
    bool computePercentile;     // true if to compute the distribution percentiles
    
    
    //! functions for the user
    Statistics();               // constructor
    void newData(T *newDataPtr, size_t dataSize);   // add a new set of data
    void computeStatistics();   // computes the statistics
    
    
    //! private functions
    private:
        void bootstrapSamples(size_t seed, std::vector<T> &outputData);  // returns a bootstraped sample of the data
        T meanComputation(std::vector<T> &inData, T *Variance=NULL);     // computes the mean and variance
        T medianComputation(std::vector<T> &inData, bool *sorted);       // computes the median
        void percentileComputation(std::vector<T> &inData, bool *sorted, T *lowerP, T * upperP); // computes the percentiles for the given input distribution
};


// Constructors
template <typename T> Statistics<T>::Statistics()
{
    size = 0;
    noBootstrapSamples = 100;
    maxSizeForBootstrap = 1000000;  //1e6
    lowerPercentileValue = T(16.);
    upperPercentileValue = T(84.);
    computeMean = true;
    computeMedian = true;
    computeVariance = true;
    computePercentile = true;
}



// add new data
template <typename T> void Statistics<T>::newData(T *newDataPtr, size_t dataSize)
{
    size = dataSize;
    data.assign( size, T(0.) );
    for (size_t i=0; i<dataSize; ++i)
        data[i] = newDataPtr[i];
}


// compute the statistics
template <typename T> void Statistics<T>::computeStatistics()
{
    // get the values for the main data
    bool sorted = false;
    if ( computeMean or computeVariance )
        mean = meanComputation( data, &variance );
    if ( computeMedian )
        median = medianComputation( data, &sorted );
    if ( computePercentile )
        percentileComputation( data, &sorted, &lowerPercentile, &upperPercentile );
    
    // estimate errors for the mean and median
    if ( noBootstrapSamples<=0 )    //return if no bootstrap is required
        return;
    if ( maxSizeForBootstrap<=size )    // return if the number of sample is large enough to get an error estimate without the need for bootstraping
    {
        if ( not computePercentile )
            percentileComputation( data, &sorted, &lowerPercentile, &upperPercentile );
        meanErrorLower = mean + (lowerPercentile-mean) / std::sqrt( size );
        meanErrorUpper = mean + (upperPercentile-mean) / std::sqrt( size );
        varianceErrorLower = variance + meanErrorLower - mean;
        varianceErrorUpper = variance + meanErrorUpper - mean;
        medianErrorLower = median + meanErrorLower - mean;
        medianErrorUpper = median + meanErrorUpper - mean;
        return;
    }
    
    // do the actual bootstrap computations
    std::vector<T> bootstrapData( size, T(0.) );  //stores a sample of bootstraped data
    std::vector<T> meanSet( noBootstrapSamples, T(0.) );  //stores the mean for each sample of bootstraped data
    std::vector<T> varianceSet( noBootstrapSamples, T(0.) );//stores the variance for each sample of bootstraped data
    std::vector<T> medianSet( noBootstrapSamples, T(0.) );  //stores the median for each sample of bootstraped data
    size_t seed = 10000;
    for (int i=0; i<noBootstrapSamples; ++i)
    {
        seed += 1000*i + 200*i + 30*i + 4*i*i + 4321;
        bootstrapSamples( seed, bootstrapData );
        
        sorted = false;
        if ( computeMean ) meanSet[i] = meanComputation( bootstrapData );
        if ( computeVariance ) meanComputation( bootstrapData, &(varianceSet[i]) );
        if ( computeMedian ) medianSet[i] = medianComputation( bootstrapData, &sorted );
    }
    
    // get the errors in the mean and median
    sorted = false;
    if ( computeMean )
        percentileComputation( meanSet, &sorted, &meanErrorLower, &meanErrorUpper );
    sorted = false;
    if ( computeVariance )
        percentileComputation( varianceSet, &sorted, &varianceErrorLower, &varianceErrorUpper );
    sorted = false;
    if ( computeMedian )
        percentileComputation( medianSet, &sorted, &medianErrorLower, &medianErrorUpper );
}


//! Private functions
// returns a bootstraped sample of the data
template <typename T> void Statistics<T>::bootstrapSamples(size_t seed, std::vector<T> &outputData)
{
    base_generator_type gen( seed );
    boost::random::uniform_int_distribution<> uni(0, size-1);
    
    for (size_t i=0; i<size; ++i)
        outputData[i] = data[ uni(gen) ];
}

// computes the mean of a data set
template <typename T> T Statistics<T>::meanComputation(std::vector<T> &inData, T *Variance)
{
    double sum = 0.;
    size_t noElem = inData.size();
    for (size_t i=0; i<noElem; ++i)
        sum += double( inData[i] );
    T meanValue = sum / noElem;
    
    // compute the variance
    if ( computeVariance and Variance )
    {
        sum = 0;
        for (size_t i=0; i<noElem; ++i)
            sum += double( (inData[i]-meanValue) * (inData[i]-meanValue) );
        *Variance = sum / noElem;
    }
    
    return meanValue;
}


// computes the median of a data set
template <typename T> T Statistics<T>::medianComputation(std::vector<T> &inData, bool *sorted)
{
    if ( not *sorted )
    {
        sort( inData.begin(), inData.end() );
        *sorted = true;
    }
    return inData[ inData.size()/2 ];
}

// computes the percentiles of a data set
template <typename T> void Statistics<T>::percentileComputation(std::vector<T> &inData, bool *sorted, T *lowerP, T * upperP)
{
    if ( not *sorted )
    {
        sort( inData.begin(), inData.end() );
        *sorted = true;
    }
    size_t noElem = inData.size();
    size_t i1 = size_t( lowerPercentileValue * noElem / 100 );
    size_t i2 = size_t( upperPercentileValue * noElem / 100 );
    *lowerP = inData[i1];
    *upperP = inData[i2];
}







template <typename T, typename T2> 
class Histogram
{
    public:
    std::vector<T> data;// pointer to the data location
    size_t size;        // number of elements in the data array
    T minRange;         // minimum range value
    T maxRange;         // maximum range value
    size_t noBins;      // number of bins
    
    std::vector<T2> histogram;          // the histogram
    std::vector<T2> histogramLower;     // lower boostrap allowed histogram values
    std::vector<T2> histogramUpper;     // upper boostrap allowed histogram values
    
    std::vector<T2> histogramCum;       // the cumulative histogram
    std::vector<T2> histogramCumLower;
    std::vector<T2> histogramCumUpper;
    
    int noBootstrapSamples;     // can give the number of of bootstrap samples to be used
    float lowerPercentileValue;   // can give a different value for the lower percentile than 16%
    float upperPercentileValue;   // can give a different value for the upper percentile than 84%
    
    Histogram();                // the constructor
    void newData(T *newDataPtr, size_t dataSize, int const numberBins, T const minValue, T const maxValue);   // add a new set of data
    void computeHistogram();    // computes the histogram
    void computeCumulativeHistogram();    // computes the cumulative histogram
    
    private:
        void bootstrapSamples(size_t seed, std::vector<T> &outputData);        // returns a bootstraped sample of the data
        void histogramComputation(std::vector<T> &input, T2 *output);          // computes the histogram of values using the same seting as in the main class
        void histogramCumulativeComputation(std::vector<T> &input, T2 *output);// computes the cumulative histogram of values using the same seting as in the main class
        void percentileComputation(std::vector<T2> &input, bool *sorted, T2 *lowerP, T2 * upperP); // computes the percentiles for the given input
};



// Constructors
template <typename T, typename T2> Histogram<T,T2>::Histogram()
{
    size = 0;
    noBootstrapSamples = 100;
    lowerPercentileValue = float(16.);
    upperPercentileValue = float(84.);
}

// add new data
template <typename T, typename T2> void Histogram<T,T2>::newData(T *newDataPtr, size_t dataSize, int const numberBins, T const minValue, T const maxValue)
{
    size = dataSize;
    data.assign( size, T(0.) );
    for (size_t i=0; i<dataSize; ++i)
        data[i] = newDataPtr[i];
    noBins = numberBins;
    minRange = minValue;
    maxRange = maxValue;
}

// computes the histogram
template <typename T, typename T2> void Histogram<T,T2>::computeHistogram()
{
    // get the histogram of the data
    histogram.assign( noBins, T2(0) );
    histogramLower.assign( noBins, T2(0) );
    histogramUpper.assign( noBins, T2(0) );
    this->histogramComputation( data, &(histogram[0]) );
    
    // estimate variability in the histogram
    if ( noBootstrapSamples<=0 ) return;
    std::vector<T>  bootstrapData( size, T(0) );  //stores a sample of bootstraped data
    std::vector<T2> tempHistogram( noBins*noBootstrapSamples, T2(0) );   //stores the histogram for each bootstrapped iteration
    size_t seed = 10000;
    for (int i=0; i<noBootstrapSamples; ++i)
    {
        seed += 1000*i + 200*i + 30*i + 4*i*i + 4321;
        this->bootstrapSamples( seed, bootstrapData );
        
        this->histogramComputation( bootstrapData, &(tempHistogram[i*noBins]) );
    }
    
    
    // get the lower and upper percentiles for the histogram estimation
    std::vector<T2> temp( noBootstrapSamples, T2(0) );
    for (int i=0; i<noBins; ++i)
    {
        for (int j=0; j<noBootstrapSamples; ++j)
            temp[j] = tempHistogram[j*noBins+i];
        bool sorted = false;
        percentileComputation( temp, &sorted, &(histogramLower[i]), &(histogramUpper[i]) );
    }
}

// computes the cumulative histogram
template <typename T, typename T2> void Histogram<T,T2>::computeCumulativeHistogram()
{
    // get the histogram of the data
    histogramCum.assign( noBins, T2(0) );
    histogramCumLower.assign( noBins, T2(0) );
    histogramCumUpper.assign( noBins, T2(0) );
    this->histogramCumulativeComputation( data, &(histogramCum[0]) );
    
    // estimate variability in the histogram
    if ( noBootstrapSamples<=0 ) return;
    std::vector<T>  bootstrapData( size, T(0.) );  //stores a sample of bootstraped data
    std::vector<T2> tempHistogram( noBins*noBootstrapSamples, T2(0.) );   //stores the histogram for each bootstrapped iteration
    size_t seed = 10000;
    for (int i=0; i<noBootstrapSamples; ++i)
    {
        seed += 1000*i + 200*i + 30*i + 4*i*i + 4321;
        this->bootstrapSamples( seed, bootstrapData );
        
        this->histogramCumulativeComputation( bootstrapData, &(tempHistogram[i*noBins]) );
    }
    
    
    // get the lower and upper percentiles for the histogram estimation
    std::vector<T2> temp( noBootstrapSamples, T2(0.) );
    for (size_t i=0; i<noBins; ++i)
    {
        for (int j=0; j<noBootstrapSamples; ++j)
            temp[j] = tempHistogram[j*noBins+i];
        bool sorted = false;
        percentileComputation( temp, &sorted, &(histogramCumLower[i]), &(histogramCumUpper[i]) );
    }
}



//! Private functions
// returns a bootstraped sample of the data
template <typename T, typename T2> void Histogram<T,T2>::bootstrapSamples(size_t seed, std::vector<T> &outputData)
{
    base_generator_type gen( seed );
    boost::random::uniform_int_distribution<> uni(0, size-1);
    
    for (size_t i=0; i<size; ++i)
        outputData[i] = data[ uni(gen) ];
}

// computes the percentiles of a data set
template <typename T, typename T2> void Histogram<T,T2>::percentileComputation(std::vector<T2> &input, bool *sorted, T2 *lowerP, T2 * upperP)
{
    if ( not *sorted )
    {
        sort( input.begin(), input.end() );
        *sorted = true;
    }
    size_t noElem = input.size();
    size_t i1 = size_t( lowerPercentileValue * noElem / 100 );
    size_t i2 = size_t( upperPercentileValue * noElem / 100 );
    *lowerP = input[i1];
    *upperP = input[i2];
}

// computes the histogram
template <typename T, typename T2> void Histogram<T,T2>::histogramComputation(std::vector<T> &input, T2 *output)
{
    T step = 1;
    for (size_t i=0; i<size; ++i)
    {
        int bin = (int) ( (input[i]-minRange)/step );
        if ( bin<0 or bin>=noBins ) continue;
        output[bin] += T2(1);
    }
}

// computes the cumulative histogram
template <typename T, typename T2> void Histogram<T,T2>::histogramCumulativeComputation(std::vector<T> &input, T2 *output)
{
    this->histogramComputation( input, output );
    for (size_t i=1; i<noBins; ++i)
        output[i] += output[i-1];
}

