#ifndef RANGE_HEADER
#define RANGE_HEADER

#include <cmath>

#include "throwError.h"


/* Creates an array of logarithmically spaced values between minValue to maxValue and with noValues values. */
template <typename T>
void linearRange(T minValue, T maxValue, size_t noValues,
                 T * result)
{
    
    if ( noValues<=1 ) throwError( "The number of bins in function 'linearRange' (which is 'noValues') must be larger than 1. " );
    double const step = (maxValue-minValue) / (noValues-1);
    
    result[0] = minValue;
    for (size_t i=1; i<noValues; ++i)
        result[i] = result[i-1] + step;
}


/* Creates an array of logarithmically spaced values between minValue to maxValue and with noValues values. */
template <typename T>
void logarithmicRange(T minValue, T maxValue, size_t noValues,
                      T * result)
{
    if ( minValue/maxValue<0. ) throwError( "The first two argument in function 'logarithmicRange' (which are 'minValue' and 'maxValue') must have the same sign. " );
    if ( minValue==0. or maxValue==0. ) throwError( "None of the first two argument in function 'logarithmicRange' (which are 'minValue' and 'maxValue') must be 0. " );
    if ( noValues<=1 ) throwError( "The number of bins in function 'logarithmicRange' (which is 'noValues') must be larger than 1. " );
    double const step = pow( 10., double(std::log10( maxValue/minValue ) / (noValues-1)) );
    
    result[0] = minValue;
    for (size_t i=1; i<noValues; ++i)
        result[i] = result[i-1] * step;
}


#endif
