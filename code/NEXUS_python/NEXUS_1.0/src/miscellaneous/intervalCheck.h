#ifndef INTERVALCHECK_HEADER
#define INTERVALCHECK_HEADER

#include "throwError.h"

/* Checks if a number is within a given interval. */
template <typename T1>
void intervalCheck(T1 const target,
                   T1 const minValue, T1 const maxValue,
                   std::string errorMessage)
{
    if ( target<minValue or target>maxValue )
    {
        Error error;
        error << "Some program variables failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be larger than " << minValue << " and smaller than " << maxValue << "." ;
        error.EndErrorMessage();
    }
}

/* Checks if a number is > than a lower bound. */
template <typename T1>
void lowerBoundCheck(T1 const target,
                     T1 const minValue,
                     std::string errorMessage)
{
    if ( target<minValue )
    {
        Error error;
        error << "Some program variables failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be larger or equal than " << minValue << ".";
        error.EndErrorMessage();
    }
}


/* Checks if a number is < than an upper bound. */
template <typename T1>
void upperBoundCheck(T1 const target,
                     T1 const maxValue,
                     std::string errorMessage)
{
    if ( target>maxValue )
    {
        Error error;
        error << "Some program variables failed a consitency check. The variable " << errorMessage << " has the value " << target << ", but it should be smaller or equal than " << maxValue << ".";
        error.EndErrorMessage();
    }
}

#endif
