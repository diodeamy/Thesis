#ifndef RESPONSE_HEADER
#define RESPONSE_HEADER

#include <defines.h>
#include <array.h>
#include <vector.h>


typedef Array<Real,NO_DIM>             ArrayReal3D;
typedef Vector<Real,NO_DIM>            VectorReal3D;    // real vector with 3 coordinates
typedef Array<VectorReal3D,NO_DIM>     ArrayVector3D;


void maximumResponse(Array<int,3> *scale,
                     ArrayReal3D *response,
                     ArrayReal3D &response2,
                     int const scale2);

void rescaleResponse(ArrayReal3D &response,
                         Real const radius,
                         Real const bias);

void responseMask(ArrayReal3D &response,
                  Real const responseThreshold,
                  Real const minimumObjectSize,
                  int const neighborFindingMethod,
                  Array<bool,3> *mask);

void maskResponse(ArrayReal3D *response,
                  Real const responseThreshold,
                  Real const minimumObjectSize,
                  int const neighborFindingMethod);
void maskResponse(ArrayReal3D *response,
                  Real const responseThreshold);
void maskResponseNode(ArrayReal3D *response,
                      ArrayReal3D &density,
                      Real const responseThreshold,
                      Real const minimumObjectSize,
                      int const neighborFindingMethod);


void substractMask(Array<bool,3> *mask1,
                   Array<bool,3> const &mask2);



/* This functions creates a mask for all the values of the response function above a certain response threshold.
Can what values to use in the output array (i.e. what it means to be):
    valid = set values to 'valid'
    invalid = set values to 'invalid'
*/
template <typename T>
void responseMask(ArrayReal3D &response,
                  Real const responseThreshold,
                  Array<T,NO_DIM> *mask,
                  T const valid = T(1), T const invalid = T(0))
{
    if ( mask->size()!=response.size() )
        throwError( "You are trying to do operations with two arrays of different sizes. Error in function 'responseMask'." );
    if ( valid==invalid )
        throwError( "In function 'responseMask'. Both the 'valid' and 'invalid' parameters have the same value. They should have different values to distinguish between valid and invalid grid cells." );
    
    size_t const gridSize = response.size();
    for (size_t i=0; i<gridSize; ++i)
        if ( response[i]>=responseThreshold )
            (*mask)[i] = valid;
    else
        (*mask)[i] = invalid;
}

/* Discards the values of 'response' which correspond to values of 'maskResponse' above the threshold. */
template <typename T>
void maskResponse(ArrayReal3D *response,
                  Array<T,NO_DIM> &maskResponse,
                  T const responseThreshold)
{
    if ( response->size()!=maskResponse.size() )
        throwError( "You are trying to do operations with array of different sizes. Error in function 'maskResponse'." );
    
    for (size_t i=0; i<response->size(); ++i)
        if ( maskResponse[i]>=responseThreshold )
            (*response)[i] = Real(0.);
}

// /* Discards the values of the 'inputField' using 'value' which correspond to values of 'maskResponse' above the threshold. */
// template <typename T>
// void maskInputField(ArrayReal3D *inputField,
//                     Real const value,
//                     Array<T,NO_DIM> &maskResponse,
//                     T const responseThreshold)
// {
//     if ( inputField->size()!=maskResponse.size() )
//         throwError( "You are trying to do operations with array of different sizes. Error in function 'maskInputField'." );
//     
//     for (size_t i=0; i<inputField->size(); ++i)
//         if ( maskResponse[i]>=responseThreshold )
//             (*inputField)[i] = value;
// }



#endif
