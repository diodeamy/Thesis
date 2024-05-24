#include <iostream>
#include <cmath>
#include <vector>

#include "response.h"
#include "threshold.h"
#include <miscellaneous.h>
#include <array.h>
using namespace std;


/* Selects the maximum response for each grid cell from two MMF response objects. It stores the results in the first MMF response object. It also stores the scale of the maximum response. */
void maximumResponse(Array<int,3> *scale,
                     ArrayReal3D *response,
                     ArrayReal3D &response2,
                     int const scale2)
{
    cout << "Extracting maximum response from scale '" << scale2 << "' ... " << flush;
    Int const gridSize = response->totalSize();
    
    for (Int i=0; i<gridSize; ++i)
        if ( response2[i]>(*response)[i] )
        {
            (*response)[i] = response2[i];
            (*scale)[i] = scale2;
        }
    cout << "Done.\n";
}




/* Rescales the MMF response such that all scales have the same weight. One can chnage the weight using a bias which gives more weigth to:
     smaller structures for 'bias'<1.
     larger structures for 'bias'>1.
*/
void rescaleResponse(ArrayReal3D &response,
                     Real const radius,
                     Real const bias)
{
    Real factor = pow( radius, Real(2.*bias) );
    Int const gridSize = response.totalSize();
    cout << "Rescaling response using 'bias'=" << bias << " (multiplication by " << factor << ") ... " << flush;
    
    for (Int i=0; i<gridSize; ++i)
        response[i] *= factor;
    
    cout << "Done.\n";
}




void maskResponse(ArrayReal3D *response,
                  Real const responseThreshold)
{
    cout << "Discarding response values smaller than the threshold " << responseThreshold << " ... " << flush;
    
    for (Int i=0; i<response->totalSize(); ++i)
        if ( (*response)[i]<responseThreshold ) (*response)[i] = 0.;
    cout << "Done.\n";
}



/* Substracts mask2 from mask1. */
void substractMask(Array<bool,3> *mask1,
                   Array<bool,3> const &mask2)
{
    if ( mask1->size()!=mask2.size() ) throwError( "You are trying to do operations with two masks of different sizes. Error in function 'substractMask'." );
    
    for (Int i=0; i<mask1->size(); ++i)
        if ( mask2.at(i) )
            (*mask1)[i] = false;
}



// /* This functions creates a mask for all the values of the response function above a certain response threshold.
//NOTE: This function discards from the mask the objects with size smaller than 'minimumObjectSize' as being noise. The value of 'minimumObjectSize' is in voxels. */
//void responseMask(ArrayReal3D &response,
//                  Real const responseThreshold,
//                  Real const minimumObjectSize,
//                  int const neighborFindingMethod,
//                  Array<bool,3> *mask)
//{
//    if ( mask->size()!=response.size() ) throwError( "You are trying to do operations with two arrays of different sizes. Error in function 'responseMask'." );
//    cout << "Computing response mask for threshold " << responseThreshold << " and discarding objects of sizes smaller than " << minimumObjectSize << ":\n" << flush;
    
//    // find how many unique objects are in the mask
//    Array<Int,1> size = response.axisSize<Int>();
//    Array<int,3> tempMask( size.ptrData(), size.totalSize() );
//    vector<int> objectSize; //keeps track of the number of grid cells in each object
//    computeCompactObjects( response, responseThreshold, neighborFindingMethod, &tempMask, &objectSize );
    
//    int minSize = int( minimumObjectSize );
//    if ( Real(minSize)<minimumObjectSize ) ++minSize;
    
//    for (Int i=0; i<tempMask.totalSize(); ++i)
//        if ( tempMask[i]>=0 and objectSize.at(tempMask.at(i))>=minSize )
//            (*mask)[i] = true;
//        else
//            (*mask)[i] = false;
//    cout << "Done.\n";
//}


// /* Discards the values of 'response' which correspond to values of the response below a certain threshold. */
//void maskResponse(ArrayReal3D *response,
//                  Real const responseThreshold,
//                  Real const minimumObjectSize,
//                  int const neighborFindingMethod)
//{
//    cout << "Discarding response values smaller than the threshold " << responseThreshold << " and also discarding objects of sizes smaller than " << minimumObjectSize << " grid cells ... \n" << flush;
    
//    // find how many unique objects are in the mask
//    Array<Int,1> size = response->axisSize<Int>();
//    Array<int,3> tempMask( size.ptrData(), size.totalSize() );
//    vector<int> objectSize; //keeps track of the number of grid cells in each object
//    computeCompactObjects( *response, responseThreshold, neighborFindingMethod, &tempMask, &objectSize );
    
    
//    int minSize = int( minimumObjectSize );
//    if ( Real(minSize)<minimumObjectSize ) ++minSize;
//    size_t noBigObjects = 0;
//    for (size_t i=0; i<objectSize.size(); ++i)
//        if ( objectSize[i]>=minSize ) ++noBigObjects;
    
//    for (Int i=0; i<tempMask.totalSize(); ++i)
//        if ( not (tempMask[i]>=0 and objectSize.at(tempMask.at(i))>=minSize) )
//            (*response)[i] = 0.;
//    cout << "Done.\n\t found " << noBigObjects << " objects bigger than " << minimumObjectSize << " grid cells.\n";
    
//}
//void maskResponseNode(ArrayReal3D *response,
//                      ArrayReal3D &density,
//                      Real const responseThreshold,
//                      Real const minimumObjectSize,
//                      int const neighborFindingMethod)
//{
//    cout << "Discarding response values smaller than the threshold " << responseThreshold << " and also discarding objects of masses smaller than " << minimumObjectSize << " average grid cells mass... \n" << flush;
    
//    // find how many unique objects are in the mask
//    Array<Int,1> size = response->axisSize<Int>();
//    Array<int,3> tempMask( size.ptrData(), size.totalSize() );
//    vector<int> objectSize; //keeps track of the number of grid cells in each object
//    Int noNodes = computeCompactObjects( *response, responseThreshold, neighborFindingMethod, &tempMask, &objectSize );
    
//    // find what is the mass of each blob identified above
    
//    Array<Real,1> massInNodes( noNodes );
//    computeMassInMask( density, tempMask, &massInNodes );
    
    
//    size_t noBigObjects = 0;
//    for (size_t i=0; i<objectSize.size(); ++i)
//        if ( massInNodes[i]>=minimumObjectSize ) ++noBigObjects;
    
//    for (Int i=0; i<tempMask.totalSize(); ++i)
//        if ( not (tempMask[i]>=0 and massInNodes(tempMask[i])>=minimumObjectSize) )
//            (*response)[i] = 0.;
//    cout << "Done.\n\t found " << noBigObjects << " objects bigger than " << minimumObjectSize << " grid cells.\n";
    
//}


