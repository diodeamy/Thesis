

#include "threshold.h"
using namespace std;


/* Computes what are the distinct objects by adding to the existing mask the voxels which have a response value between 'minThreshlod' <= response < 'maxThreshold'.
Input parameters to the function:
        response - the response map
        minThreshold, maxThreshold - selects the voxels that will be add on top of the existing mask
        mask - the existing mask where each individual object has a label
        neighborFindingMethod - what are the neighbors 1=6 neighbors, 2=26 neighbors
        mass - the array giving the mass in each grid cell
        minimumSize - the minimum size of a significant object (volume in grid cells)
        properties - vector that stores the return results.
        trackStatistics - 2 component int array with: [0] - number of significant objects &&  [2] - total number of objects
The return results are:
                [0] - total mass in significant objects
                [1] - total volume in significant objects
                [2] - fraction of volume of largest object compared to the total volume in significant objects
                [3] - total mass in all objects
                [4] - total volume in all objects
                [5] - fraction of volume of largest object compared to the total volume in all objects
*/
void thresholdProperties(Array<Real,3> &response,
                         Real const minThreshold,
                         Real const maxThreshold,
                         Array<int,3> *mask,
                         int const neighborFindingMethod,
                         Array<Real,3> &mass,
                         Int const minimumSize,
                         Real *properties,
                         int *trackStatistics)
{
    // get the voxels that are added by including the cells with 'minThreshlod' <= response < 'maxThreshold'
    vector<Int> newVoxels;  //vector to store the grid indices for the new voxels
    newVoxels.reserve( response.size()/100 );   //a reasonable amount to reserve for the memory requierments
    selection( response, minThreshold, maxThreshold, &newVoxels );
    
    // compute the new distinct objects using the additional voxels in 'newVoxels'
    vector<Int> objectVolume;
    identifyDistinctObjects( mask, trackStatistics[1], newVoxels, neighborFindingMethod, &objectVolume );
    
//    // relabel objects according to their size
//    relabelObjects( mask, &objectVolume, false );
    
    // get the mass associated to each object
    vector<double> objectMass;
    objectsMass( *mask, mass, &objectMass, objectVolume.size() );
    
    // loop over the objects and compute the output results
    int noSignificantObjects = 0;
    Int totalVolume = 0, significantVolume = 0, largestObject = 0;
    double totalMass = 0., significantMass = 0.;
    for (size_t i=0; i<objectVolume.size(); ++i)
    {
        totalVolume += objectVolume[i];
        totalMass   += objectMass[i];
        if ( objectVolume[i]>=minimumSize )
        {
            significantVolume += objectVolume[i];
            significantMass   += objectMass[i];
            if ( objectVolume[i]>largestObject ) largestObject = objectVolume[i];
            ++noSignificantObjects;
        }
    }
    Real volumeFraction = significantVolume==0 ? 1. : Real(largestObject)/significantVolume;   // returns the volume fraction of largest filament compared to all significant filaments
    properties[0] = significantMass;                            //total mass in significant objects
    properties[1] = significantVolume;                          //total volume in significant objects
    properties[2] = volumeFraction;             //fraction of volume of largest object compared to the total volume in significant objects
    properties[3] = totalMass;                                  //total mass in all objects
    properties[4] = totalVolume;                                //total volume in all objects
    properties[5] = Real(largestObject) / Real(totalVolume);    //fraction of volume of largest object compared to the total volume in objects
    trackStatistics[0] = noSignificantObjects;  //number of significant objects
    trackStatistics[1] = objectVolume.size();   //the total number of objects
}





/* Computes the percentage of the objects that have a density higher than the virial density. See the above function for what each function argument means.
NOTE: This function can be used to get the optimal threshold for nodes. */
Real fractionVirialNodes(Array<Real,3> &response,
                         Real const minThreshold,
                         Real const maxThreshold,
                         Array<int,3> *mask,
                         int const neighborFindingMethod,
                         Array<Real,3> &mass,
                         double const minSize,
                         double const virialDensity,
                         Real *properties,
                         int *trackStatistics)
{
    // get the voxels that are added by including the cells with 'minThreshlod' <= response < 'maxThreshold'
    vector<Int> newVoxels;  //vector to store the grid indices for the new voxels
    newVoxels.reserve( response.size()/100 );   //a reasonable amount to reserve for the memory requierments
    selection( response, minThreshold, maxThreshold, &newVoxels );
    
    
    // compute the new distinct objects using the additional voxels in 'newVoxels'
    vector<Int> objectVolume;
    identifyDistinctObjects( mask, trackStatistics[1], newVoxels, neighborFindingMethod, &objectVolume );
    
    
    // get the mass associated to each object
    vector<double> objectMass;
    objectsMass( *mask, mass, &objectMass, objectVolume.size() );
    
    
    // loop over the objects and compute the output results
    int noSignificantObjects = 0, noValidNodes = 0; //number of objects larger than 'minSize' and number of objects with averageDensity > virialDensity
    Int totalVolume = 0, significantVolume = 0;
    double totalMass = 0., significantMass = 0.;
    for (size_t i=0; i<objectVolume.size(); ++i)
    {
        totalVolume += objectVolume[i];
        totalMass   += objectMass[i];
        if ( objectMass[i]>=minSize )
        {
            significantVolume += objectVolume[i];
            significantMass   += objectMass[i];
            ++noSignificantObjects;
            double tempDensity = objectMass[i] / objectVolume[i];
            if ( tempDensity>=virialDensity )
                ++noValidNodes;
        }
    }
    Real validFraction = noSignificantObjects==0 ? 1. : (noValidNodes==0 ? 0. : Real(noValidNodes)/noSignificantObjects);   // returns the fraction of nodes with density > virialDensity (returns 0 for very small thresholds and 1 for very large thresholds)
    properties[0] = significantMass;                            //total mass in significant objects
    properties[1] = significantVolume;                          //total volume in significant objects
    properties[2] = validFraction;                              //fraction of significant nodes with average density > virialDensity
    properties[3] = totalMass;                                  //total mass in all objects
    properties[4] = totalVolume;                                //total volume in all objects
    properties[5] = Real(noSignificantObjects);                 //number of significant objects
    trackStatistics[0] = noSignificantObjects;  //number of significant objects
    trackStatistics[1] = objectVolume.size();   //the total number of objects
    
    return validFraction;
}





// /* Computes the number of different objects for a given threshold. Returns the total number of objects.
//NOTE: This function can be used to get the threshold for filaments and walls. */
//Real countObjects(Array<Real,3> &response,
//                  Real const threshold,
//                  int const neighborFindingMethod,
//                  Real const minimumSize,
//                  size_t *sizeLargestObject,
//                  size_t *numberLargeObjects,
//                  Real *volumeFraction)
//{
//    // find how many unique objects are in the mask
//    Array<Int,1> size = response.axisSize<Int>();
//    Array<int,3> mask( size.ptrData(), size.totalSize() );
//    vector<int> objectSize; //keeps track of the number of grid cells in each object
//    computeCompactObjects( response, threshold, neighborFindingMethod, &mask, &objectSize ); // returns an array with the object number for each grid cell (if !=-1 which is empty cell) and the number of grid cells in each object
    
//    // count how many of the objects have a size larger than 'minimumSize'
//    size_t noObjects = 0;
//    size_t volumeLarge = 0, volumeTotal = 0;
//    for (size_t i=0; i<objectSize.size(); ++i)
//    {
//        volumeTotal += objectSize[i];
//        if ( Real(objectSize[i])>=minimumSize )
//        {
//            ++noObjects;
//            volumeLarge += objectSize[i];
//        }
//    }
//    *numberLargeObjects = noObjects;
//    *sizeLargestObject = Max( &(objectSize[0]), objectSize.size() );
//    *volumeFraction = Real(volumeTotal) / ( size[0]*size[1]*size[2] ) ;
//    return (*sizeLargestObject)/Real(volumeLarge);
//}





// /* Finds the number of particles in each 'object' (here 'object' is defined as an area of the mask with all the cells connected). */
//void computeParticlesInMask(Array<Particle_pm,1> &p,
//                            Real const boxSize,
//                            Array<int,3> &mask,
//                            Array<float,1> *massInBlobs)
//{
//    if ( NO_DIM!=3 ) throwError( "The function 'computeParticlesInMask' gives correct results only for 3 dimensions." );
////     cout << "Computing the number of particles in MMF objects for the given threshold ... " << flush;
//    // define some constants
//    Real const dx[3] = { boxSize / mask.getSize(0), boxSize / mask.getSize(1), boxSize / mask.getSize(2) };
//    Int const nZ = mask.getSize(2);
//    Int const nYnZ = mask.getSize(1) * nZ;
//    massInBlobs->assign(0.); //reset number of particles in each blob to 0
    
    
//    for (Int i=0; i<p.totalSize(); ++i)
//    {
//        int temp[NO_DIM];
//        for (int j=0; j<NO_DIM; ++j)
//            temp[j] = int( p[i].pos[j] / dx[j] );
//        Int index = nYnZ*temp[0] + nZ*temp[1] + temp[2];
        
//        if ( mask(index)>=0 )
//            (*massInBlobs)[ mask(index) ] += 1.;//p[i].mass;
//    }
//    Real max2 = 0.;
//    for (size_t i=0; i<massInBlobs->size(); ++i)
//        max2 += massInBlobs->at(i);
    
////     cout << "Done.\n";
////     cout << "\t there are " << max2 << " particles in MMF objects which represent " << max2/p.size()*100. << "\% of the total number of particles.\n";
//}









