#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include <MMF_computations.h>
#include <arrayProperties.h>
#include "raw.h"
using namespace std;


void densityMMFObjects(Array<int,3> &mask,
                       Array<Real,3> &density,
                       vector<int> objectSize,
                       Array<Real,1> *averageDensity,
                       Array<Real,1> *medianDensity);



/* Counts what is the size and properties of each different MMF object. */
void outputPropertiesMMFObjects(Array<Real,3> &response,
                              Real const responseThreshold,
                              Real const minimumObjectSize,
                              int const neighborFindingMethod,
                              Array<Real,3> &density,
                              Array<Particle_pm,1> &p,
                              Real const boxSize,
                              string rootOutputFile,
                              string const programOptions)
{
    cout << "\n    Computing the properties of the MMF objects:\n" << flush;
    int minSize = int( floor(minimumObjectSize) );
    
    // Find the MMF objects above a given threshold
    Array<Int,1> size = response.axisSize<Int>();
    Array<int,3> mask( size.ptrData(), size.totalSize() );
    vector<int> objectSize; //keeps track of the number of grid cells in each object
    computeCompactObjects( response, responseThreshold, neighborFindingMethod, &mask, &objectSize ); // finds the number of distinct objects
    
    // Label the objects according to their volume size
    relabelMMFObjects( &mask, &objectSize );
    
    // Compute the average density of the MMF objects
    Array<Real,1> averageDensity( objectSize.size() );
    Array<Real,1> medianDensity( objectSize.size() );
    densityMMFObjects( mask, density, objectSize, &averageDensity, &medianDensity );
    
    // Compute the number of particles in each MMF object
    Array<Real,1> noParticlesInObjects( objectSize.size() );
    computeParticlesInMask( p, boxSize, mask, &noParticlesInObjects );
    
    
    // Compute properties of the object volume distribution
    ArrayProperties<int,size_t> propVolume( &(objectSize[0]), objectSize.size() );
    double noValidCells = propVolume.sum();     // returns the number of cells with response >= responseThreshold
    double gridCellVolume = pow( boxSize/MPC_UNIT, NO_DIM ) / response.size();
    
    
    // now 'objectSize' stores the size, 'averageDensity' the average density while 'noParticlesInObjects' stores the number of particles in each MMF object
    // Let us output this properties for each object
    string output1 = rootOutputFile + ".MMFprop";
    cout << "Writting the properties of the MMF objects to the ASCII file '" << output1 << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, output1 );   //open the output file
    
    outputFile << "# The file contains the properties of MMF objects for a given response threshold:\n"
            << "#\t 1st column: the MMF object label.\n"
            << "#\t 2nd column: the MMF object volume as a number of grid cells.\n"
            << "#\t 3rd column: the MMF object volume as a physical volume (in Mpc^3).\n"
            << "#\t 4th column: the MMF object average density.\n"
            << "#\t 5th column: the MMF object median density.\n"
            << "#\t 6th column: the number of particles in the MMF object.\n"
            << "#\t 7th column: ration b/a when fitting an ellipsoid to the MMF object.\n"
            << "#\t 8th column: ration c/a when fitting an ellipsoid to the MMF object.\n"
            << "# The results were obtained using the command:  " << programOptions << "\n"
            << "# Number of MMF objects: " << objectSize.size() << "\n"
            << "# Number of cells with response larger than threshold: " << noValidCells << " which repsesents " << setprecision(4) << noValidCells/response.size()*100. <<  setprecision(6) <<"\% of the total number of grid cells.\n"
            << "\n";
    
    for (size_t i=0; i<objectSize.size(); ++i)
        if ( objectSize[i]>=minSize )
            outputFile << i << "\t" << objectSize[i] << "\t" <<  objectSize[i]*gridCellVolume << "\t" << averageDensity[i] << "\t" << /*medianDensity[i] << "\t" <<*/ noParticlesInObjects[i] << "\n";
    
    cout << "Done.\n";
    outputFile.close();
}


/* Outputs the different MMF objects labeled according to their volume in descending order. The output is a binary file with the same grid size as the response grid. */
void outputMMFObjects(Array<Real,3> &response,
                      Real const responseThreshold,
                      Real const minimumObjectSize,
                      int const neighborFindingMethod,
                      string nameOutputFile,
                      string const programOptions)
{
    cout << "\n    Outputing the MMF objects labeled according to their volume:\n" << flush;
    int minSize = int( floor(minimumObjectSize) );
    
    // Find the MMF objects above a given threshold
    Array<Int,1> size = response.axisSize<Int>();
    Array<int,3> mask( size.ptrData(), size.totalSize() );
    vector<int> objectSize; //keeps track of the number of grid cells in each object
    computeCompactObjects( response, responseThreshold, neighborFindingMethod, &mask, &objectSize ); // finds the number of distinct objects
    
    // Label the objects according to their volume size
    relabelMMFObjects( &mask, &objectSize );
    
    // keep only objects >= minSize
    int maxObjects = -1;  //maximum number of outputed objects (i.e. number of objects larger than minSize)
    for (size_t i=0; i<objectSize.size(); ++i)
        if ( objectSize[i]<minSize )
        {
            maxObjects = i;
            break;
        }
    if ( maxObjects>=0 )
        for (Int i=0; i<mask.size(); ++i)
            if ( mask[i]>=maxObjects ) mask[i] = -1;
    
    
    // print the result to a binary file
    outputRaw( mask.ptrData(), mask.size(), nameOutputFile );
}






/* Finds the average density of MMF objects. */
void densityMMFObjects(Array<int,3> &mask,
                       Array<Real,3> &density,
                       vector<int> objectSize,
                       Array<Real,1> *averageDensity,
                       Array<Real,1> *medianDensity)
{
    if ( objectSize.size()!=averageDensity->size() ) throwError( "The 'objectSize' and 'averageDensity' arguments of function 'densityMMFObjects' have different sizes. This must not be the case." );
    if ( mask.size()!=density.size() ) throwError( "The 'mask' and 'density' arguments of function 'densityMMFObjects' have different sizes. This must not be the case." );
    
    cout << "Computing the average and median density of each MMF object ... " << flush;
    // allocate memory for storing the density for each object
    Array<int,1> count( objectSize.size() );
    count.assign(0);
//     Real **objectDensity;
//     objectDensity = new Real*[objectSize.size()];
//     for (size_t i=0; i<objectSize.size(); ++i)
//         objectDensity[i] =  new Real[objectSize[i]];
    
    
    averageDensity->assign( Real(0.) );
    for (size_t i=0; i<mask.size(); ++i)
        if ( mask[i]>=0 )
        {
            (*averageDensity)[ mask[i] ] += density[i];
//             objectDensity[ mask[i] ][ count[mask[i]] ] = density[i];
            ++count[mask[i]];
        }
    
    // normalize density by volume of each MMF object
    for (size_t i=0; i<objectSize.size(); ++i)
        if ( objectSize[i]!=0 )
            (*averageDensity)[i] /= objectSize[i];
    
//     // compute the median density for each object
//     for (size_t i=0; i<objectSize.size(); ++i)
//     {
//         ArrayProperties<Real,size_t> prop( objectDensity[i], objectSize[i] );
//         (*medianDensity)[i] = prop.medianApproximative( 100, LOGARITHMIC_BIN );
//         delete[] objectDensity[i];
//     }
//     delete[] objectDensity;
    
    cout << "Done.\n";
}



