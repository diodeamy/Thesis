
#ifndef MMFCONTRACTION_HEADER
#define MMFCONTRACTION_HEADER


#include <vector>
#include <utility>
#include <iostream>

#include <defines.h>
#include <box.h>
#include <misc.h>
#include <array.h>
#include <box.h>

using namespace std;



/* Convergence parameters and other setings for the filament/wall contraction functions. */
struct ContractionOptions
{
    // overall parameters
    Box<Real,NO_DIM> box;       // keeps track of the dimensions of the box encompassing the data
    Real    length[NO_DIM];     // the length of the box along each dimension
    int     grid[NO_DIM];       // the size of the grid along each dimension
    bool    periodic;           // true if the data is in a periodic box - will apply periodic boundaries
    Real    radius;             // the radius around which to search for neighbors and do the contraction process (should contain at least a few grid cells)
    int     feature;            // 3=filament, 2=wall - tells the function what environment they are contracting 
    
    
    // parameters for the contraction of the objects once each cell has a filament/wall direction
    Real    convergeFraction;       // stop the iteration process to contract the volume to a spine/wall when at most this fraction did not converge to a final position
    Real    distanceThreshold;      // points converged when will move less than this distance in the next step
    Real    eigenvalueThreshold;    // points converged if the shape of the cloud they are imbedded in is very anisotropic (i.e. lambda_2/lambda_3 < eigenvalueThreshold for filaments and lambda_1/lambda_2 < eigenvalueThreshold for walls - we have lambda_1 <= lambda_2 <= lambda_3 and the fila direction is given by eigenvector_3 while the wall normal is given by eigenvector_1 )
    int     maxLoop;                // maximum loop - automatically stop after this many iterations indifferently if convergence is achieved
    int     maxIterations;          // the number of maximum times to iterate the procedure for direction identification (if convergence not reached before)
    Real    cosDirectionThreshold;  // point direction converged when the filament/wall direction associated to them changed by less than this angle from the previous step
    
    
    ContractionOptions()
    {
        periodic = false;
        radius   = 1.;
        feature  = -1;
        
        convergeFraction = 0.01;
        distanceThreshold = 0.1;
        eigenvalueThreshold = 0.3;
        maxLoop = 60;
        maxIterations = 20;
        cosDirectionThreshold = 0.996;  // within 5 degrees of expected direction
    }
    
    template <typename T>
    void addBoxCoordinates(T *coords)
    {
        for (int i=0; i<2*NO_DIM; ++i)
            box[i] = Real(coords[i]);
        if ( not box.validBox() )
            throwError( "The box coordinates supplied to function 'ContractionOptions::addBoxCoordinates' are not valid. They failed the 'Box<T,N>::validBox()' test." );
        for (int i=0; i<NO_DIM; ++i)
            length[i] = box[2*i+1] - box[2*i];
    }
};


// Structure to keep track of the data that has to be changed between all the different functions
struct Datas
{
    Real *spinePosition;
    Real *spineDirection;
    int  objectIndex;
    vector<Real>  neighbors;
    vector<Real>  directions;
    bool updateDirection;
    
    Datas() {objectIndex=-1; updateDirection=false;}
};






// returns the total number of valid 'object' cells
size_t noObjectCells(Array<int,1> &objectSize);
size_t noObjectCells(Array<int,3> &mask);
size_t noObjectCells(Array<shortInt,3> &mask);


// returns the cell indices and positions corresponding to valid object grid cells
void validObjectCell(Array<shortInt,3> &mask,
                     Array<int,2> &cellIndices,
                     Array<Real,2> &spinePosition,
                     Box<Real,NO_DIM> &box);
void validObjectCell(Array<shortInt,3> &mask,
                     Array<Real,2> &spinePosition,
                     Box<Real,NO_DIM> &box);




/* This function computes the inertia tensor of the cloud of points around a given point. */
size_t directionDetection(Array<Real,2> &spinePosition,
                          Array<Real,2> &oldDirections,
                          Array<Real,2> &directions,
                          ContractionOptions &options,
                          bool const VERBOSE = true);


/* Contracts the filaments/walls to their spine/principal wall. */ 
void MMFContraction(Array<Real,2> &spinePosition,
                    Array<Real,2> &directionVector,
                    ContractionOptions options,
                    Real const convergeFraction,
                    Real const distanceThreshold,
                    Real const eigenvalueThreshold,
                    int const maxLoop);
void MMFContraction(Array<Real,2> &spinePosition,
                    ContractionOptions options,
                    Real const convergeFraction,
                    Real const distanceThreshold,
                    Real const eigenvalueThreshold,
                    int const maxLoop);


/* Uses an iterative method to get the direction of the filament/wall. */
void MMFObjectDirection(Array<Real,2> &spinePosition,
                        Array<Real,2> &directionVector,
                        ContractionOptions &options,
                        bool hasInitialDirection = false);



// function that computes different quantities along the filament and plane of the wall
void environmentProperties(Array<Real,1> &density,
                           Array<Real,2> &spinePosition,
                           ContractionOptions &options,
                           Array<Real,2> &properties,
                           Array<double,2> &sizeData,
                           Array<double,2> &massData,
                           Array<Real,1> &averageData);

#endif
