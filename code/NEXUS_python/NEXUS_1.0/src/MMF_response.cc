#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <omp.h>

#include <defines.h>
#include "densityFile/densityFile.h"
#include "densityFile/density_header.h"
#include "MMF_file/MMF_file.h"
#include "MMF_file/MMF_header.h"
#include "miscellaneous/miscellaneous.h"
#include "MMF_computations/MMF_computations.h"
#include <array.h>
#include <vector.h>
#include <k_space.h>
#include "classes/FFTW.h"
using namespace std;



#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;


struct Environment
{
    bool   envOn;         // true if the user asked to compute this environment
    string outputFile;    // name of the output file
    int    feature;       // cosmic web feature (4=blob, 3=filament, 2=wall)
    int    filter;        // filter type to be applied (4x for blobs, 3x for filaments and 2x for walls)
};

struct User_options
{
    string inputFile;     // name of input density file
    string scale;         // scales for which to compute the MMF response function
    vector<int>   scales; // the scales inserted by the user
    
    Environment node;     // keeps track of node information
    Environment fila;     // keeps track of filament information
    Environment wall;     // keeps track of wall information
    
    bool   radiaOn;       // true if user inserted radia for each scale
    vector<float> radia;  // values of radia inserted by the user
    float  radius0;       // the value of the filter radius for scale 0 (in Mpc)
    float  base;          // the value of the base which raised to the scale 'n' gives the filter radius at that scale (radius = radius0 * base^n)
    
    float  bias;          // the value of the bias, if any
    bool   biasOn;        // true if the user supplied a bias value
    
    string nodeFile;      // the clean MMF response of nodes used for filament and wall computations
    string filaFile;      // the clean MMF response of filaments used for wall computations
    bool   mask;          // true if to use a mask for the fillament and wall computations
    
    bool   densityOn;     // true if to compute the MMF response using the density field
    bool   logarithmOn;   // true if to compute the MMF response using log(density)
    bool   gravitationOn; // true if to compute the MMF response using the gravitational field
    bool   logarithmFiltering;     //true if to compute the MMF response using gaussian filtering in logarithm space
    bool   velocityOn;    // true if to compute the MMF response using the velocity divergence
    bool   velocityPotentialOn;    // true if to compute the MMF response using the velocity potential
    
    bool   eigenOn;       // true if to compute the eigenvalues and eigenvectors corresponding to the maximum response
    bool   allEigen;      // true if to output all the eigenvectors and eigenvalues
    bool   singleEigen;   // true if to output only the eigenvector corresponding to the given LSS feature (works only for filaments and walls)
    
    bool   maximumScaleOn;// true if to output also the scale at which the maximum response was found
    
    bool boxCoordinatesOn;// true if the user suplied the box coordinates along the 3 directions
    vector<Real> boxCoordinates; //stores xMin, xMax, yMin, yMax, zMin and zMax for the box in question
    
    
    string programOptions;// string that stores the 'argv' values for the users information
    
    
    User_options()
    { 
    node.envOn = false; node.feature = 4; node.filter=40;
    fila.envOn = false; fila.feature = 3; fila.filter=30;
    wall.envOn = false; wall.feature = 2; wall.filter=20;
    radiaOn=false; radius0=0.5; base=sqrt(2.);
    bias = 1.; biasOn = false;
    mask = false;
    densityOn=false; logarithmOn=false; gravitationOn=false;  logarithmFiltering=false; velocityOn=false; velocityPotentialOn=false;
    eigenOn=false; allEigen = false; singleEigen = false;
    maximumScaleOn = false;
    boxCoordinatesOn = false;
    programOptions = "";}
};
void readOptions(int argc, char *argv[],
                 User_options *userOptions);            // read the options inserted by the user

void logarithmFilteringResponse(User_options &userOptions,
                                ArrayReal3D &density,
                                Density_header &densityHeader,
                                MMF_header &mmfHeader); // computes the response using the logarithm filtering method

void outputResults(User_options &userOptions,
                   MMF_header  &mmfHeader,
                   ArrayReal3D &nodeResponse,
                   ArrayReal3D &filaResponse,
                   ArrayReal3D &wallResponse,
                   Array<int,3>  &maxScale,
                   ArrayVector3D &eigenValues,
                   ArrayVector9D &eigenVectors);    // outputs the results to file






int main( int argc, char *argv[] )
{
    // read the user options
    User_options userOptions;
    readOptions( argc, argv, &userOptions );
    
    
    //! Read the input data
    // read the density file header
    Density_header densityHeader;
    readDensityHeader( &densityHeader, userOptions.inputFile );
    
    // check that the input data is compatible with the program options
    if ( densityHeader.fileType!=DENSITY_FILE and densityHeader.fileType!=VELOCITY_DIVERGENCE_FILE )
        throwError( "The program expected as input a density map or velocity divergence map (file type 1 and 13), but it recieved as input a file type '", densityHeader.fileType, "'. Cannot continue without the correct file type!" );
    if ( densityHeader.fileType==DENSITY_FILE and (userOptions.velocityOn or userOptions.velocityPotentialOn) )
        if ( userOptions.velocityOn ) throwError( "The input file contains a density map. You cannot compute the MMF response using the velocity divergence field or related quantity using the given input data." );
    else if ( densityHeader.fileType==VELOCITY_DIVERGENCE_FILE and ( userOptions.densityOn or userOptions.logarithmOn or userOptions.gravitationOn or userOptions.logarithmFiltering ) )
        throwError( "The input file contains a velocity divergence map. You cannot compute the MMF response using the density field or related quantity using the given input data." );
    
    // read the data in the file
    ArrayReal3D density( densityHeader.gridSize, NO_DIM );
    readDensityData( &density, userOptions.inputFile );
    for (size_t i=0; i<density.size(); ++i)
        if ( not boost::math::isfinite(density[i]) )
        {
            density[i] = Real(0.);
            cout << "Nan or inf at grid cell =" << i << "  replaced with value =" << density[i] << "\n";
        }
    
    
    // define the header for the output file
    MMF_header mmfHeader;
    mmfHeader.copyDensityHeader( densityHeader );      //copy common information to MMF header
    mmfHeader.overwriteObservations( argv, argc );
    mmfHeader.method = MMF_METHOD_DENSITY;
    // update data box if suplied by the user
    if ( userOptions.boxCoordinatesOn )
        for (int i=0; i<6; ++i)
            mmfHeader.box[i] = userOptions.boxCoordinates[i];
    Real boxLength[3];
    for (int i=0; i<3; ++i)
        boxLength[i] = mmfHeader.box[2*i+1] - mmfHeader.box[2*i];
    
    
    
    //! Apply the node and filament masks (if required)
    Real minValue = Real(0.);//Min( density.ptrData(), density.size() );
    if ( userOptions.mask and (userOptions.fila.envOn or userOptions.wall.envOn) )  //apply node mask
        maskInputField( &density, mmfHeader, minValue, userOptions.nodeFile, "node" );
    if ( userOptions.mask and userOptions.wall.envOn )  //apply filament mask
        maskInputField( &density, mmfHeader, minValue, userOptions.filaFile, "filament" );
    
    
    
    //! Compute from the input file the quantities requested by the user
    // call an external function for the logarithm filtering method
    if ( userOptions.logarithmFiltering )
    {
        logarithmFilteringResponse(userOptions, density, densityHeader, mmfHeader );
        return 0;   //exit the program
    }
    
    // if the input field is a velocity divergence, multiply by -1 since the density maxima correspond to minima in the velocity divergence
    if ( densityHeader.fileType==VELOCITY_DIVERGENCE_FILE and userOptions.velocityOn )
    {
        cout << "Computing the absolute value of the velocity divergence. This is used when computing the environments from the velocity divergence.\n" << flush;
        for (size_t i=0; i<density.size(); ++i)
            density[i] *= Real(-1.);
//            if ( density[i]<Real(0.) ) density[i] *= Real(-1.);
        mmfHeader.method = MMF_METHOD_VELOCITY;
    }
    
    
    // if logarithm options, first take the logarithm of the density
    if ( userOptions.logarithmOn and densityHeader.fileType==DENSITY_FILE )
    {
        densityLogarithm( &density, Real(1.e-3) ); // computes the density logarithm and outputs the values to the 'density' array
        mmfHeader.method = MMF_METHOD_DENSITY_LOGARITHM;
    }
    
    
    // Compute the gravitational field, if computing the MMF using the gravitational field and not the density field
    if ( userOptions.gravitationOn and densityHeader.fileType==DENSITY_FILE  )
    {
        if ( not userOptions.biasOn )
        {
            userOptions.biasOn = true;
            userOptions.bias = Real(0.);
            cout << "No bias detected. Switching to the default bias for NEXUS gravitational field computations: bias = 0.\n" << flush;
        }
        
        //the negative of the gravitational field to keep the sign of the eigenvalues the same as for the density Hessian
        solvePoisson( &density, boxLength, -1., "gravitational field from the density field" );
        mmfHeader.method = MMF_METHOD_GRAVITY;
    }
    // Compute the velocity potential, if computing the MMF using the velocity potential
    else if ( userOptions.velocityPotentialOn and densityHeader.fileType==VELOCITY_DIVERGENCE_FILE )
    {
        if ( not userOptions.biasOn )
        {
            userOptions.biasOn = true;
            userOptions.bias = Real(0.);
            cout << "No bias detected. Switching to the default bias for NEXUS velocity potential field computations: bias = 0.\n" << flush;
        }
        
         //the "+" sign comes such that the eigenvalues of the velocity potential have the same sign as that of the density
        solvePoisson( &density, boxLength, +1., "velocity potential field from the velocity divergence field" );
        mmfHeader.method = MMF_METHOD_VELOCITY_POTENTIAL;
    }
    
    
    
    //! Compute the DFT of the field used for feature identification
    ArrayComplex3D denComplex( densityHeader.gridSize[0], densityHeader.gridSize[1], densityHeader.gridSize[2]/2+1 );
    computeFFTW( density, denComplex, FFTW_class::R2C, omp_get_max_threads() );
    density.freeMemory();   // free the memory allocated to the density data - don't need it anymore
    
    
    
    //! define the output quantities
    ArrayReal3D nodeResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for nodes
    if ( not userOptions.node.envOn ) nodeResponse.freeMemory();//release the memory if no node computations
    ArrayReal3D filaResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for filaments
    if ( not userOptions.fila.envOn ) filaResponse.freeMemory();//release the memory if no filament computations
    ArrayReal3D wallResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for walls
    if ( not userOptions.wall.envOn ) wallResponse.freeMemory();//release the memory if no wall computations
    
    Array<int,3> scale( densityHeader.gridSize, NO_DIM );       // will keep track of the scale corresponding to the maximum response value
    scale.assign(userOptions.scales[0]);                        // initialize all the values to the current value of the scale
    ArrayVector3D eigenvalues( densityHeader.gridSize, NO_DIM );       // will store the eigenvalues corresponding to the maximum response
    ArrayVector9D eigenvectors( densityHeader.gridSize, NO_DIM );      // will store the eigenvectors corresponding to the maximum response (will store only the first NO_DIM-1 eigenvectors)
    
    // define additional variable to keep track of the response, eigenvalues and eigenvectors for each scale
    ArrayReal3D tempResponse( densityHeader.gridSize, NO_DIM );  // will store response for each scale != scale[0]
    ArrayVector3D *eValues;         // pointer to eigenvalues for each scale != scale[0]
    ArrayVector9D *eVectors;        // pointer to eigenvectors for each scale != scale[0]
    if ( userOptions.eigenOn )      // if the user asked for the eigenvalues and eigenvectors corresponding to the maximum response, need to allocate additional memory for that
    {
        eValues = new ArrayVector3D( densityHeader.gridSize, NO_DIM );
        eVectors = new ArrayVector9D( densityHeader.gridSize, NO_DIM );
    }
    else
    {
        eValues = &eigenvalues;
        eVectors = &eigenvectors;
    }
    
    
    
    //! loop over the scales - compute the MMF response and output the results
    // compute the quantities for the 1st scale
    cout << "\nComputations for scale: " << userOptions.scales[0] << "\n";
    
    // get the eigenvalues and eigenvectors of the Hessian matrix 
    hessianMatrix( denComplex, eigenvalues, eigenvectors, userOptions.radia[0], boxLength, userOptions.bias, userOptions.eigenOn );
    
    //compute the response for the given environments
    if ( userOptions.node.envOn )
        computeResponse( eigenvalues, nodeResponse, userOptions.node.filter );
    if ( userOptions.fila.envOn )
        computeResponse( eigenvalues, filaResponse, userOptions.fila.filter );
    if ( userOptions.wall.envOn )
        computeResponse( eigenvalues, wallResponse, userOptions.wall.filter );
    
    //compute the quantities for the rest of the scales
    for (size_t i=1; i<userOptions.scales.size(); ++i)
    {
        cout << "\nComputations for scale: " << userOptions.scales[i] << "\n";
        
        // find the eignevalues and eigenvectors of the hessian matrix
        hessianMatrix( denComplex, *eValues, *eVectors, userOptions.radia[i], boxLength, userOptions.bias, userOptions.eigenOn );
        
        // compute response for the given scale and the maximum response
        if ( userOptions.node.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.node.filter );
            maximumResponse( &scale, &nodeResponse, tempResponse, userOptions.scales[i] );
        }
        if ( userOptions.fila.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.fila.filter );
            maximumResponse( &scale, &filaResponse, tempResponse, userOptions.scales[i] );
        }
        if ( userOptions.wall.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.wall.filter );
            maximumResponse( &scale, &wallResponse, tempResponse, userOptions.scales[i] );
        }
        
        // get the maximum eigenvalues and eigenvectors if asked so - this is correct since there can be only one environment computation when the user asks to output eigenvalues and/or eigenvectors - this is enforced when reading the program options
        if ( userOptions.eigenOn )
            eigenMaximumResponse( scale, &eigenvalues, &eigenvectors, *eValues, *eVectors, userOptions.scales[i] );
    }
    
    
    
    //! Output the results
    // free the dynamically allocated memory
    if ( userOptions.eigenOn )
    {
        delete eValues;
        delete eVectors;
    }
    
    outputResults( userOptions, mmfHeader, nodeResponse, filaResponse, wallResponse, scale, eigenvalues, eigenvectors );
}


/* The function used to output the results to file/files. */
void outputResults(User_options &userOptions,
                   MMF_header  &mmfHeader,
                   ArrayReal3D &nodeResponse,
                   ArrayReal3D &filaResponse,
                   ArrayReal3D &wallResponse,
                   Array<int,3>  &maxScale,
                   ArrayVector3D &eigenValues,
                   ArrayVector9D &eigenVectors)
{
    // output the maximum response
    // update some values in the MMF header
    mmfHeader.updateMaximumResponse();
    if ( userOptions.biasOn ) mmfHeader.updateBias( userOptions.bias );
    
    // write the maximum MMF response file
    if ( userOptions.node.envOn )
    {
        string mmfFileName = mmfHeader.maxResponseFilename( userOptions.node.outputFile );
        mmfHeader.updateFeature( userOptions.node.filter );
        if ( not userOptions.maximumScaleOn )
        {
            mmfHeader.fileType = MMF_MAX_RESPONSE;
            writeMMF_response( nodeResponse, mmfHeader, mmfFileName );
        }
    }
    if ( userOptions.fila.envOn )
    {
        string mmfFileName = mmfHeader.maxResponseFilename( userOptions.fila.outputFile );
        mmfHeader.updateFeature( userOptions.fila.filter );
        if ( not userOptions.maximumScaleOn )
        {
            mmfHeader.fileType = MMF_MAX_RESPONSE;
            writeMMF_response( filaResponse, mmfHeader, mmfFileName );
        }
    }
    if ( userOptions.wall.envOn )
    {
        string mmfFileName = mmfHeader.maxResponseFilename( userOptions.wall.outputFile );
        mmfHeader.updateFeature( userOptions.wall.filter );
        if ( not userOptions.maximumScaleOn )
        {
            mmfHeader.fileType = MMF_MAX_RESPONSE;
            writeMMF_response( wallResponse, mmfHeader, mmfFileName );
        }
    }
    
    
    // write the rest of the files, if the user asked so - these files are written for single environment computations
    string outputFile = "";
    ArrayReal3D *response;
    int feature = -1;
    if ( userOptions.node.envOn )
    {
        outputFile = userOptions.node.outputFile;
        response   = &nodeResponse;
        feature    = userOptions.node.feature;
    }
    if ( userOptions.fila.envOn )
    {
        outputFile = userOptions.fila.outputFile;
        response   = &filaResponse;
        feature    = userOptions.fila.feature;
    }
    if ( userOptions.wall.envOn )
    {
        outputFile = userOptions.wall.outputFile;
        response   = &wallResponse;
        feature    = userOptions.wall.feature;
    }
    
    // write the maximum scales file
    if ( userOptions.maximumScaleOn )
    {
        string mmfFileName = mmfHeader.maxResponseFilename( outputFile );
        mmfHeader.fileType = MMF_MAX_RESPONSE_SCALE;
        writeMMF_maxResponse( *response, maxScale, mmfHeader, mmfFileName );
    }    
    
    // write the eigenvalues and eigenvectors
    string mmfFileName = mmfHeader.maxEigenFilename( outputFile );
    if ( userOptions.allEigen )
    {
        mmfHeader.fileType = MMF_MAX_EIGEN;
        writeMMF_eigen( eigenValues, eigenVectors, mmfHeader, mmfFileName );
    }
    else if ( userOptions.singleEigen )
    {
        mmfHeader.fileType = MMF_MAX_EIGENVECTOR;
        featureEigenvector( eigenVectors, &eigenValues, feature );
        writeMMF( eigenValues, mmfHeader, mmfFileName, "feature eigenvector" );
    }
}


/* Returns the average density in the box. */
double averageDensity(ArrayReal3D &density)
{
    double res = 0.;
    for (Int i=0; i<density.size(); ++i)
        res += density[i];
    return res / density.size();
}
/* Rescales the density such that the average density is 'avgDen'. */
void rescaleDensity(ArrayReal3D &density, double avgDen)
{
    double res = 0.;
    for (Int i=0; i<density.size(); ++i)
        res += density[i];
    res /= density.size();
    Real factor = avgDen / res;
    for (Int i=0; i<density.size(); ++i)
        density[i] *= factor;
}




/* This function computes the MMF response using gaussian filtering in logarithmic space. It is given separately since involve a lot of steps different from the previous methods. */
void logarithmFilteringResponse(User_options &userOptions,
                                ArrayReal3D &density,
                                Density_header &densityHeader,
                                MMF_header &mmfHeader)
{
    // write some values in the output header
    mmfHeader.method = MMF_METHOD_LOGARITHM_FILTERING;
    Real boxLength[3];
    for (int i=0; i<3; ++i)
        boxLength[i] = mmfHeader.box[2*i+1] - mmfHeader.box[2*i];
    
    //! define the output quantities
    ArrayReal3D nodeResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for nodes
    if ( not userOptions.node.envOn ) nodeResponse.freeMemory();//release the memory if no node computations
    ArrayReal3D filaResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for filaments
    if ( not userOptions.fila.envOn ) filaResponse.freeMemory();//release the memory if no filament computations
    ArrayReal3D wallResponse( densityHeader.gridSize, NO_DIM ); // will keep track of the maximum response values for walls
    if ( not userOptions.wall.envOn ) wallResponse.freeMemory();//release the memory if no wall computations
    
    Array<int,3> scale( densityHeader.gridSize, NO_DIM );       // will keep track of the scale corresponding to the maximum response value
    scale.assign(userOptions.scales[0]);                    // initialize all the values to the current value of the scale
    ArrayVector3D eigenvalues( densityHeader.gridSize, NO_DIM );       // will store the eigenvalues corresponding to the maximum response
    ArrayVector9D eigenvectors( densityHeader.gridSize, NO_DIM );      // will store the eigenvectors corresponding to the maximum response (will store only the first NO_DIM-1 eigenvectors)
    
    // define additional variables to keep track of the response, eigenvalues and eigenvectors for each scale
    ArrayReal3D tempResponse( densityHeader.gridSize, NO_DIM );  // will store response for each scale != scale[0]
    ArrayVector3D *eValues;         // pointer to eigenvalues for each scale != scale[0]
    ArrayVector9D *eVectors;        // pointer to eigenvectors for each scale != scale[0]
    if ( userOptions.eigenOn )      // if the user asked for the eigenvalues and eigenvectors corresponding to the maximum response, need to allocate additional memory for that
    {
        eValues = new ArrayVector3D( densityHeader.gridSize, NO_DIM );
        eVectors = new ArrayVector9D( densityHeader.gridSize, NO_DIM );
    }
    else
    {
        eValues = &eigenvalues;
        eVectors = &eigenvectors;
    }
    
    
    
    //! loop over the scales - compute the MMF response and output the results
    // compute the quantities for the 1st scale
    cout << "\nComputations for scale: " << userOptions.scales[0] << "\n";
    ArrayReal3D tempDensity( densityHeader.gridSize, NO_DIM );
    double avgDen = averageDensity( density );
    
    // filter the density in logarithmic space
    densityLogarithm( &density, Real(1.e-3), false );
    gaussianFilter( density, userOptions.radia[0], boxLength, &tempDensity );
    densityExponential( &tempDensity, Real(1.), false );
    rescaleDensity( tempDensity, avgDen );
    
    // compute the hessian and hessian eigenvalues
    hessianMatrix( tempDensity, eigenvalues, eigenvectors, userOptions.radia[0], boxLength, userOptions.bias, userOptions.eigenOn );
    
    //compute the response for the given environments
    if ( userOptions.node.envOn )
        computeResponse( eigenvalues, nodeResponse, userOptions.node.filter );
    if ( userOptions.fila.envOn )
        computeResponse( eigenvalues, filaResponse, userOptions.fila.filter );
    if ( userOptions.wall.envOn )
        computeResponse( eigenvalues, wallResponse, userOptions.wall.filter );
    
    
    //compute the quantities for the rest of the scales
    for (size_t i=1; i<userOptions.scales.size(); ++i)
    {
        cout << "\nComputations for scale: " << userOptions.scales[i] << "\n";
        // filter the density in logarithmic space
        //densityLogarithm( density, Real(1.e-3), &tempDensity, false );
        gaussianFilter( density, userOptions.radia[i], boxLength, &tempDensity );
        densityExponential( &tempDensity, Real(1.), false );
        rescaleDensity( tempDensity, avgDen );
        
        // find the eignevalues and eigenvectors of the hessian matrix
        hessianMatrix( tempDensity, *eValues, *eVectors, userOptions.radia[i], boxLength, userOptions.bias, userOptions.eigenOn );
        
        // compute response for the given scale and the maximum response
        if ( userOptions.node.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.node.filter );
            maximumResponse( &scale, &nodeResponse, tempResponse, userOptions.scales[i] );
        }
        if ( userOptions.fila.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.fila.filter );
            maximumResponse( &scale, &filaResponse, tempResponse, userOptions.scales[i] );
        }
        if ( userOptions.wall.envOn )
        {
            computeResponse( *eValues, tempResponse, userOptions.wall.filter );
            maximumResponse( &scale, &wallResponse, tempResponse, userOptions.scales[i] );
        }
        
        // get the maximum eigenvalues and eigenvectors if asked so - this is correct since there can be only one environment computation when the user asks to output eigenvalues and/or eigenvectors - this is enforced when reading the program options
        if ( userOptions.eigenOn )
            eigenMaximumResponse( scale, &eigenvalues, &eigenvectors, *eValues, *eVectors, userOptions.scales[i] );
    }
    
    
    
    //! Output the results
    // free the dynamically allocated memory
    if ( userOptions.eigenOn )
    {
        delete eValues;
        delete eVectors;
    }
    
    outputResults( userOptions, mmfHeader, nodeResponse, filaResponse, wallResponse, scale, eigenvalues, eigenvectors );
}







/* Use the boost program options library to store the available user options.
The options are stored in two different structures:
    allOptions - stores all available options
    visibleOptions - stores the options visible to the user
*/
void addOptions(po::options_description &allOptions,
                po::options_description &visibleOptions,
                po::positional_options_description &p,
                User_options &userOptions)
{
    visibleOptions.add_options()
            ("help,h", "produce help message")
            ("node", po::value<string>(&(userOptions.node.outputFile)), "compute the response for nodes. Give the name of the output file for nodes.")
            ("fila", po::value<string>(&(userOptions.fila.outputFile)), "compute the response for filaments. Give the name of the output file for filaments.")
            ("wall", po::value<string>(&(userOptions.wall.outputFile)), "compute the response for walls. Give the name of the output file for walls.")
            
            ("filter,f", po::value<vector<int> >()->multitoken(), "specify the filter to be used to compute the response (e.g. '-f 41'). The entry format is:\n\t nodes = '4x' with x=0->3\n\t filaments = '3x' with x=0->3\n\t walls = '2x' with x=0->3. If you give this option you must supply a filter value for each Cosmic Web environment in the computation.")
            ("radius,r", po::value<float>(&(userOptions.radius0)), "filter radius for scale 0 in Mpc units ([DEFAULT] value 0.5 Mpc)")
            ("radia", po::value<vector<float> >(&(userOptions.radia))->multitoken(), "insert space sepparated real values giving the filter radius in Mpc for each scale used in the program (alternative way to enter 'radius'). It must have the same number of entries as the number of scales (e.g. for 3 scales '-radia 1.1 1.5 2.7').")
            ("base,b", po::value<float>(&(userOptions.base)), "filter radius at scale n is given by 'radius0 * base^n', where 'base' is the value inserted using this option and 'radius0' is the filter radius at scale 0. [DEFAULT] 'base' = sqrt(2).")
            ("bias", po::value<float>(&(userOptions.bias)), "insert a bias value - for bias>1 give larger weight to larger scales, while bias<1 gives larger weight to smaller scales. [DEFAULT] bias=1 (all scales have the same weight).")
            
            ("den", "compute the response using the density field. [DEFAULT]")
            ("logFilter", "compute the response using gaussian filtering in logarithm space.")
            ("log", "compute the response using the logarithm of the density field.")
            ("grav", "compute the response using the gravitational field (computed from the density distribution using the Poisson equation).")
            ("vel", "compute the MMF response using the divergence of the velocity field (input data must be a velocity divergence file).")
            ("velPot", "compute the MMF response using the velocity potential field (input data must be a velocity divergence file).")
            
            ("mask", "use a node mask for filament computation and use a node+filament mask for walls computations. With this option you also need to give node and/or filament clean files that will be used as masks - see the next two options. This option is recommended for the density and gravity methods.")
            ("nodeFile,N", po::value<string>(&(userOptions.nodeFile)), "specify the clean MMF response file for nodes used to substract the nodes from the image. This option is used for computing the response for filaments and walls.")
            ("filaFile,F", po::value<string>(&(userOptions.filaFile)), "specify the clean MMF response file for filaments used to substract the filaments from the image. This option is used for cleaning the response for walls.")
            
            ("allEigen", "choose this option if you would like to output a file with the eigenvalues and eigenvectors corresponding to the maximum response.")
            ("eigenvector", "choose this option if you would like to output a file with the eigenvector corresponding to the maximum response for filaments (eigenvector 3 which gives filament direction) or walls (eigenvector 1 which gives the normal to the wall plane).")
            ("box",po::value<vector<float> >(&(userOptions.boxCoordinates))->multitoken(), "insert the box coordinates (xMin, xMax, yMin, yMax, zMin and zMax) if not specified in the input data file.")
            ("maxScale", "output also the scale for the corresponding maximum MMF response value. This is an int value saved in the same file as the maximum MMF response, after the values of the response.")
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
             ("inputFile", po::value<string>(&(userOptions.inputFile)), "name of the input Gadget snapshot")
             ("scales", po::value<string>(&(userOptions.scale)), "the scales used to compute the maximum MMF response")
    ;
    
    allOptions.add(visibleOptions).add(hidden);
    
    // now add the hidden options to 'positional_options_description'
    p.add( "inputFile", 1 );
    p.add( "scales", 1 );
}


//print help information to the user
void helpInformation( po::options_description &visibleOptions, char *argv[] )
{
    cout << "Use this program to compute for a set of scales the maximum response function to cosmic web environment (node, filament or wall).\n";
    cout << "Usage:    " << argv[0] << "  name_density_file  scales(1:10,15,20) 'options - see below' \n";
    cout << "On top of the above, the user can specify the following additional options:\n";
    cout << visibleOptions << "\n";
    exit(0);
}


/* Print to the user what options the program will use. */
void printOptions( User_options &userOptions )
{
    cout << "RUNNING:   " << userOptions.programOptions << "\n\n";
    cout << "The program will compute the MMF response for a set of scales using the following input parameters:\n"
        << "\t input data file             : " << userOptions.inputFile << "\n"
        << "\t scales                      : " << userOptions.scale << "\n";
    if ( userOptions.node.envOn )
        cout << "\t node environment info       :   output file : " << userOptions.node.outputFile << "    and filter : " << userOptions.node.filter << "\n";
    if ( userOptions.fila.envOn )
        cout << "\t filament environment info   :   output file : " << userOptions.fila.outputFile << "    and filter : " << userOptions.fila.filter << "\n";
    if ( userOptions.wall.envOn )
        cout << "\t wall environment info       :   output file : " << userOptions.wall.outputFile << "    and filter : " << userOptions.wall.filter << "\n";
    if ( userOptions.radiaOn )
    {
        cout << "\t radia in Mpc                : ";
        for (size_t i=0; i<userOptions.radia.size(); ++i)
            cout << userOptions.radia[i] << "  ";
        cout << "\n";
    }
    else
        cout << "\t radius in Mpc               : " << userOptions.radius0 << "\n"
            << "\t base                        : " << userOptions.base << "\n";
    if ( userOptions.biasOn )
        cout << "\t bias                        : " << userOptions.bias << "\n";
    
    
    if ( not userOptions.mask )
        cout << "\t NO mask will be applied for the computation.\n";
    else if (userOptions.fila.envOn or userOptions.wall.envOn)
    {
        cout << "\t node mask file              : " << userOptions.nodeFile << "\n";
        if (userOptions.wall.envOn)
            cout << "\t filament mask file          : " << userOptions.nodeFile << "\n";
    }
    
    
    if ( userOptions.densityOn )
        cout << "\t computing the response using the density field.\n";
    else if ( userOptions.logarithmFiltering )
        cout << "\t computing the response using gaussian filtering in logarithmic space.\n";
    else if ( userOptions.logarithmOn )
        cout << "\t computing the response using the logarithm of the density field.\n";
    else if ( userOptions.gravitationOn )
        cout << "\t computing the response using the gravitational field.\n";
    else if ( userOptions.velocityOn )
        cout << "\t computing the response using the velocity divergence field.\n";
    else if ( userOptions.velocityPotentialOn )
        cout << "\t computing the response using the velocity potential field.\n";
    
    
    if ( userOptions.allEigen )
        cout << "\t outputting to files: 1= the maximum response and 2= the eigenvalues and eigenvectors corresponding to the maximum response.\n";
    else if ( userOptions.singleEigen )
        cout << "\t outputting to files: 1= the maximum response and 2= the eigenvector giving the filament direction or wall normal.\n";
    else
        cout << "\t outputting to file: 1= maximum response for over all scales.\n";
    
    if ( userOptions.boxCoordinatesOn )
        cout << "The data sits in a box of extension: [" << userOptions.boxCoordinates[0] << " to " << userOptions.boxCoordinates[1] << "], [" << userOptions.boxCoordinates[2] << " to " << userOptions.boxCoordinates[3] << "], [" << userOptions.boxCoordinates[4] << " to " << userOptions.boxCoordinates[5] << "] Mpc\n";
    
    if ( userOptions.maximumScaleOn )
        cout << "The program will output also the scales of the maximum response in the output file.\n";
    
    cout << "\nThe program will use the following scales to compute the response:\n" << "\tscale \tfilter radius (Mpc)\n";
    for (size_t i=0; i<userOptions.scales.size(); ++i)
        cout << "\t" << userOptions.scales[i] << "\t" << userOptions.radia[i] << "\n";
    cout << "\n";
}


/* Read the user supplied options and check that they satify some restrictions. */
void readOptions(int argc, char *argv[],
                 User_options *userOptions)
{
    po::options_description visibleOptions("Allowed options"), allOptions("All options");
    po::positional_options_description p;
    
    addOptions( allOptions, visibleOptions, p, *userOptions );  //read the options available to the program
    
    po::variables_map vm;
    po::store( po::command_line_parser(argc, argv).options(allOptions).positional(p).run(), vm );
    po::notify(vm);
    
    if ( not vm.count("help") and not vm.count("inputFile") ) cout << "~~~ERROR~~~ No input file detected.\n";
    if ( not vm.count("help") and not vm.count("scales") ) cout << "~~~ERROR~~~ No input scales detected.\n";
    
    if ( vm.count("help") or not vm.count("inputFile") or not vm.count("scales") )
        helpInformation( visibleOptions, argv ); // print available options and exit
    
    const char *optionsList[6] = { "den", "log", "logFilter", "grav", "vel", "velPot" };
    conflicting_options(vm, 6, optionsList);
    conflicting_options(vm, "allEigen", "eigenvector");
    
    
    
    int noEnvironments = 0;
    if ( vm.count("node") )
    {
        userOptions->node.envOn = true;
        ++noEnvironments;
    }
    if ( vm.count("fila") )
    {
        userOptions->fila.envOn = true;
        ++noEnvironments;
    }
    if ( vm.count("wall") )
    {
        userOptions->wall.envOn = true;
        ++noEnvironments;
    }
    if ( noEnvironments==0 )
    {
        cout << "~~~ERROR~~~ No cosmic web feature selected.\n";
        helpInformation( visibleOptions, argv );
    }
    
    if ( vm.count("mask") )
    {
        userOptions->mask = true;
        option_dependency( vm, "fila", "nodeFile" );
        option_dependency( vm, "wall", "nodeFile" );
        option_dependency( vm, "wall", "filaFile" );
        if ( noEnvironments>1 )
            throwError( "When using the '--mask' option you can compute the maximum response for up to an environment at a time. You selected to compute the response for several environments. Please select only one environment." );
    }
    
    int filterLimit[] = { 4,4,4};   //filter limits for nodes, filaments and walls
    if ( vm.count("filter") )
    {
        int const temp = vm["filter"].as< std::vector<int> >().size();
        if ( temp!=noEnvironments )
            throwError( "The number of arguments supplied to the '--filter' option must be equal to the number of environment's response that the program computes. You asked to do computations for ", noEnvironments, " Cosmic Web environments while you inserted ", temp, " values for the '--filter' option." );
        for (int i=0; i<temp; ++i)
        {
            int filter = vm["filter"].as< std::vector<int> >().at(i);
            int feature = filter / 10;
            int temp2 = filter % 10;
            
            // for node environments
            if ( feature==4 )
            { 
                if ( not userOptions->node.envOn )
                    throwError( "You supplied a filter value for the node environment but the node response will not be computed by the program." );
                else if ( temp2>filterLimit[0] and userOptions->node.envOn )
                    throwError( "The node environment filter can have values between 40 to ", userOptions->node.feature*10+filterLimit[0], ". You inserted a value of ", filter, " for the node filter." );
                else if (userOptions->node.envOn)
                    userOptions->node.filter = filter;
            }
            // for filament environments
            else if ( feature==3 )
            {
                if ( not userOptions->fila.envOn )
                    throwError( "You supplied a filter value for the filament environment but the filament response will not be computed by the program." );
                else if ( temp2>filterLimit[1] and userOptions->fila.envOn )
                    throwError( "The filament environment filter can have values between 30 to ", userOptions->fila.feature*10+filterLimit[1], ". You inserted a value of ", filter, " for the filament filter." );
                else if (userOptions->fila.envOn)
                    userOptions->fila.filter = filter;
            }
            // for wall environments
            else if ( feature==2 )
            {
                if ( not userOptions->wall.envOn )
                    throwError( "You supplied a filter value for the wall environment but the wall response will not be computed by the program." );
                else if ( temp2>filterLimit[2] and userOptions->wall.envOn )
                    throwError( "The wall environment filter can have values between 20 to ", userOptions->wall.feature*10+filterLimit[2], ". You inserted a value of ", filter, " for the wall filter." );
                else if (userOptions->wall.envOn)
                    userOptions->wall.filter = filter;
            }
            else
                throwError( "Unrecognized value of ", filter, "for the '--filter' option. Use '-h' to get the allowed values." );
        }
    }
    
    if ( vm.count("radia") ) userOptions->radiaOn = true;
    if ( vm.count("radius") and vm.count("radia") )
       throwError( "You need to insert the filter radius using only one of the available options ('--radius' to insert the radius for scale 0 or '--radia' to insert the radius for each scale), not more than one." );
    if ( vm.count("bias") ) userOptions->biasOn = true;
    
    
    if ( vm.count("den") ) userOptions->densityOn = true;
    if ( vm.count("logFilter") ) userOptions->logarithmFiltering = true;
    if ( vm.count("log") ) userOptions->logarithmOn = true;
    if ( vm.count("grav") ) userOptions->gravitationOn = true;
    if ( vm.count("vel") ) userOptions->velocityOn = true;
    if ( vm.count("velPot") ) userOptions->velocityPotentialOn = true;
    if ( not (userOptions->densityOn or userOptions->logarithmFiltering or userOptions->logarithmOn or userOptions->gravitationOn or userOptions->velocityOn or userOptions->velocityPotentialOn) )
        userOptions->densityOn = true;
    
    
    if ( vm.count("allEigen") and noEnvironments==1 )
        userOptions->allEigen = true;
    else if ( vm.count("allEigen") )
        throwError( "When using the '--allEigen' option you can compute the maximum response for up to an environment at a time. You selected to compute the response for several environments. Please select only one environment." );
    if ( vm.count("eigenvector") and noEnvironments==1 )
        userOptions->singleEigen = true;
    else if ( vm.count("eigenvector") )
        throwError( "When using the '--eigenvector' option you can compute the maximum response for up to an environment at a time. You selected to compute the response for several environments. Please select only one environment." );
    if ( userOptions->node.envOn ) userOptions->singleEigen = false;  // no eigenvector direction for clusters
    if ( userOptions->allEigen or userOptions->singleEigen ) userOptions->eigenOn = true;  // compute eigenvalues and eigenvectors corresponding to maximum response
    if ( vm.count("maxScale") and noEnvironments==1 )
        userOptions->maximumScaleOn = true;
    else if ( vm.count("maxScale") )
        throwError( "When using the '--maxScale' option you can compute the maximum response for up to an environment at a time. You selected to compute the response for several environments. Please select only one environment." );
    
    
    if ( vm.count("box") )
    {
        userOptions->boxCoordinatesOn = true;
        if ( userOptions->boxCoordinates.size()!=6 )
            throwError( "The program option '--box' must be followed by 6 values that give the box coordinates: xMin, xMax, yMin, yMax, zMin and zMax." );
        for (int i=0; i<3; ++i)
            if ( userOptions->boxCoordinates[2*i+1]<=userOptions->boxCoordinates[2*i] )
                throwError( "The values suplied to the program option '--box' must be such that: 'xMin<xMax', 'yMin<yMax' and 'zMin<zMax'." );
    }
    
    lowerBoundCheck( userOptions->base, float(0.), "value of the '--base' option" );
    
    // read the 'scales' supplied by the user
    getIntegerVector( userOptions->scale, &(userOptions->scales), '-', ',' );
    if ( duplicateEntries(&(userOptions->scales[0]), userOptions->scales.size()) )
        throwError( "There are one or more entries in the 'scale' array that have the same value. This is not permitted since the 2nd entry in the 'scale' array with the same value will overwrite the result from the first entry in the array." );
    // read/compute the values of the filter radia
    if ( not userOptions->radiaOn  )
    {
        userOptions->radia.clear();
        for (size_t i=0; i<userOptions->scales.size(); ++i)
            userOptions->radia.push_back( MMF_header::filterRadius( userOptions->radius0, userOptions->base, userOptions->scales[i] ) );
    }
    else if ( userOptions->scales.size()!=userOptions->radia.size() )
        throwError( "The number of filter radia values inserted using the option '--radia' does not match with the number of inserted scales." );
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        userOptions->programOptions += string( argv[i] ) + " ";
    
    
    printOptions( *userOptions );
}




