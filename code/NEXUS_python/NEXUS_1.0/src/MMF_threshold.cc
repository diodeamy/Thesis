#include <iostream>
#include <cmath>

#include "densityFile/densityFile.h"
#include <density_header.h>
#include "MMF_file/MMF_file.h"
#include <MMF_header.h>
#include "miscellaneous/miscellaneous.h"
#include "MMF_computations/MMF_computations.h"
#include "MMF_computations/objects.h"
#include "classes/FFTW.h"
#include "output/MMF_thresholdOutput.h"
#include <array.h>
#include <linear_fit.h>
#include <arrayProperties.h>

using namespace std;
typedef float   Real;



#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct User_options
{
    string inputFile;     // name of input density file
    string outputFile;    // root name of output threshold behavior files
    string densityFile;   // name of the input density file for nodes
    string nodeFile;      // the MMF max response of nodes used for filament and wall computations
    string filaFile;      // the MMF max response of filaments used for wall computations
    
    int    feature;       // cosmic web feature (4=blob, 3=filament, 2=wall)
    vector<float> thresholds; // the range of thresholds for which to compute the environment behavior on the response threshold
    bool   rangeOn;       // true if the user supplied a range to the program
    vector<float> range;  // store the range details inserted by the user
    
    float  minimumSize;   // the minimum size of the objects that will be counted
    int    neighborType;  // can be 1 or 2 - described what is a neighbor 1= take only the 6 closest neighbors, 2= take only the 26th closest neighbors
    float  virialDensity; // specify what is the minimal average density of a blob to be considered a real detection - if not supplied, the code computes the virial density
    bool   mask;          // true if to use a mask for the filament and wall computations
    bool   cleanOutputOn; // true if to output the clean file too
    string cleanOutputFile; // the name of the clean file where to output the clean response
    
    string programOptions;// string that stores the 'argv' values for the users information
    
    
    //default values for the range parameter
    vector<float> rangeDefault;
    vector<float> minSizeDefault;
    Real rhoCritical;    //the critical density in units of (M0/h) / (Mpc/h)^3
    
    
    User_options()
    { feature=-1; rangeOn = false; minimumSize = 0.; neighborType = 1; virialDensity = 350.;
    mask = true; cleanOutputOn = false;
    rangeDefault.push_back(9.e-1); rangeDefault.push_back(1.e-5); rangeDefault.push_back(50.);
    minSizeDefault.push_back(10.); minSizeDefault.push_back(10.); minSizeDefault.push_back(5.e13);
    rhoCritical = 27.7538e+10;}
    
    double MinimumSize(Real gridCellVolume, Real Omega0)
    {
        if (feature==4)
        {
            Real cellMass = rhoCritical * gridCellVolume * Omega0;
            minimumSize /= cellMass;
            //cout << "\n Minimum node mass = " << minimumSize << " in average mass units inside a given voxel.\n";
        }
        else minimumSize /= gridCellVolume;
        return minimumSize;
    }
};
void readOptions(int argc, char *argv[],
                 User_options *userOptions);            // read the options inserted by the user








int main( int argc, char *argv[] )
{
    // read the user options
    User_options userOptions;
    readOptions( argc, argv, &userOptions );
    
    
    //! read from the MMF file the maximum response
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, userOptions.inputFile );
    if ( mmfHeader.fileType!=MMF_MAX_RESPONSE ) throwWarning( "The input file is not a maximum MMF response file." );
    
    ArrayReal3D response( mmfHeader.gridSize, NO_DIM );
    readMMF_response( &response, userOptions.inputFile );
    double minSize = userOptions.MinimumSize( mmfHeader.gridCellVolume(), mmfHeader.Omega0 );    // minimum size as a number of grid cell units
    
    
    
    //! Read the masking files, if any, and apply the masks
    if ( userOptions.mask and (userOptions.feature==3 or userOptions.feature==2) ) //apply node mask to filament & wall computations
    {
        // open node MMF response
        cout << "\nReading the node clean MMF response:" << flush;
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.nodeFile );
        mmfHeader.compatible( tempHeader );// check that the two MMF files are compatible (i.e. have the same grid dimensions)
        if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE )
            throwError( "The input file '" + userOptions.nodeFile + "' is not a 'clean' MMF response file." );

        Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
        readMMF_cleanResponse( &tempResponse, userOptions.nodeFile );
        
        //discard the values of 'response' in the pixels where 'tempResponse' a valid mask
        cout << "Masking the maximum MMF response using the node clean MMF response ... " << flush;
        maskResponse( &response, tempResponse, shortInt(1) );
        cout << "Done.\n";
    }
    if ( userOptions.mask and userOptions.feature==2 ) //apply filament mask to wall computations
    {
        // open node MMF response
        cout << "\nReading the filament clean MMF response:" << flush;
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.filaFile );
        mmfHeader.compatible( tempHeader );// check that the two MMF files are compatible (i.e. have the same grid dimensions)
        if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE )
            throwError( "The input file '" + userOptions.filaFile + "' is not a 'clean' MMF response file." );

        Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
        readMMF_cleanResponse( &tempResponse, userOptions.filaFile );
        
        //discard the values of 'response' in the pixels where 'tempResponse' a valid mask
        cout << "Masking the maximum MMF response using the filament clean MMF response ... " << flush;
        maskResponse( &response, tempResponse, shortInt(1) );
        cout << "Done.\n";
    }
    
    
    //! Read the density file
    cout << "\nReading the density data:" << flush;
    Density_header densityHeader;
    readDensityHeader( &densityHeader, userOptions.densityFile );
    
    ArrayReal3D density( densityHeader.gridSize, NO_DIM );
    readDensityData( &density, userOptions.densityFile );
    
    
    
    //! Call the individual functions that compute the optimal threshold for each environment
    // define array to keep track of threshold and threshold response
    Real maximumResponseValue = Max( response.ptrData(), response.size() );
    size_t const noThresholds = userOptions.thresholds.size();
    Array<Real,1> threshold( noThresholds );
    for (size_t i=0; i<noThresholds; ++i)
        threshold[i] = userOptions.thresholds[i] * maximumResponseValue;
    
    Real optimalThreshold = 0.;
    switch ( userOptions.feature )
    {
        case 4:     //nodes
        {
            optimalThreshold =  outputNodeThreshold( density, response, userOptions.neighborType, userOptions.virialDensity, minSize,  threshold, userOptions.outputFile, userOptions.programOptions ); //writes the results to an ASCII file
            break;
        }
        case 3:     // filaments
        {
            optimalThreshold = outputWallThreshold( density, response, userOptions.neighborType, int(minSize+.99), threshold, userOptions.outputFile, userOptions.programOptions ); //writes the results to an ASCII file
            break;
        }
        case 2:     //walls
        {
            optimalThreshold = outputWallThreshold( density, response, userOptions.neighborType, int(minSize+.99), threshold, userOptions.outputFile, userOptions.programOptions ); //writes the results to an ASCII file
            break;
        }
    }
    
    
    //! Output the clean file if the user requested so
    if ( userOptions.cleanOutputOn )
    {
        Array<shortInt,NO_DIM> cleanData( mmfHeader.gridSize, NO_DIM );
        if ( userOptions.feature==4 )
            cleanResponse( response, optimalThreshold, userOptions.neighborType, density, minSize, &cleanData, true );
        else
            cleanResponse( response, optimalThreshold, userOptions.neighborType, int(minSize+.99), &cleanData, true );
        
        mmfHeader.fileType = MMF_CLEAN_RESPONSE;
        mmfHeader.updateObservations( ";   ");
        mmfHeader.updateObservations( argv, argc );
        writeMMF( cleanData, mmfHeader, userOptions.cleanOutputFile, "clean MMF response");
    }
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
            ("node,4", "compute the MMF response for nodes.")
            ("fila,3", "compute the MMF response for filaments.")
            ("wall,2", "compute the MMF response for walls.")
            ("densityFile,D", po::value<string>(&(userOptions.densityFile)), "specify the density binary file from which the MMF response was computed - used for optimal node threshold detection.")
            ("nodeFile,N", po::value<string>(&(userOptions.nodeFile)), "specify the clean MMF response file for nodes used to substract the nodes from the image. This option is used for optimal filament and wall threshold detection.")
            ("filaFile,F", po::value<string>(&(userOptions.filaFile)), "specify the clean MMF response file for filaments used to substract the filaments from the image. This option is used for optimal wall threshold detection.")
            ("range,r", po::value<vector<float> >(&(userOptions.range))->multitoken(), "insert space separated real values giving the minimum and maximum values of the threshold parameter as well as the number of threshold points. The maximum and minimum must be given as a fraction with respect to the maximum MMF response (i.e. --range '1.e-3 1.e-1 20' to compute the behavior for thresholds from 1.e-3 to 1.e-1 of the maximum response value).")
            ("minSize", po::value<float>(&(userOptions.minimumSize)), "give the minimum size of what means a significant object. This should be a mass value for nodes (in units of Msolar/h) and a volume (in units of (Mpc/h)^3) for filaments and walls. Default values: 5.e13 Msolar/h for nodes and 10 (Mpc/h)^3 for filaments/walls.")
            ("neighbor", po::value<int>(&(userOptions.neighborType))->default_value(userOptions.neighborType), "specify what constitutes a neighbor for a pixel/voxel - all neighbors are grouped as the same object. There are two different types of neigbors:\n\t 1 - take as neighbors only the closest 6 points in 3D (4 in 2D) as voxel neighbors.\n\t 2 - take as neighbors only the closest 26 points in 3D (8 in 2D) as voxel neighbors.")
            ("virDen", po::value<float>(&(userOptions.virialDensity)), "specify the virial density used for optimal threshold detection for nodes. Feature used only for the node threshold computation.")
            ("noMask", "do not use a node mask for filament threshold computation and do not use a node+filament mask for walls threshold computations.")
            ("cleanFile", po::value<string>(&(userOptions.cleanOutputFile)), "Specify the name of the output clean file for the current environments. If given, the code will output the clean response for the optimal threshold.")
            ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("inputFile", po::value<string>(&(userOptions.inputFile)), "name of the input maximum MMF response file")
            ("outputFile", po::value<string>(&(userOptions.outputFile)), "root name of the output ASCII file/files")
            ;
    
    allOptions.add(visibleOptions).add(hidden);
    
    // now add the hidden options to 'positional_options_description'
    p.add( "inputFile", 1 );
    p.add( "outputFile", 1 );
}


//print help information to the user
void helpInformation( po::options_description &visibleOptions, char *argv[] )
{
    cout << "Use this program to compute the optimal threshold value for different cosmic web environments (node, filament or wall). The program will output a '.count' file which for nodes gives the number of objects larger than the given size as a function of threshold while for filaments/walls it outputs the volume fraction taken by the largest filament/wall.\n";
    cout << "Usage:    " << argv[0] << "  name_MMF_max_response_file  root_name_output_ASCII_file  'options - see below' \n";
    cout << "On top of the above, the user can specify the following additional options (many options are available only for a given feature of the cosmic web):\n";
    cout << visibleOptions << "\n";
    exit(0);
}


/* Print to the user what options the program will use. */
void printOptions( User_options &userOptions )
{
    cout << "RUNNING:   " << userOptions.programOptions << "\n\n";
    cout << "The program will compute the optimal threshold value using the following input parameters:\n"
        << "\t maximum MMF response file   : " << userOptions.inputFile << "\n"
        << "\t root name output ASCII file : " << userOptions.outputFile << "\n"
        << "\t density file                : " << userOptions.densityFile << "\n"
        << "\t feature                     : ";
    if ( userOptions.feature==4 )
    {
        cout << "node\n"
        << "\t virial density              : ";
        if ( userOptions.virialDensity>0. ) cout <<  userOptions.virialDensity << "\n";
        else cout << "value computed later using analytical formula" << "\n";
    }
    else if ( userOptions.feature==3 )
    {
        cout << "filament\n";
        if ( userOptions.mask )
            cout << "\t node clean MMF response file: " << userOptions.nodeFile << "\n";
    }
    else if ( userOptions.feature==2 )
    {
        cout << "wall\n";
        if ( userOptions.mask )
            cout << "\t node clean MMF response file: " << userOptions.nodeFile << "\n"
                << "\t filament clean response file: " << userOptions.filaFile << "\n";
    }
    if ( userOptions.minimumSize>0. )
        cout << "\t minimum object size         : " << userOptions.minimumSize << (userOptions.feature==4 ? " Msolar/h" : " (Mpc/h)^3") << "\n";
    cout << "\t number of cell neighbors    : " << ((userOptions.neighborType==1)?6:26) << "\n";
    cout << "\t There will be " << userOptions.range[2] << " threshold points in the interval " << userOptions.range[1] << " to " << userOptions.range[0] << " of the maximum response value.\n";
    
    if ( userOptions.cleanOutputOn )
        cout << "\t the clean response will be writen to the file '" << userOptions.cleanOutputFile << "'.\n";
    else
        cout << "\t the clean response will NOT be computed.\n";
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
    if ( not vm.count("help") and not vm.count("outputFile") ) cout << "~~~ERROR~~~ No output file detected.\n";
    
    if ( vm.count("help") or not vm.count("inputFile") or not vm.count("outputFile") )
        helpInformation( visibleOptions, argv ); // print available options and exit
    
    conflicting_options( vm, "node", "fila" );
    conflicting_options( vm, "node", "wall" );
    conflicting_options( vm, "fila", "wall" );
    
    if ( not vm.count("noMask") )
    {
        option_dependency( vm, "fila", "nodeFile" );
        option_dependency( vm, "wall", "nodeFile" );
        option_dependency( vm, "wall", "filaFile" );
    }
    
    superfluous_options( vm, "virDen", "node" );
    superfluous_options( vm, "nodeFile", vm.count("fila") or vm.count("wall"), "'filament' or 'wall'" );
    superfluous_options( vm, "filaFile", "wall" );
    
    
    
    if ( vm.count("node") ) userOptions->feature = 4;
    if ( vm.count("fila") ) userOptions->feature = 3;
    if ( vm.count("wall") ) userOptions->feature = 2;
    if ( userOptions->feature<1 )
        throwError( "No cosmic web feature selected. Use one of the '--node','--fila' or '--wall' to select a Cosmic Web environment." );
    if ( not vm.count("densityFile") )
        throwError( "No density file supplied to the program. Use the option '--densityFile' to give the name of the density file used to compute the environmental response." );
    if ( not vm.count("minSize") ) userOptions->minimumSize = userOptions->minSizeDefault[userOptions->feature-2];
    
    if ( vm.count("noMask") ) userOptions->mask = false;
    if ( vm.count("cleanFile") ) userOptions->cleanOutputOn = true;
    if ( vm.count("range") )
    {
        userOptions->rangeOn = true;
        intervalCheck( userOptions->range[0], float(0.), float(1.), "first argument of program option 'range'" );
        intervalCheck( userOptions->range[1], float(0.), float(1.), "second argument of program option 'range'" );
        lowerBoundCheck( userOptions->range[2], float(5.), "third argument of program option 'range'" );
        float rMin = Min<float>( userOptions->range[0], userOptions->range[1] );
        float rMax = Max<float>( userOptions->range[0], userOptions->range[1] );
        userOptions->range[0] = rMax;
        userOptions->range[1] = rMin;
    }
    else userOptions->range = userOptions->rangeDefault;
    
    intervalCheck( userOptions->neighborType, 1, 2, "program option 'neighbor'" );
    
    
    // compute the values for the threshold
    size_t const noThresholds = size_t( userOptions->range[2] );
    userOptions->thresholds.reserve( noThresholds );
    userOptions->thresholds.push_back( userOptions->range[0] );
    Real step = log10( userOptions->range[1] / userOptions->range[0] ) / (noThresholds-1);
    for (size_t i=1; i<noThresholds; ++i)
        userOptions->thresholds.push_back( userOptions->range[0] * pow( Real(10.), Real(step*i) ) );
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        userOptions->programOptions += string( argv[i] ) + " ";
    
    
    printOptions( *userOptions );
}





