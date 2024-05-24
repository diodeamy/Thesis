#include <iostream>
#include <cmath>

#include <defines.h>
#include "densityFile/densityFile.h"
#include <density_header.h>
#include <array.h>
#include "MMF_file/MMF_file.h"
#include <MMF_header.h>
#include "miscellaneous/miscellaneous.h"
#include "MMF_computations/MMF_computations.h"
#include "MMF_computations/objects.h"
#include "classes/FFTW.h"

using namespace std;



#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct User_options
{
    string inputFile;     // name of input density file
    string outputFile;    // root name of output threshold behavior files
    string densityFile;   // name of the input density file for nodes
    string nodeFile;      // the clean MMF response of nodes used for filament and wall computations
    string filaFile;      // the clean MMF response of filaments used for wall computations
    
    int    feature;       // cosmic web feature (4=blob, 3=filament, 2=wall)
    float  threshold;     // threshold values used to get the clean MMF response
    
    float  minimumSize;   // the minimum size of the objects that will be counted
    int    neighborType;  // can be 1 or 2 - described what is a neighbor 1= take only the 6 closest neighbors, 2= take only the 26th closest neighbors
    bool   mask;          // true if to use a mask for the fillament and wall computations
    
    string programOptions;// string that stores the 'argv' values for the users information
    Real rhoCritical;     // the critical density in units of (M0/h) / (Mpc/h)^3
    vector<float> minSizeDefault; //default values for the 'minimumSize' option
    
    
    
    User_options()
    { feature=-1; threshold = -1.;  minimumSize = 0.; neighborType = 1;
    mask = true;
    rhoCritical = 27.7538e+10;
    minSizeDefault.push_back(10.); minSizeDefault.push_back(10.); minSizeDefault.push_back(5.e13);}
    
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
    
    
    //! Compute the clean MMF response
    Array<shortInt,NO_DIM> cleanData( mmfHeader.gridSize, NO_DIM );
    double minSizeGridUnits = userOptions.MinimumSize( mmfHeader.gridCellVolume(), mmfHeader.Omega0 );
    if ( userOptions.feature==4 )   // for nodes
    {
        // read the density file
        cout << "\nReading the density data for node identification:" << flush;
        Density_header densityHeader;
        readDensityHeader( &densityHeader, userOptions.densityFile );
        ArrayReal3D density( densityHeader.gridSize, NO_DIM );
        readDensityData( &density, userOptions.densityFile );
        
        // compute the clean response
        cleanResponse( response, userOptions.threshold, userOptions.neighborType, density, minSizeGridUnits, &cleanData, true );
    }
    else
        cleanResponse( response, userOptions.threshold, userOptions.neighborType, int(minSizeGridUnits+.99), &cleanData, true );
    
    mmfHeader.fileType = MMF_CLEAN_RESPONSE;
    mmfHeader.updateObservations( ";   ");
    mmfHeader.updateObservations( argv, argc );
    writeMMF( cleanData, mmfHeader, userOptions.outputFile, "clean MMF response");
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
            ("node,4", "clean the MMF response for nodes.")
            ("fila,3", "clean the MMF response for filaments.")
            ("wall,2", "clean the MMF response for walls.")
            ("densityFile,D", po::value<string>(&(userOptions.densityFile)), "specify the density binary file from which the MMF response was computed - used for optimal node threshold detection.")
            ("nodeFile,N", po::value<string>(&(userOptions.nodeFile)), "specify the clean MMF response file for blobs used to substract the blobs from the image. This option is used for cleaning the response for filaments and walls.")
            ("filaFile,F", po::value<string>(&(userOptions.filaFile)), "specify the clean MMF response file for filaments used to substract the filaments from the image. This option is used for cleaning the response for walls.")
            ("minSize", po::value<float>(&(userOptions.minimumSize)), "give the minimum size of what means a significant object. This should be a mass value for nodes (in units of Msolar/h) and a volume (in units of (Mpc/h)^3) for filaments and walls. Default values: 5.e13 Msolar/h for nodes and 10 (Mpc/h)^3 for filaments/walls.")
            ("neighbor,n", po::value<int>(&(userOptions.neighborType))->default_value(userOptions.neighborType), "specify what constitutes a neighbor for a pixel/voxel - all neighbors are grouped as the same object. There are two different types of neigbors:\n\t 1 - take as neighbors only the closest 6 points in 3D (4 in 2D) as voxel neighbors[DEFAULT].\n\t 2 - take as neighbors only the closest 26 points in 3D (8 in 2D) as voxel neighbors.")
            ("noMask", "do not use a node mask for filament computation and do not use a node+filament mask for walls computations.")
            ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("inputFile", po::value<string>(&(userOptions.inputFile)), "name of the input MMF maximum response file")
            ("outputFile", po::value<string>(&(userOptions.outputFile)), "root name of the output MMF file")
            ("threshold", po::value<float>(&(userOptions.threshold)), "threhold value for the MMF response")
            ;
    
    allOptions.add(visibleOptions).add(hidden);
    
    // now add the hidden options to 'positional_options_description'
    p.add( "inputFile", 1 );
    p.add( "outputFile", 1 );
    p.add( "threshold", 1 );
}


//print help information to the user
void helpInformation( po::options_description &visibleOptions, char *argv[] )
{
    cout << "Use this program to compute the clean MMF reponse for the cosmic web features (nodes, filaments or walls). The program outputs a 'clean.MMF' file which contains 'short int' values for each grid point with value=1 if the respective grid cell is part of the cosmic web feature of interest.\n";
    cout << "Usage:    " << argv[0] << "  name_MMF_max_response_file  root_name_output_file  threshold_value  'options - see below' \n";
    cout << "On top of the above, the user can specify the following additional options (many options are available only for a given feature of the cosmic web):\n";
    cout << visibleOptions << "\n";
    exit(0);
}


/* Print to the user what options the program will use. */
void printOptions( User_options &userOptions )
{
    string feature = "Unknown";
    if ( userOptions.feature==4 ) feature = "nodes";
    else if ( userOptions.feature==3 ) feature = "filaments";
    else if ( userOptions.feature==2 ) feature = "walls";
    
    cout << "RUNNING:   " << userOptions.programOptions << "\n\n";
    cout << "The program will compute the clean MMF response for the following input parameters:\n"
        << "\t maximum MMF response file   : " << userOptions.inputFile << "\n"
        << "\t root name output file       : " << userOptions.outputFile << "\n"
        << "\t feature                     : " << feature << "\n";
    if ( userOptions.minimumSize>0. )
        cout << "\t minimum object size         : " << userOptions.minimumSize << (userOptions.feature==4 ? " Msolar/h" : " (Mpc/h)^3") << "\n";
    cout << "\t number of cell neighbors    : " << ((userOptions.neighborType==1)?6:26) << "\n";
    if ( not userOptions.mask and (userOptions.feature==3 or userOptions.feature==2) )
        cout << "\t NO mask will be applied for the computation.\n";
    else if (userOptions.feature==3 or userOptions.feature==2)
    {
        cout << "\t node mask MMF file          : " << userOptions.nodeFile << "\n";
        if (userOptions.feature==2)
            cout << "\t filament mask MMF file      : " << userOptions.nodeFile << "\n";
    }
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
    if ( not vm.count("help") and not vm.count("threshold") ) cout << "~~~ERROR~~~ No threshold value detected.\n";
    
    if ( vm.count("help") or not vm.count("inputFile") or not vm.count("outputFile") or not vm.count("threshold") )
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
    
    
    if ( vm.count("node") ) userOptions->feature = 4;
    if ( vm.count("fila") ) userOptions->feature = 3;
    if ( vm.count("wall") ) userOptions->feature = 2;
    if ( userOptions->feature<1 )
        throwError( "No cosmic web feature selected. Use one of the '--node','--fila' or '--wall' to select a Cosmic Web environment." );
    if ( vm.count("node") and not vm.count("densityFile") )
        throwError( "Computing the clean response for node environments needs as input a density file. Use the option '--densityFile' to give the name of the density file used to compute the environmental response." );
    
    if ( not vm.count("minSize") ) userOptions->minimumSize = userOptions->minSizeDefault[userOptions->feature-2];
    if ( vm.count("noMask") ) userOptions->mask = false;
    intervalCheck( userOptions->neighborType, 1, 2, "program option 'neighbor'" );
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        userOptions->programOptions += string( argv[i] ) + " ";
    
    
    printOptions( *userOptions );
}





