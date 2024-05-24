#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <omp.h>

#include <defines.h>
#include <box.h>
#include "MMF_file/MMF_file.h"
#include "MMF_file/MMF_header.h"
#include "miscellaneous/miscellaneous.h"
#include "MMF_contraction/MMF_contraction.h"
using namespace std;


#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;



struct User_options
{
    string inputFile;     // name of input density file
    string outputFile;    // root name of output MMF response files
    
    ContractionOptions options; // options for the contraction procedure - see the structure in file 'MMF_contraction.h' 
    string nodeFile;
    string filaFile;
    bool   nodeFileOn;
    bool   filaFileOn;
    
    
    bool   useInitialGuess; // true if to use an initial guess for the directions
    string directionFile;   // input a file giving the directions that will be used as initial guess for the computation
    bool   noIteration;     // true if the user did not request an iteration - only works with the '--initialGuess' option
    
    string programOptions;  // string that stores the 'argv' values for the users information
    
    
    User_options()
    {
        options.radius = 1.;
        options.periodic = true;
        options.convergeFraction    = 0.01;
        options.distanceThreshold   = 0.05;     // in Mpc/h
        options.eigenvalueThreshold = 0.3;      // ratio lambda1/lambda2 or lambda3/lambda2 for filaments and respectively walls
        options.maxLoop = 30;
        options.maxIterations = 20;
        options.cosDirectionThreshold = 0.996;  // less than 5 degrees
        
        nodeFileOn = false;
        filaFileOn = false;
        
        useInitialGuess = false;
        noIteration = false;
        programOptions = "";
    }
};
void readOptions(int argc, char *argv[],
                 User_options *userOptions);            // read the options inserted by the user







int main( int argc, char *argv[] )
{
    // read the user options
    User_options userOptions;
    readOptions( argc, argv, &userOptions );
    
    
    //! read the input data
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, userOptions.inputFile );
    if ( mmfHeader.fileType!=MMF_CLEAN_RESPONSE ) throwWarning( "The input file is not a single environment clean response file." );
    
    Array<shortInt,NO_DIM> response( mmfHeader.gridSize, NO_DIM );
    readMMF( &response, userOptions.inputFile );
    
    // read some values from the input file
    userOptions.options.box.assign( mmfHeader.box, 2*NO_DIM );
    userOptions.options.feature = mmfHeader.feature;
    
    
    // read the node and fila data and use it to add the masked regions
    if ( userOptions.nodeFileOn and userOptions.options.feature<4 )
    {
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.nodeFile );
        if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE ) throwWarning( "The input file is not a single environment clean response file." );
        
        Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
        readMMF( &tempResponse, userOptions.nodeFile );
        
        for (Int i=0; i<tempResponse.size(); ++i)
            if ( tempResponse[i]>0 ) response[i] = 1;
    }
    else if (userOptions.options.feature<4 )
        throwWarning( "No node clean response file was supplied to the program. Because of this the program needs to compute the directional information in the absence of the node information." );
    
    if ( userOptions.filaFileOn and userOptions.options.feature<3 )
    {
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.filaFile );
        if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE ) throwWarning( "The input file is not a single environment clean response file." );
        
        Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
        readMMF( &tempResponse, userOptions.filaFile );
        
        for (Int i=0; i<tempResponse.size(); ++i)
            if ( tempResponse[i]>0 ) response[i] = 1;
    }
    else if (userOptions.options.feature<3 )
        throwWarning( "No filament clean response file was supplied to the program. Because of this the program needs to compute the directional information in the absence of the filament information." );
    
    
    
    
    //! reserve memory for the output array
    size_t const noValidCells = noObjectCells( response );    // returns the number of valid filament/wall cells
    Array<Real,2> spinePosition( noValidCells, NO_DIM );
    Array<Real,2> direction( noValidCells, NO_DIM );
    Array<int,2>  cellPosition( noValidCells, NO_DIM );
    
    
    //! if the user specified an initial direction guess file, use it
    if ( userOptions.useInitialGuess )
    {
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.directionFile );
        if ( tempHeader.fileType!=MMF_MAX_EIGENVECTOR ) throwWarning( "The input file is not an environment eigenvector direction file." );
        
        Array<Real,NO_DIM+1> initialDirections( tempHeader.gridSize[0], tempHeader.gridSize[1],  tempHeader.gridSize[2], NO_DIM );
        readMMF( &initialDirections, userOptions.directionFile );
        
        // write the directions at valid pixels to the direction array
        Int count = 0;
        for (Int i=0; i<response.size(); ++i)
            if ( response[i]==1 )
            {
                for (int j=0; j<NO_DIM; ++j)
                    direction( count,j ) = initialDirections[ i*NO_DIM+j ];
                ++count;
            }
    }
    
    
    
    //! compute the geometrical direction of filaments/walls at each voxel
    if ( not userOptions.noIteration )
    {
        validObjectCell( response, cellPosition, spinePosition, userOptions.options.box );
        
        Real dx = (mmfHeader.box[1]-mmfHeader.box[0]) / mmfHeader.gridSize[0];
        userOptions.options.distanceThreshold = dx/5.;
        
        MMFObjectDirection( spinePosition, direction, userOptions.options, userOptions.useInitialGuess );
    }
    
    
    
    //! Output the results
    mmfHeader.fileType = MMF_DIRECTIONS;
    mmfHeader.gridSize[0] = direction.axisSize(0);
    mmfHeader.gridSize[1] = direction.axisSize(1);
    mmfHeader.gridSize[2] = 1; 
    mmfHeader.totalGrid = direction.size();
    mmfHeader.updateObservations( ";   ");
    mmfHeader.updateObservations( argv, argc );
    writeMMF( direction, mmfHeader, userOptions.outputFile, "filament/wall direction file");
//    return 0;
    
    
//    //! output in a text file the results
//    cout << "There are " << noValidCells << " points\n";
//    Array<Real,2> output( noValidCells, 3*NO_DIM );
//    int Nx = response.axisSize(0), Ny = response.axisSize(1), Nz = response.axisSize(2);
//    Box<Real,NO_DIM> box( userOptions.options.box );
//    Real dx[NO_DIM] = { (box[1]-box[0])/Nx, (box[3]-box[2])/Ny, (box[5]-box[4])/Nz };
//    int noCells = 0;
//    for (int i=0; i<noValidCells; ++i)
//        if ( spinePosition(i,0)+dx[0]/2>=60. and spinePosition(i,0)+dx[0]/2<=80. )
//        {
//            for (int j=0; j<NO_DIM; ++j)
//            {
//                output(noCells,0*NO_DIM+j) = box[2*j] + cellPosition(i,j)*dx[j]+dx[j]/2;
//                output(noCells,1*NO_DIM+j) = spinePosition(i,j)+dx[j]/2;
//                output(noCells,2*NO_DIM+j) = direction(i,j);
//            }
//            ++noCells;
//        }
//    fstream out;
//    openOutputTextFile( out, userOptions.outputFile );
//    cout << "There are " << noCells << " points\n";
//    for (int i=0; i<noCells; ++i)
//    {
//        out << output(i,0);
//        for (int j=1; j<output.axisSize(1); ++j)
//            out << "\t" << output(i,j);
//        out << "\n";
//    }
//    out.close();
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
            ("nodeFile",po::value<string>(&(userOptions.nodeFile)), "node clean response file")
            ("filaFile",po::value<string>(&(userOptions.filaFile)), "filament clean response file")
            ("radius",po::value<Real>(&(userOptions.options.radius)), "radius in Mpc used for the smoothing and neighbor search")
            ("initialGuess,i",po::value<std::string>(&(userOptions.directionFile)), "Give a direction file to be used as starting point for the iteration (optional).")
            ("convergenceFraction",po::value<Real>(&(userOptions.options.convergeFraction)), "stop when the change in convergence fraction is less than this value")
            ("eigenvalueThreshold",po::value<Real>(&(userOptions.options.eigenvalueThreshold)), "thresholds for eigenvalues for when convergence is achieved")
            ("distanceThreshold",po::value<Real>(&(userOptions.options.distanceThreshold)), "thresholds for displacement value for when convergence is achieved")
            ("directionThreshold",po::value<Real>(&(userOptions.options.cosDirectionThreshold)), "thresholds for when filament/wall direction has converged")
            ("maxLoop",po::value<int>(&(userOptions.options.maxLoop)), "number of maximum loops for each iteration")
            ("maxIteration",po::value<int>(&(userOptions.options.maxIterations)), "number of maximum iterations")
            ("nonPeriodic","specify this option if the data is non-periodic")
            ("noIteration","specify this option if the program should output the directions directly from input initial guess file - no iteration is performed")
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("inputFile", po::value<string>(&(userOptions.inputFile)), "name of the input MMF file")
            ("outputFile", po::value<string>(&(userOptions.outputFile)), "name of the output file")
    ;
    
    allOptions.add(visibleOptions).add(hidden);
    
    // now add the hidden options to 'positional_options_description'
    p.add( "inputFile", 1 );
    p.add( "outputFile", 1 );
}


//print help information to the user
void helpInformation( po::options_description &visibleOptions, char *argv[] )
{
    cout << "Use this program to compute the geometrical direction of filaments and walls.\n";
    cout << "Usage:    " << argv[0] << "  input_clean_file  name_output_file 'options - see below' \n";
    cout << "On top of the above, the user can specify the following additional options:\n";
    cout << visibleOptions << "\n";
    exit(0);
}


/* Print to the user what options the program will use. */
void printOptions( User_options &userOptions )
{
    cout << "RUNNING:   " << userOptions.programOptions << "\n\n";
    cout << "The program will compute the MMF response for a set of scales using the following input parameters:\n"
        << "\t input clean file            : " << userOptions.inputFile << "\n"
        << "\t output direction file       : " << userOptions.outputFile << "\n"
        << "\t clean node file (masking)   : " << (userOptions.nodeFileOn ? userOptions.nodeFile : "none given") << "\n"
        << "\t clean filament file(masking): " << (userOptions.filaFileOn ? userOptions.filaFile : "none given") << "\n";
    if ( userOptions.useInitialGuess )
        cout << "\t initial direction guess     : " << userOptions.directionFile << "\n";
    cout << "\t filter radius               : " << userOptions.options.radius << "  Mpc/h\n"
        << "\t convergence fraction        : " << userOptions.options.convergeFraction << "\n"
        << "\t distance threshold          : " << userOptions.options.distanceThreshold << "  Mpc/h\n"
        << "\t eigenvalue threshold        : " << userOptions.options.eigenvalueThreshold << "\n"
        << "\t direction cosinus threshold : " << userOptions.options.cosDirectionThreshold << "\n"
        << "\t boundary conditions         : " << (userOptions.options.periodic ? "periodic" : "free (non-periodic)") << "\n"
        << "\t maximum loops               : " << userOptions.options.maxLoop << "\n"
        << "\t maximum iterations          : " << userOptions.options.maxIterations << "\n";
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
    if ( not vm.count("help") and not vm.count("outputFile") ) cout << "~~~ERROR~~~ No output file detected.\n";
    
    if ( vm.count("help") or not vm.count("inputFile") or not vm.count("outputFile") )
        helpInformation( visibleOptions, argv ); // print available options and exit
    
    if ( vm.count("nodeFile") ) userOptions->nodeFileOn = true;
    if ( vm.count("filaFile") ) userOptions->filaFileOn = true;
    if ( vm.count("nonPeriodic") ) userOptions->options.periodic = false;
    if ( vm.count("initialGuess") ) userOptions->useInitialGuess = true;
    if ( vm.count("noIteration") ) userOptions->noIteration = true;
    if ( vm.count("noIteration") and not vm.count("initialGuess") )
        throwError( "The option '--noIteration' can only be used only with the '--initialGuess' option." );
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        userOptions->programOptions += string( argv[i] ) + " ";
    
    
    printOptions( *userOptions );
}




