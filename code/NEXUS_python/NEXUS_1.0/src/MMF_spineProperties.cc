#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <omp.h>

#include <defines.h>
#include <box.h>
#include "densityFile/densityFile.h"
#include "densityFile/density_header.h"
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
    string inputFile;     // name of input spine direction file
    string outputFile;    // name of output text file
    
    ContractionOptions options; // options for the contraction procedure - see the structure in file 'MMF_contraction.h' 
    string nodeFile;
    string filaFile;
    bool   nodeFileOn;
    bool   filaFileOn;
    
    string densityFile;   // the file that has the density data
    string cleanFile;     // the clean file corresponding to the input spine direction file
    
    string programOptions;// string that stores the 'argv' values for the users information
    
    
    User_options()
    {
        options.radius = 1.;
        options.periodic = true;
        options.convergeFraction    = 0.005;
        options.distanceThreshold   = 0.05;     // in Mpc/h
        options.eigenvalueThreshold = 0.3;      // ratio lambda1/lambda2 or lambda3/lambda2 for filaments and respectively walls
        options.maxLoop = 60;
        
        nodeFileOn = false;
        filaFileOn = false;
        
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
    // read the spine direction data
    cout << "Reading the spine direction data:\n";
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, userOptions.inputFile );
    if ( mmfHeader.fileType!=MMF_DIRECTIONS ) throwWarning( "The input file is not spine direction binary file." );
    
    Array<Real,2> direction( mmfHeader.gridSize, 2 );
    readMMF( &direction, userOptions.inputFile );
    // read some values from the input file
    userOptions.options.addBoxCoordinates( mmfHeader.box );
    userOptions.options.feature = mmfHeader.feature;
    
    
    // read the clean response data
    cout << "Reading the clean response data:\n";
    MMF_header mmfHeader2;
    readMMF_Header( &mmfHeader2, userOptions.cleanFile );
    if ( mmfHeader2.fileType!=MMF_CLEAN_RESPONSE ) throwWarning( "The input file is not a single environment clean response file." );
    
    Array<shortInt,NO_DIM> response( mmfHeader2.gridSize, NO_DIM );
    readMMF( &response, userOptions.cleanFile );
    for (int i=0; i<NO_DIM; ++i) userOptions.options.grid[i] = mmfHeader2.gridSize[i];
    
    
    // read the density data
    cout << "Reading the density field data:\n";
    Density_header densityHeader;
    readDensityHeader( &densityHeader, userOptions.densityFile );
    
    ArrayReal3D density( densityHeader.gridSize, NO_DIM );
    readDensityData( &density, userOptions.densityFile );
    mmfHeader2.compatible( densityHeader.gridSize );    //check if both files have the same grid size
    
    
    // read the node and fila data and use it to add the masked regions
    if ( userOptions.nodeFileOn and userOptions.options.feature<4 )
    {
        MMF_header tempHeader;
        readMMF_Header( &tempHeader, userOptions.nodeFile );
        if ( tempHeader.fileType!=MMF_CLEAN_RESPONSE ) throwWarning( "The input file is not a single environment clean response file." );
        
        Array<shortInt,NO_DIM> tempResponse( tempHeader.gridSize, NO_DIM );
        readMMF( &tempResponse, userOptions.nodeFile );
        
        for (Int i=0; i<tempResponse.size(); ++i)
            if ( tempResponse[i]>0 )
            {
                response[i] = 1;
//                density[i] = 0.;
            }
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
            if ( tempResponse[i]>0 )
            {
                response[i] = 1;
//                density[i] = 0.;
            }
    }
    else if (userOptions.options.feature<3 )
        throwWarning( "No filament clean response file was supplied to the program. Because of this the program needs to compute the directional information in the absence of the filament information." );
    
    
    
    
    //! get the density values corresponding to valid objects
    Array<Real,1> mass( mmfHeader.gridSize[0] );
    size_t count = 0;
    Real averageMass = densityHeader.massInCell();
    for (size_t i=0; i<response.size(); ++i)
        if ( response[i]>shortInt(0) )
            mass[count++] = density[i] * averageMass;
    if ( count!=mmfHeader.gridSize[0] )
        throwError( "There was an incompatibility between the spine direction file and the clean response file. The valid voxels in the two files are not the same." );
    
    
    
    //! reserve memory for the output array
    size_t const noBins = 70;
    Real const lengthBoundaries[] = {0., 10.};
    Real const massBoundaries[] = {1.e8, 1.e15};
    Real const dL = (lengthBoundaries[1]-lengthBoundaries[0])/noBins, dM = log10(massBoundaries[1]/massBoundaries[0])/noBins;
    Array<double,2> sizeData(noBins+1,5);
    sizeData.assign( 0. );
    Array<double,2> massData(noBins+1,3);
    massData.assign( 0. );
    Array<Real,1> averageData(2);
    for (size_t i=0; i<noBins+1; ++i)
    {
        sizeData(i,0) = lengthBoundaries[0] + i*dL;
        massData(i,0) = massBoundaries[0] * pow(10.,i*dM);
    }
    
    
    
    
    //! contract the filaments/walls to the main spine/plane and compute the properties of each object
    cout << "\n";
    Array<Real,2> spinePosition( mmfHeader.gridSize, 2 );
    validObjectCell( response, spinePosition, userOptions.options.box );  //returns in 'spinePosition' the positions of all valid voxels
    
    Real dx = (mmfHeader.box[1]-mmfHeader.box[0]) / mmfHeader.gridSize[0];
    userOptions.options.distanceThreshold = dx/5;
    MMFContraction( spinePosition, direction, userOptions.options, userOptions.options.convergeFraction, userOptions.options.distanceThreshold, userOptions.options.eigenvalueThreshold, userOptions.options.maxLoop);  // contracts the filaments to a spine and the walls to a plane
    
    Array<Real,2> spineProperties( spinePosition.axisSize(0), 2 );
    environmentProperties( mass, spinePosition, userOptions.options, spineProperties, sizeData, massData, averageData );  // computes the properties
    
    
    
    
    //! output the environment properties into a binary MMF file
    mmfHeader.fileType = MMF_PROPERTIES;
    mmfHeader.gridSize[0] = spineProperties.axisSize(0);
    mmfHeader.gridSize[1] = spineProperties.axisSize(1);
    mmfHeader.gridSize[2] = 1; 
    mmfHeader.totalGrid = spineProperties.size();
    mmfHeader.updateObservations( ";   ");
    mmfHeader.updateObservations( argv, argc );
    string outFile = userOptions.outputFile + "_properties.MMF";
    writeMMF( spineProperties, mmfHeader, outFile, "filament/wall properties file");
    
    
    
    
    //! output in a text file the results
    string sizeFile = userOptions.outputFile + ".size";
    string massFile = userOptions.outputFile + ".mass";
    cout << "Writing the object size information to the '" << sizeFile << "' file and the object mass density information to the '" << massFile << "'\n.";
    
    fstream out;
    openOutputTextFile( out, sizeFile );
    out << "# This file contains average properties for Cosmic Web " << (userOptions.options.feature==3 ? "filaments" : "walls") << "\n.";
    if ( userOptions.options.feature==3 )
        out << "#\tAverage cross-section diameter: " << averageData[0] << "  (Mpc/h)\n"
            << "#\tAverage cross-section area    : " << averageData[0]*averageData[0]*3.14/4. << "  (Mpc/h)^2\n"
            << "#\tAverage length                : " << averageData[1] << "  Mpc/h per (Mpc/h)^3\n";
    else if ( userOptions.options.feature==2 )
        out << "#\tAverage wall height           : " << averageData[0] << "  (Mpc/h)\n"
            << "#\tAverage wall area             : " << averageData[1] << "  (Mpc/h)^2 per (Mpc/h)^3\n\n";
    out << "# The following columns contain the filament diameter / wall thickness information.\n"
        << "# The 5 columns give: 1st=the diameter/thickness value, 2nd=the length of filaments/ area of walls with the given value, 3rd=normalized 2nd column to the total filament length / wall area, 4th=cummulative mass as function of diameter/thickness and 5th=cummulative volume function.\n"
        << "#The results were obtained using the commands: " << userOptions.programOptions << "\n\n\n";
    
    for (size_t i=0; i<noBins; ++i)
    {
        out << (sizeData(i,0) + sizeData(i+1,0) ) / 2.;
        for (int j=1; j<sizeData.axisSize(1); ++j)
            out << "\t" << sizeData(i,j);
        out << "\n";
    }
    out.close();
    
    
    openOutputTextFile( out, massFile );
    out << "# The following columns contain the filament linear mass density / wall surface mass density information.\n"
        << "# The 3 columns give: 1st=the linear/surface mass density value, 2nd=the length of filaments/ area of walls with the given value and 3rd=normalized 2nd column to the total filament length / wall area.\n"
        << "#The results were obtained using the commands: " << userOptions.programOptions << "\n\n\n";
    
    for (size_t i=0; i<noBins; ++i)
    {
        out << sqrt(massData(i,0) * massData(i+1,0) );
        for (int j=1; j<massData.axisSize(1); ++j)
            out << "\t" << massData(i,j);
        out << "\n";
    }
    out.close();
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
            ("densityFile", po::value<string>(&(userOptions.densityFile)), "give the name of the file that stores the density field data")
            ("cleanFile", po::value<string>(&(userOptions.cleanFile)), "give the associated clean response file to the input spine direction file")
            ("nodeFile",po::value<string>(&(userOptions.nodeFile)), "node clean response file")
            ("filaFile",po::value<string>(&(userOptions.filaFile)), "filament clean response file")
            ("radius",po::value<Real>(&(userOptions.options.radius)), "radius in Mpc used for the smoothing and neighbor search")
            ("maxLoop",po::value<int>(&(userOptions.options.maxLoop)), "number of maximum loops for the contraction procedure")
            ("convergenceFraction",po::value<Real>(&(userOptions.options.convergeFraction)), "stop when the change in convergence fraction is less than this value")
            ("eigenvalueThreshold",po::value<Real>(&(userOptions.options.eigenvalueThreshold)), "thresholds for eigenvalues for when convergence is achieved")
            ("distanceThreshold",po::value<Real>(&(userOptions.options.distanceThreshold)), "thresholds for displacement value for when convergence is achieved")
            ("nonPeriodic","specify this option if the data is non-periodic")
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("inputFile", po::value<string>(&(userOptions.inputFile)), "name of the input spine direction file")
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
    cout << "Use this program to compute the properties of filaments and walls.\n";
    cout << "Usage:    " << argv[0] << "  input_direction_file  name_output_file 'options - see below' \n";
    cout << "On top of the above, the user can specify the following additional options:\n";
    cout << visibleOptions << "\n";
    exit(0);
}


/* Print to the user what options the program will use. */
void printOptions( User_options &userOptions )
{
    cout << "RUNNING:   " << userOptions.programOptions << "\n\n";
    cout << "The program will compute the MMF response for a set of scales using the following input parameters:\n"
        << "\t input direction file        : " << userOptions.inputFile << "\n"
        << "\t input clean response file   : " << userOptions.cleanFile << "\n"
        << "\t input density file          : " << userOptions.densityFile << "\n"
        << "\t output results text file    : " << userOptions.outputFile << "\n"
        << "\t clean node file (masking)   : " << (userOptions.nodeFileOn ? userOptions.nodeFile : "none given") << "\n"
        << "\t clean filament file(masking): " << (userOptions.filaFileOn ? userOptions.filaFile : "none given") << "\n";
    cout << "\t filter radius               : " << userOptions.options.radius << "  Mpc/h\n"
        << "\t distance threshold          : " << userOptions.options.distanceThreshold << "  Mpc/h\n"
        << "\t eigenvalue threshold        : " << userOptions.options.eigenvalueThreshold << "\n"
        << "\t boundary conditions         : " << (userOptions.options.periodic ? "periodic" : "free (non-periodic)") << "\n"
        << "\t maximum loops               : " << userOptions.options.maxLoop << "\n";
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
    
    if (not vm.count("densityFile") )
        throwError( "The program needs a density file. Specify the density field data via the '--densityFile' program option." );
    if (not vm.count("cleanFile") )
        throwError( "The program needs a clean response file. Specify the density field data via the '--cleanFile' program option." );
    
    if ( vm.count("nodeFile") ) userOptions->nodeFileOn = true;
    if ( vm.count("filaFile") ) userOptions->filaFileOn = true;
    if ( vm.count("nonPeriodic") ) userOptions->options.periodic = false;
    
    
    // set values to userOptions->programOptions
    for (int i=0; i<argc; ++i)
        userOptions->programOptions += string( argv[i] ) + " ";
    
    
    printOptions( *userOptions );
}




