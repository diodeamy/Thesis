#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>

#include <density_header.h>
#include <miscellaneous.h>
#include "MMF_header.h"
using namespace std;


MMF_header::MMF_header()
{
    for (int i=0; i<3; ++i)
    {
        gridSize[i] = size_t(0);
        MMF_fileGrid[i] = 1;
    }
    totalGrid = size_t(0);
    radius = -10.;
    scale = -10;
    feature = -10;
    fileType = -10;
    filter = -10;
    bias = 1.;
    noMMF_files = 1;
    indexMMF_file = -1;
    method = -1;
    for (int i=0; i<6; ++i)
        box[i] = 0.;
    
    
    for (int i=0; i<6; ++i)
    {
        npartTotal[i] = size_t(0);
        mass[i] = 0.;
    }
    time = 0.;
    redshift = 0.;
    BoxSize = 0.;
    Omega0 = 0.;
    OmegaLambda = 0.;
    HubbleParam = 0.;
    
    FILE_ID = 10;
    initializeCharArray( fill, sizeof(fill) );
}




/* This function prints the contents of a density header object. */
void MMF_header::print()
{
    string featureName, fileTypeName, biasValue;
    switch(feature)
    {
        case 4: featureName = "4 - blob"; break;
        case 3: featureName = "3 - filament"; break;
        case 2: featureName = "2 - wall"; break;
        default: featureName = "unknown";
    }
    switch(fileType)
    {
        case 1: fileTypeName = "1 - MMF response for a single scale"; break;
        case 2: fileTypeName = "2 - Hessian eigenvalues and eigenvectors for a given scale"; break;
        case 3: fileTypeName = "3 - maximum MMF response over a given range of scales"; break;
        case 4: fileTypeName = "4 - Hessian eigenvalues and eigenvectors corresponding to the maximum response"; break;
        default: fileTypeName = "unknown";
    }
    if ( scale!= -1 )
        biasValue = "not used";
    else
    {
        ostringstream buffer;
        buffer << bias ;
        biasValue = buffer.str();
    }
    
    cout << "\nThe header of the MMF output file contains the following info:\n"
            "1) Information about the actual MMF file:\n"
            << "gridSize      = " << gridSize[0] << "  " << gridSize[1] << "  " << gridSize[2] << "\n"
            << "totalGrid     = " << totalGrid << "\n"
            << "feature       = " << featureName << "\n"
            << "scale         = " << scale << "\n"
            << "radius        = " << radius << "\n"
            << "bias          = " << biasValue << "\n"
            << "filter        = " << filter << "\n"
            << "file type     = " << fileTypeName << "\n"
            << "no MMF files  = " << noMMF_files << "\n";
    if ( noMMF_files!=1 )
        cout << "MMF file grid = " << MMF_fileGrid[0] << "  " << MMF_fileGrid[1] << "  " << MMF_fileGrid[2] << "\n"
                << "index MMF file= " << indexMMF_file << "\n";
    cout << "box coords    = " << box[0] << "  " << box[1] << "  " << box[2] << "  " << box[3] << "  " << box[4] << "  " << box[5] << "\n";
    
    
    cout << "\n2) Information about the snapshot file used to compute this data:\n"
            << "npartTotal[6] =  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "\n"
            << "mass[6]       =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time          =  " << time << "\n"
            << "redshift      =  " << redshift << "\n"
            << "BoxSize       =  " << BoxSize << "\n"
            << "Omega0        =  " << Omega0 << "\n"
            << "OmegaLambda   =  " << OmegaLambda << "\n"
            << "HubbleParam   =  " << HubbleParam << "\n";
    
    
    cout << "\n3) Additional information contained in the MMF header:\n"
            << "fill          = " << fill << "\n";
}


/* Checks if two MMF headers are compatible: the two represent data with the same grid size. */
void MMF_header::compatible(MMF_header &header)
{
    if ( gridSize[0]!=header.gridSize[0] or gridSize[1]!=header.gridSize[1] or gridSize[2]!=header.gridSize[2] or totalGrid!=header.totalGrid )
    {
        cout << "~~~ ERROR ~~~ The two MMF files do not have grids of the same dimension. Because of that they cannot both be used to compute the maximum MMF response. The program ended unsuccesfully!\n";
        exit( EXIT_FAILURE );
    }
}
void MMF_header::compatible(size_t *grid)
{
    if ( gridSize[0]!=grid[0] or gridSize[1]!=grid[1] or gridSize[2]!=grid[2] )
    {
        cout << "~~~ ERROR ~~~ In function 'MMF_header::compatible'. The two grids do not have the same dimension. The program ended unsuccesfully!\n";
        exit( EXIT_FAILURE );
    }
}
/* Checks if the MMF header scale is the same as the argument. */
void MMF_header::checkScale(int const scaleNo)
{
    if ( this->scale!=scaleNo )
    {
        cout << "~~~ ERROR ~~~ Error found when checking the scale of an MMF header. The expected scale of " <<scaleNo << " is different from the scale found in the header " << this->scale << ". The program ended unsuccesfully!\n";
        exit( EXIT_FAILURE );
    }
}

/* Get the total number of grid cells of the response. */
Int MMF_header::totalSize()
{
    return gridSize[0] * gridSize[1] * gridSize[2];    //compute it like this to work for very large grids
}
Int MMF_header::dataSize()
{
    return Int(gridSize[0] * gridSize[1] * gridSize[2]);    //compute it like this to work for very large grids
}




//! Functions that write values into the MMF header
/* Copy the shared information from the density header into the MMF one. */
void MMF_header::copyDensityHeader( Density_header const &densityHeader )
{
    for (int i=0; i<6; ++i)
    {
        mass[i] = densityHeader.mass[i];
        npartTotal[i] = densityHeader.npartTotal[i];
        box[i] = densityHeader.box[i];
    }
    time = densityHeader.time;
    redshift = densityHeader.redshift;
    BoxSize = densityHeader.BoxSize;
    Omega0 = densityHeader.Omega0;
    OmegaLambda = densityHeader.OmegaLambda;
    HubbleParam = densityHeader.HubbleParam;
    
    for (int i=0; i<3; ++i)
        gridSize[i] = densityHeader.gridSize[i];
    totalGrid = densityHeader.totalGrid;
}
/* Return the value of the filter radius for the given scale. */
float MMF_header::filterRadius(float const radius0,
                               float const base,
                               int const scaleNo)
{
    return radius0 * pow( base, scaleNo );
}
/* Update 'scale' and 'radius'. */
void MMF_header::updateScale(int const scaleNo,
                             float const radius0,
                             float const base)
{
    scale = scaleNo;
    radius = this->filterRadius( radius0, base, scale );
}
/* Update 'scale' and 'radius'. */
void MMF_header::updateFeature(int const filterType)
{
    feature = int(filterType/10);
    filter = filterType;
}
/* Update 'fileType'. */
void MMF_header::updateFileType(int const file_type)
{
    fileType = file_type;
}
/* Update 'bias'. */
void MMF_header::updateBias(float const biasValue)
{
    bias = biasValue;
}
/* Update 'MMF_fileGrid' and 'noMMF_files' variables. */
void MMF_header::updateMMF_fileGrid(int const grid[],
                                    int const size)
{
    if ( size!=3 ) throwError( "The size of the input array for function 'MMF_header::updateMMF_fileGrid' must be 3." );
    noMMF_files = 1;
    for (int i=0; i<3; ++i)
    {
        MMF_fileGrid[i] = grid[i];
        noMMF_files *= grid[i];
    }
}
void MMF_header::updateMMF_fileGrid(int const nX,
                                    int const nY,
                                    int const nZ)
{
    noMMF_files = nX * nY * nZ;
    MMF_fileGrid[0] = nX;
    MMF_fileGrid[1] = nY;
    MMF_fileGrid[2] = nZ;
}
void MMF_header::updateMMF_fileIndex(int const fileIndex)
{
    lowerBoundCheck( noMMF_files, 1, "'noMMF_files' variable of class 'MMF_header' in function 'MMF_header::updateMMF_fileIndex'" );
    intervalCheck( fileIndex, 0, noMMF_files-1, "'fileIndex' argument of function 'MMF_header::updateMMF_fileIndex'" );
    indexMMF_file = fileIndex;
}
void MMF_header::updateMaximumResponse()
{
    scale = -1;
    radius = -1.;
    fileType = MMF_MAX_RESPONSE;
}
/* Returns the volume of the box. */
double MMF_header::boxVolume()
{
    double result = 1.;
    for (int i=0; i<3; ++i)
        result *= (box[2*i+1] - box[2*i]);
    if (result<0.)
        throwError( "Negative volume in function 'MMF_header::boxVolume()'. Check the box coordinates for inconsistencies." );
    return result;
}
/* Returns the volume of a grid cell. */
double MMF_header::gridCellVolume()
{
    return this->boxVolume() / double(this->dataSize());
}
/* Returns true if there are no values set for the box coordinates for the given density data. */
bool MMF_header::emptyBox()
{
    if ( this->boxVolume()==double(0.) )
        return true;
    return false;
}




//! Functions returning file names
/* This gives the name of the MMF response file for a given scale. */
string MMF_header::responseFilename(string const &rootName,
                                    int const scale,
                                    int const fileIndex)
{
    if ( noMMF_files==1 )
    {
        ostringstream buffer;
        buffer << rootName << "_scale=" << scale << ".MMF";
        return buffer.str();
    }
    if ( noMMF_files<1 ) throwError( "The 'noMMF_files' variable in the MMF header wasn't initialized to the correct value (>=1). Error in function 'MMF_header::responseFilename'." );
    ostringstream buffer;
    buffer << rootName << "_scale=" << scale << ".MMF." << fileIndex;
    return buffer.str();
}
/* This gives the name of the MMF maximum response file. */
string MMF_header::maxResponseFilename(string const &rootName,
                                       int const fileIndex)
{
    if ( noMMF_files==1 )
    {
        ostringstream buffer;
        buffer << rootName << "_maxResponse.MMF";
        return buffer.str();
    }
    if ( noMMF_files<1 ) throwError( "The 'noMMF_files' variable in the MMF header wasn't initialized to the correct value (>=1). Error in function 'MMF_header::maxResponseFilename'." );
    ostringstream buffer;
    buffer << rootName << "_maxResponse.MMF." << fileIndex;
    return buffer.str();
}
/* This gives the name of the Hessian eigenvalue and eigenvectors file for a given scale. */
string MMF_header::eigenFilename(string const &rootName,
                                 int const scale,
                                 int const fileIndex)
{
    if ( noMMF_files==1 )
    {
        ostringstream buffer;
        buffer << rootName << "_scale=" << scale << "_eigen.MMF";
        return buffer.str();
    }
    if ( noMMF_files<1 ) throwError( "The 'noMMF_files' variable in the MMF header wasn't initialized to the correct value (>=1). Error in function 'MMF_header::eigenFilename'." );
    ostringstream buffer;
    buffer << rootName << "_scale=" << scale << "_eigen.MMF." << fileIndex;
    return buffer.str();
}
/* This gives the name of the MMF maximum response file. */
string MMF_header::maxEigenFilename(string const &rootName,
                                    int const fileIndex)
{
    if ( noMMF_files==1 )
    {
        ostringstream buffer;
        buffer << rootName << "_maxEigen.MMF";
        return buffer.str();
    }
    if ( noMMF_files<1 ) throwError( "The 'noMMF_files' variable in the MMF header wasn't initialized to the correct value (>=1). Error in function 'MMF_header::maxEigenFilename'." );
    ostringstream buffer;
    buffer << rootName << "_maxEigen.MMF." << fileIndex;
    return buffer.str();
}

/* This function checks that all the MMF response files corresponding to the scales in the scales array exist. */
void MMF_header::checkResponseFiles(string const &mmfRootName,
                                    int scale[],
                                    int const noScales)
{
    for (int i=0; i<noScales; ++i)
    {
        string filename = this->responseFilename( mmfRootName, scale[i], 0 );
        existentFile( filename );
    }
}
/* This function checks that all the Hessian eigenvalues and eigenvectors files corresponding to the scales in the scales array exist. */
void MMF_header::checkEigenFiles(string const &mmfRootName,
                                 int scale[],
                                 int const noScales)
{
    for (int i=0; i<noScales; ++i)
    {
        string filename = this->eigenFilename( mmfRootName, scale[i], 0 );
        existentFile( filename );
    }
}
/* This function checks that all the files of an output saved in multiple files exist. */
void MMF_header::checkFiles(string const &mmfRootName)
{
    if ( noMMF_files<1 ) throwError( "The 'noMMF_files' variable in the MMF header wasn't initialized to the correct value (>=1). Error in function 'MMF_header::checkFiles'." );
    for (int i=0; i<noMMF_files; ++i)
    {
        string filename;
        switch(fileType)
        {
            case MMF_RESPONSE:
                filename = this->responseFilename( mmfRootName, scale, i );
                break;
            case MMF_EIGEN:
                filename = this->eigenFilename( mmfRootName, scale, i );
                break;
            case MMF_MAX_RESPONSE:
                filename = this->maxResponseFilename( mmfRootName, i );
                break;
            case MMF_MAX_EIGEN:
                filename = this->maxEigenFilename( mmfRootName, i );
                break;
            default:
                throwError( "Unrecognized value for the MMF header variable 'fileType' in function 'MMF_header::checkFiles'." );
        }
        existentFile( filename );
    }
}



//! Observations
/* Updates the 'observations'. */
void MMF_header::updateObservations(string newObservations)
{
    string oldObs( fill );
    oldObs += newObservations;
    copyStringEnd( oldObs, fill, sizeof(fill), "The updated observations in the MMF header had a too large size so only the end part was kept. Some observations are lost." );
}
void MMF_header::updateObservations(char **obs, int const size)
{
    string temp = "";
    for (int i=0; i<size; ++i)
        temp += string(obs[i]) + " ";
    this->updateObservations( temp );
}

/* Replaces the 'observations'. */
void MMF_header::overwriteObservations(string newObservations)
{
    copyStringFront( newObservations, fill, sizeof(fill), "The new observations in the MMF header had a too large size so only the beginning part was kept. Some observations are lost." );
}
void MMF_header::overwriteObservations(char **obs, int const size)
{
    string temp = "";
    for (int i=0; i<size; ++i)
        temp += string(obs[i]) + " ";
    this->overwriteObservations( temp );
}


