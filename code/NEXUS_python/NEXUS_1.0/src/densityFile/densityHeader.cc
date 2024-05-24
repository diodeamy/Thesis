#include <iostream>
#include <string>
#include <sstream>
#include <cstring>

#include <gadget_header.h>
#include <miscellaneous.h>
#include "density_header.h"
using namespace std;


/* Constructor for density header. */
Density_header::Density_header()
{
    for (int i=0; i<3; ++i)
    {
        gridSize[i] = size_t(0);
        densityFileGrid[i] = 1;
    }
    totalGrid = size_t(0);
    fileType = UNKNOW_FILE;
    noDensityFiles = 1;
    indexDensityFile = -1;
    for (int i=0; i<6; ++i)
        box[i] = 0.;
    
    time = -1.;
    redshift = -1.;
    
    method = UNKNOW_METHOD;
    FILE_ID = 1;
    initializeCharArray( fill, sizeof(fill) );
}


/* This function prints the contents of a density header object.
*/
void Density_header::print()
{
    string densityMethod;
    if ( method==DTFE_METHOD ) densityMethod = "DTFE";
    else if ( method==TSC_METHOD ) densityMethod = "TSC";
    else if ( method==SPH_METHOD ) densityMethod = "SPH";
    else densityMethod = "unknown";
    
    string fileTypeName = "unknown file type";
    if ( fileType==DENSITY_FILE ) fileTypeName = "the file stores a density field";
    else if ( fileType==VELOCITY_FILE ) fileTypeName = "the file stores a velocity field";
    else if ( fileType==VELOCITY_GRADIENT_FILE ) fileTypeName = "the file stores the gradient of a velocity field";
    else if ( fileType==VELOCITY_DIVERGENCE_FILE ) fileTypeName = "the file stores a velocity divergence";
    else if ( fileType==VELOCITY_SHEAR_FILE ) fileTypeName = "the file stores a velocity shear";
    else if ( fileType==VELOCITY_VORTICITY_FILE ) fileTypeName = "the file stores a velocity vorticity";
    else if ( fileType==SCALAR_FIELD_FILE ) fileTypeName = "the file stores a scalar field";
    else if ( fileType==SCALAR_FIELD_GRADIENT_FILE ) fileTypeName = "the file stores the gradient of a scalar field";
    
    
    cout << "\nThe header of the density file contains the following info:\n" <<
            "1) Information about the actual density computations:\n"
            << "gridSize      = " << gridSize[0] << "  " << gridSize[1] << "  " << gridSize[2] << "\n"
            << "totalGrid     = " << totalGrid << "\n"
            << "file type     = " << fileTypeName << "\n"
            << "# density file= " << noDensityFiles << "\n";
    if ( noDensityFiles>1 )
        cout << "file grid size= " << densityFileGrid[0] << "  " << densityFileGrid[1] << "  " << densityFileGrid[3] << "\n"
            << "file index    = " << indexDensityFile << "\n";
    cout << "box coords    = " << box[0] << "  " << box[1] << "  " << box[2] << "  " << box[3] << "  " << box[4] << "  " << box[5] << "\n";
            
    cout << "\n2) Information about the snapshot file used to compute the density:\n"
            << "npartTotal[6] =  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "\n"
            << "mass[6]       =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time          =  " << time << "\n"
            << "redshift      =  " << redshift << "\n"
            << "BoxSize       =  " << BoxSize << "\n"
            << "Omega0        =  " << Omega0 << "\n"
            << "OmegaLambda   =  " << OmegaLambda << "\n"
            << "HubbleParam   =  " << HubbleParam << "\n";
    
    
    cout << "\n3) Information about files and additional remarks:\n"
            << "method          = " << densityMethod << "\n"
            << "fill            = " << fill << "\n\n";
}


/* This function copies some of the contents of the gadget header into the density header.
*/
void Density_header::copyGadgetHeader(Gadget_header const & gadgetHeader)
{
    for (int i=0; i<6; ++i)
    {
        npartTotal[i] = gadgetHeader.npartTotal[i];
        mass[i] = gadgetHeader.mass[i];
    }
    time = gadgetHeader.time;
    redshift = gadgetHeader.redshift;
    BoxSize = gadgetHeader.BoxSize;
    Omega0 = gadgetHeader.Omega0;
    OmegaLambda = gadgetHeader.OmegaLambda;
    HubbleParam = gadgetHeader.HubbleParam;
}

/* Check if two density headers have the same grid size. */
void Density_header::checkSimilar(Density_header &other)
{
    for (int i=0; i<4; ++i)
        if ( gridSize[i]!=other.gridSize[i] )
            throwError( "The two density headers are not similar since they do not have the same grid dimensions." );
}
bool Density_header::isSimilar(Density_header &other)
{
    for (int i=0; i<4; ++i)
        if ( gridSize[i]!=other.gridSize[i] )
            return false;
    return true;
}

/* Update density grid size information. */
void Density_header::updateGridSize(size_t const sizeX,
                                    size_t const sizeY,
                                    size_t const sizeZ)
{
    totalGrid = sizeX * sizeY * sizeZ;
    gridSize[0] = sizeX;
    gridSize[1] = sizeY;
    gridSize[2] = sizeZ;
}

/* Returns the total number of data pointed by the density file. */
size_t Density_header::dataSize()
{
    return size_t(gridSize[0]) * size_t( gridSize[1]*gridSize[2] );
}
/* Returns the volume of the box. */
double Density_header::boxVolume()
{
    double result = 1.;
    for (int i=0; i<3; ++i)
        result *= (box[2*i+1] - box[2*i]);
    if (result<0.)
        throwError( "Negative volume in function 'Density_header::boxVolume()'. Check the box coordinates for inconsistencies." );
    return result;
}
/* Returns true if there are no values set for the box coordinates for the given density data. */
bool Density_header::emptyBox()
{
    if ( this->boxVolume()==double(0.) )
        return true;
    return false;
}
/* Returns the mass corresponding to a density=1 cell. */
double Density_header::massInCell()
{
    double result = 1.;
#ifndef EXTERNAL_COMPILATION
    result = RHO_CRITICAL * Omega0 * this->boxVolume() / totalGrid;
#endif
    return result;
}


void Density_header::noneFileGrid()
{
    for (int i=0; i<3; ++i)
        densityFileGrid[i] = 1;
    noDensityFiles = 1;
}



/* Updates the 'observations'. */
void Density_header::updateObservations(string newObservations)
{
    string oldObs( fill );
    oldObs += newObservations;
    copyStringFront( oldObs.c_str(), oldObs.length(), fill, sizeof(fill), "The updated observations in the density header had a too large size so only the beginning part was kept. Some observations are lost." );
}


/* Updates the 'observations'. */
void Density_header::updateObservations(char **obs, int const size)
{
    for (int i=0; i<size; ++i)
        this->updateObservations( string(obs[i])+" ;  " );
}


/* Replaces the 'observations'. */
void Density_header::overwriteObservations(string newObservations)
{
    copyStringFront( newObservations.c_str(), newObservations.length(), fill, sizeof(fill), "The new observations in the density header had a too large size so only the beginning part was kept. Some observations are lost." );
}



//! Functions used for writting/reading the density into multiple files
/* Outputs the name of the 'i'-th file for density saved in multiple files. */
string Density_header::filename(string &rootName, int const i)
{
    ostringstream buffer;
    buffer << rootName << i;
    return buffer.str();
}

/* Returns the secondary density array size for file 'fileNo'. */
size_t Density_header::secondaryTotalSize(int const fileNo)
{
    size_t indices[3], temp = 1;
    this->secondarySize( fileNo, indices );
    for (int i=0; i<3; ++i)
        temp *= indices[i];
    return temp;
}

/* Get boundary values in spatial units for multiple files. */
void Density_header::boundariesUnits(int const fileNo, Real *ptr)
{
    if ( this->emptyBox() ) throwError( "Cannot return the boundaries of the secondary box in function 'Density_header::boundariesUnits' since the 'Density_header::box' member wasn't initialized to the coordinates of the box whose density data is represented here." );
    Real temp[6];
    this->boundariesBoxLength( fileNo, temp );
    for (int i=0; i<3; ++i)
    {
        double length = box[2*i+1] - box[2*i];
        ptr[2*i] = box[2*i] + length*temp[2*i];
        ptr[2*i+1] = box[2*i] + length*temp[2*i+1];
    }
}

/* Get boundary values relative to boxLength for multiple files. */
void Density_header::boundariesBoxLength(int const fileNo, Real *ptr)
{
    size_t indices[6];
    this->boundariesIndices( fileNo, indices );
    for (int i=0; i<3; ++i)
    {
        ptr[2*i] = indices[2*i] / Real(gridSize[i]);
        ptr[2*i+1] = indices[2*i+1] / Real(gridSize[i]);
    }
}

/* Copy the density stored in file 'fileNo' to the main density array. */
void Density_header::copyDensityToMain(Real *mainArray,
                                       Real const *secondaryArray,
                                       int const fileNo)
{
    size_t indices[6], secIndex[3];
    this->boundariesIndices( fileNo, indices );
    this->secondarySize( fileNo, secIndex );
    
    for (size_t i1=indices[0]; i1<indices[1]; ++i1)
        for (size_t i2=indices[2]; i2<indices[3]; ++i2)
            for (size_t i3=indices[4]; i3<indices[5]; ++i3)
            {
                size_t temp1 = i1*gridSize[1]*gridSize[2] + i2*gridSize[2] + i3;
                size_t temp2 = (i1-indices[0])*secIndex[1]*secIndex[2] + (i2-indices[2])*secIndex[2] + i3-indices[4];
                mainArray[temp1] = secondaryArray[temp2];
            }
}

/* Copy from the main density array to the array that will be stored in file 'fileNo'. */
void Density_header::copyDensityToSecondary(Real *secondaryArray,
                                            int const fileNo,
                                            Real const *mainArray)
{
    size_t indices[6], secIndex[3];
    this->boundariesIndices( fileNo, indices );
    this->secondarySize( fileNo, secIndex );
    
    for (size_t i1=0; i1<secIndex[0]; ++i1)
        for (size_t i2=0; i2<secIndex[1]; ++i2)
            for (size_t i3=0; i3<secIndex[2]; ++i3)
            {
                size_t temp1 = i1*secIndex[1]*secIndex[2] + i2*secIndex[2] + i3;
                size_t temp2 = (i1+indices[0])*gridSize[1]*gridSize[2] + (i2+indices[2])*gridSize[2] + i3+indices[4];
                secondaryArray[temp1] = mainArray[temp2];
            }
}





/* Update density grid size information. */
void Density_header::updateGridSize(size_t size[3])
{
    totalGrid = 1;
    for (int i=0; i<3; ++i)
    {
        gridSize[i] = int(size[i]);
        totalGrid *= gridSize[i];
    }
}

/* Update the variable 'densityFilename'. */
void Density_header::updateDensityFilename(std::string densityFileName, size_t const fileGrid[])
{
    for (int i=0; i<3; ++i)
        densityFileGrid[i] = int(fileGrid[i]);
    noDensityFiles = densityFileGrid[0] * densityFileGrid[1] * densityFileGrid[2];
}

/* Returns the secondary density array size for file 'fileNo'. */
void Density_header::secondarySize(int const fileNo, size_t *ptr)
{
    size_t indices[6];
    this->boundariesIndices( fileNo, indices );
    for (int i=0; i<3; ++i)
        ptr[i] = indices[2*i+1] - indices[2*i];
}

/* Get boundary indices for multiple files. */
void Density_header::boundariesIndices(int const fileNo, size_t *ptr)
{
    if ( fileNo>=noDensityFiles ) throwError( "You asked for a density file number that is larger than the total number of multiple density files. Error found in function 'Density_header::boundariesIndices'." );
    int temp[3];
    temp[0] = fileNo / (densityFileGrid[1]*densityFileGrid[2]);
    temp[1] = (fileNo - temp[0]*densityFileGrid[1]*densityFileGrid[2]) / densityFileGrid[2];
    temp[2] = fileNo - temp[0]*densityFileGrid[1]*densityFileGrid[2] - temp[1]*densityFileGrid[2];
    
    for (int i=0; i<3; ++i)
    {
        ptr[2*i] = size_t( gridSize[i]/densityFileGrid[i] * temp[i] );
        ptr[2*i+1] = size_t( gridSize[i]/densityFileGrid[i] * (temp[i]+1) );
        if ( temp[i]==densityFileGrid[i]-1 )
            ptr[2*i+1] = size_t( gridSize[i] );
    }
}




