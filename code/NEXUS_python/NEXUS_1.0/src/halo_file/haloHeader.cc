#include <iostream>
#include <string>
#include <sstream>
#include <cstring>

#include <miscellaneous.h>
#include "halo_header.h"
using namespace std;


/* Constructor for halo header. */
Halo_header::Halo_header()
{
    noHalos = 0;
    noColumnsIntegers = 0;
    noColumnsFloats = 0;
    noColumns = 0;
    
    mpcUnit = 1;
    for (int i=0; i<6; ++i) box[i] = 0.;
    for (int i=0; i<3; ++i) positionColumns[i] = -1;
    
    massUnit = 1;
    for (int i=0; i<2; ++i) massRange[i] = 0.;
    massColumn = -1;
    
    noFiles = 1;
    initializeCharArray( fill, sizeof(fill) );
    FILE_ID = 100;
}

/* This function prints the contents of a halo header object. */
void Halo_header::print()
{
    cout << "\nThe header of the halo file contains the following info:\n" 
        << "noHalos       = " << noHalos << "\n"
        << "noColumnsIntegers = " << noColumnsIntegers << "\n"
        << "noColumnsFloats   = " << noColumnsFloats << "\n"
        << "noColumns     = " << noColumns << "\n"
        << "mpcUnit       = " << mpcUnit << "\n"
        << "box coords    = " << box[0] << "  " << box[1] << "  " << box[2] << "  " << box[3] << "  " << box[4] << "  " << box[5] << "\n"
        << "positionColumns   = " << positionColumns[0] << "  " << positionColumns[1] << "  " << positionColumns[2] << "\n"
        << "massUnit      = " << massUnit << "\n"
        << "massRange     = " << massRange[0] << "  " << massRange[1] << "\n"
        << "massColumn    = " << massColumn << "\n"
        << "noFiles       = " << noFiles << "\n"
        << "fill          = " << fill << "\n\n";
}

/* Returns the volume of the box. */
double Halo_header::boxVolume()
{
    double result = 1.;
    for (int i=0; i<3; ++i)
        result *= (box[2*i+1] - box[2*i]);
    if (result<0.)
        throwError( "Negative volume in function 'halo_header::boxVolume()'. Check the box coordinates for inconsistencies." );
    return result;
}

/* Updates the 'observations'. */
void Halo_header::updateObservations(string newObservations)
{
    string oldObs( fill );
    oldObs += newObservations + " ;  ";
    copyStringEnd( oldObs.c_str(), oldObs.length(), fill, sizeof(fill), "The updated observations in the halo header had a too large size so only the ending part was kept. Some observations are lost." );
}

/* Updates the 'observations'. */
void Halo_header::updateObservations(char **obs, int const size)
{
    ostringstream buffer;
    buffer << getFilename( string(obs[0]) );
    for (int i=1; i<size; ++i)
        buffer << " " << obs[i];
    this->updateObservations( buffer.str() );
}

/* Replaces the 'observations'. */
void Halo_header::overwriteObservations(string newObservations)
{
    newObservations += " ;  ";
    copyStringFront( newObservations.c_str(), newObservations.length(), fill, sizeof(fill), "The new observations in the halo header had a too large size so only the beginning part was kept. Some observations are lost." );
}


// update a new column name at the position given by the char pointer
void Halo_header::updateColumnNames(char *columNames,
                                    std::string newName)
{
    initializeCharArray( columNames, columnNameLength );
    copyStringFront( newName.c_str(), newName.length(), columNames, columnNameLength, "The new column name of halo properties had a too large size so only the beginning part was kept." );
}

