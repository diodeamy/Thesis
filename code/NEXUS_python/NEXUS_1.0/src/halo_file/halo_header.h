#ifndef HALO_HEADER
#define HALO_HEADER

#include <string>
#include <defines.h>

static int const fillSize3 = 1024 - 4*8 - 10*8 - 4*8 - 2*8 ;
static int const columnNameLength = 16;
static int const maximumColumns = 100;

// data structure that stores information in the header of a halo binary file output
struct Halo_header
{
    // Information about the halo file
    size_t  noHalos;            // the number of halos in the file
    size_t  noColumnsIntegers;  // the number of halos in the file
    size_t  noColumnsFloats;    // the number of halos in the file
    size_t  noColumns;          // the number of halos in the file
    
    double  mpcUnit;            // the value of 1Mpc in halo position units units
    double  box[6];             // keep track of the box coordinates (xMin, xMax, yMin, yMax, zMin, zMax)
    size_t  positionColumns[3]; // gives the columns in the floating point data set which stores the positions x,y,z of the halos
    
    double  massUnit;           // the mass unit in Msolar/h
    double  massRange[2];       // min and max mass of the halos in the file
    size_t  massColumn;         // gives the column in the floating point data set which gives the halo mass
    
    size_t  noFiles;            // keeps track of the number of files
    char    fill[fillSize3];    // fill to 1024 bytes - left 760 - used to keep track of information on how the file was obtained
    size_t  FILE_ID;            // unique ID for this file type
    
    
    Halo_header();              // constructor - initializes to 0 or to non-assigned value (=-1) 
    void print();               // print the content of the halo header to the standard output
    double boxVolume();         // returns the volume of the box for the given data
    
    void updateObservations(char **obs, int const size);       // update the observation field
    void updateObservations(std::string newObservations);      // update the observation field
    void overwriteObservations(std::string newObservations);   // overwrite in the observation field
    void updateColumnNames(char *columNames, std::string newName);   // update a new column name at the position given by the char pointer
    
};


#endif
