#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "output.h"
#include <miscellaneous.h>
using namespace std;


/* This function outputs the 3D density map to a binary file plus a small ASCII file which describes the contents of the binary file. */
void outputVisivo(float *delta,
                  Int const gridSize[3],
                  string rootName)
{
    //! write the description file
    string outputDataFile = rootName + ".raw";      //name of output binary file
    string descriptionFile = rootName + ".desc";    //name of description file
    
    cout << "Writing the result in a VisiVo accessible format. Writing the description ASCII file '" << descriptionFile << "' and binary raw file '" << outputDataFile << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, descriptionFile );   //open the description file
    
    outputFile << "rawGridsDesc\n"
            "density\n"
            "3\n"
            "Float\n" <<
            gridSize[0] << "\n" <<
            gridSize[1] << "\n" <<
            gridSize[2] << "\n"
            "time\n"
            "l\n"
            "1.0\t" << outputDataFile;
    outputFile.close();
    
    
    //! write the binary data file
    Int const totalSize = gridSize[0]*gridSize[1]*gridSize[2];
    outputRaw( delta, totalSize, outputDataFile );
    cout << "Done.\n";
}


/* This function takes the logarithm of a field: log(delta/delta0). */
void logarithm(Real *delta,
               Int const totalSize,
               Real const delta0,
               Real const minimumDelta)
{
    Real const temp = log10( minimumDelta/delta0 );
    for (Int i=0; i<totalSize; ++i)
    {
        if ( delta[i]>=minimumDelta ) 
            delta[i] = log10( delta[i]/delta0 );
        else 
            delta[i] = temp;
    }
}

