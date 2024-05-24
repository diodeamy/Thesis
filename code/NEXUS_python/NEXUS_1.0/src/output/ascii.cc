#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include <defines.h>
#include "output.h"
#include <miscellaneous.h>
#include <hopGroup.h>
using namespace std;



/* Writes a 2D density map to an ASCI file. */
void outputASCII_matrix(Real *densityMatrix,
                        Int const nGrid[],
                        string outputFileName,
                        double const limitDown,
                        double const limitUp,
                        string const inputFileName)
{
    cout << "Writing a 2D matrix into the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    int direction, nx, ny;
    if ( nGrid[0]==1 )
    {
        direction = 0;
        nx = 1;
        ny = 2;
    }
    else if ( nGrid[1]==1 )
    {
        direction = 1;
        nx = 0;
        ny = 2;
    }
    else if ( nGrid[2]==1 )
    {
        direction = 2;
        nx = 0;
        ny = 1;
    }
    else throwError( "When writting the 2D slice in function 'outputASCII_matrix'. At least one of the grid sizes supplied to the function via the argument 'nGrid' must be one to give the direction of the slice." );
    
    // print helpfull information for the user
    outputFile << "# This file gives the matter density distribution along a 2D slice. The 2D slice is given by a " << nGrid[0] << " * " << nGrid[1] << " * " << nGrid[2] << " grid.\n"
            "# The data was loaded from the file '" << inputFileName << "'.\n"
            "# This slice was taken along the " << direction << " direction (x=0, y=1 and z=2), has a thickness of " << (limitUp-limitDown) << " Mpc and has the center at " << (limitDown+limitUp)/2. << " Mpc.\n\n";
    
    // print the output
    Int index;
    for (Int i1=0; i1<nGrid[nx]; ++i1)
    {
        for (Int i2=0; i2<nGrid[ny]; ++i2)
        {
            index = i1 * nGrid[ny] + i2;
            outputFile << densityMatrix[index] << "\t";
        }
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}



/* Outputs the haloes position, velocities and mass to an ASCII file. */
void outputASCII_Haloes(list<HopGroup> &blobHaloes,
                        string outputFileName,
                        string descriptionOfData,
                        string programOptions)
{
    cout << "Writing halo information to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    outputFile << "# This file contains " << descriptionOfData << " data. The first columns gives the halo tag, the second gives the number of particles while the 3rd gives the halo mass (in 10^10 solar masses), the 4th to 6th columns give the halo position (in Mpc) along the x, y and z axis while the 7th to 9th column give the halo CM velocities (in km/s).\n"
            << "# The data was obtain via the command: " << programOptions << "\n\n"
            << "# The next line gives the number of haloes in the file: \n"
            << blobHaloes.size() << "\n\n";
    
    for (list<HopGroup>::iterator it=blobHaloes.begin(); it!=blobHaloes.end(); ++it)
    {
        outputFile << it->id << "\t" << it->noParticles << "\t" << it->mass;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << "\t" << it->pos[i]/MPC_UNIT;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << "\t" << it->vel[i];
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}
void outputASCII_Haloes(HopGroup *blobHaloes,
                        Int const noGroups,
                        string outputFileName,
                        string descriptionOfData,
                        string programOptions)
{
    cout << "Writing halo information to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    outputFile << "# This file contains " << descriptionOfData << " data. The first columns gives the halo tag, the second gives the number of particles while the 3rd gives the halo mass (in 10^10 solar masses), the 4th to 6th columns give the halo position (in Mpc) along the x, y and z axis while the 7th to 9th column give the halo CM velocities (in km/s).\n"
            << "# The data was obtain via the command: " << programOptions << "\n\n"
            << "# The next line gives the number of haloes in the file: \n"
            << noGroups << "\n\n";
    
    for (Int j=0; j<noGroups; ++j)
    {
        outputFile << blobHaloes[j].id << "\t" << blobHaloes[j].noParticles << "\t" << blobHaloes[j].mass;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << "\t" << blobHaloes[j].pos[i]/MPC_UNIT;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << "\t" << blobHaloes[j].vel[i];
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}
void outputASCII_Haloes(list<HopGroup> &blobHaloes,
                        string outputFileName)
{
    cout << "Writing halo information to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    outputFile << "id, noParticles, mass, posX, posY, posZ, velX, velY, velZ\n";
    
    for (list<HopGroup>::iterator it=blobHaloes.begin(); it!=blobHaloes.end(); ++it)
    {
        outputFile << it->id << ", " << it->noParticles << ", " << it->mass;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << ", " << it->pos[i]/MPC_UNIT;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << ", " << it->vel[i];
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}
void outputASCII_Haloes(HopGroup *blobHaloes,
                        Int const noGroups,
                        string outputFileName)
{
    cout << "Writing halo information to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    outputFile << "id, noParticles, mass, posX, posY, posZ, velX, velY, velZ\n";
    
    for (Int j=0; j<noGroups; ++j)
    {
        outputFile << blobHaloes[j].id << ", " << blobHaloes[j].noParticles << ", " << blobHaloes[j].mass;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << ", " << blobHaloes[j].pos[i]/MPC_UNIT;
        for (int i=0; i<NO_DIM; ++i)
            outputFile << ", " << blobHaloes[j].vel[i];
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}




/* Outputs the haloes position to an ASCII file. */
void outputASCII_HaloesPositions(list<HopGroup> &blobHaloes,
                                 string outputFileName)
{
    cout << "Writing halo positions to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    for (list<HopGroup>::iterator it=blobHaloes.begin(); it!=blobHaloes.end(); ++it)
    {
        outputFile << it->pos[0]/MPC_UNIT;
        for (int i=1; i<NO_DIM; ++i)
            outputFile << "\t" << it->pos[i]/MPC_UNIT;
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}
void outputASCII_HaloesPositions(HopGroup *blobHaloes,
                                 Int const noGroups,
                                 string outputFileName)
{
    cout << "Writing halo positions to the ASCII file '" << outputFileName << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFileName );   //open the output file
    
    for (Int j=0; j<noGroups; ++j)
    {
        outputFile << blobHaloes[j].pos[0]/MPC_UNIT;
        for (int i=1; i<NO_DIM; ++i)
            outputFile << "\t" << blobHaloes[j].pos[i]/MPC_UNIT;
        outputFile << "\n";
    }
    outputFile.close();
    cout << "Done.\n";
}


