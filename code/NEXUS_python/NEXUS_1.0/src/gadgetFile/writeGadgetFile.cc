#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;


#include "gadgetFile.h"
#include <miscellaneous.h>







/* This function writes a Gadget2 snapshot.
*/
void writeGadgetFile(string &outputFileName,
		     Gadget_header &gadgetHeader,
		     Particle p[])
{
#define WRITE_BUFFER outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) )
	
    string const del(3, static_cast<char>( 8) );	// string that deletes 3 characters from output
    cout << "Writing the Gadget snapshot file. Done:  0%" << flush;
    int buffer;	//variable storing how many bytes has each block of data
    
    // Open the file for writting
    gadgetHeader.num_files = 1;
    fstream outputFile;
    openOutputBinaryFile( outputFile, outputFileName );
    
    // write the header
    buffer = 256;	//Gadget2 specifically tests that buffer=256 when reading the header
    WRITE_BUFFER;
    outputFile.write( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
    WRITE_BUFFER;
    
    // write the position of the particles
    int noParticles = gadgetHeader.totalParticles();
    buffer = noParticles*3*sizeof(float);	//number particles * 3 floats * 4 bytes/float
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
    WRITE_BUFFER;
    cout << del << setw(2) << int(3./7.*100.) << "%" << flush;
    
    // write the velocity of the particles
    buffer = noParticles*3*sizeof(float);
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
    WRITE_BUFFER;
    cout << del << setw(2) << int(6./7.*100.) << "%" << flush;
    
    // write the tags of the particles
    buffer = noParticles*sizeof(int);
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(&p[i].id), sizeof(int) );
    WRITE_BUFFER;
    
    
    // write masses if different
    if ( gadgetHeader.differentMass() )
    {
        noParticles = gadgetHeader.totalParticlesDifferentMass();
        buffer = noParticles*sizeof(float);
        WRITE_BUFFER;
    
    // now write the masses of particles with different masses
        int temp = 0;
        for (int i=0; i<6; ++i)
        {
            if ( gadgetHeader.differentMass(i) )
                for (int j=temp; j<temp+gadgetHeader.npart[i]; ++j)
                    outputFile.write( reinterpret_cast<char *>(&p[j].mass), sizeof(float) );
            temp += gadgetHeader.npartTotal[i];
        }
    
        WRITE_BUFFER;
    }
    
    outputFile.close();
    cout << del << "100%\n";
}


/* This function writes only the particle positions in a truncated Gadget2 snapshot file (file contains only the Gadget header and particle positions).
*/
void writeGadgetFile(string &outputFileName,
                     Gadget_header &gadgetHeader,
                     Particle_pm p[])
{
    cout << "Writing the particle positions in a truncated Gadget snapshot file (contains only particle positions)... " << flush;
    int buffer;	//variable storing how many bytes has each block of data
	
	// Open the file for writting
    fstream outputFile;
    openOutputBinaryFile( outputFile, outputFileName );
	
	// write the header
    buffer = 256;	//Gadget2 specifically tests that buffer=256 when reading the header
    WRITE_BUFFER;
    outputFile.write( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
    WRITE_BUFFER;
	
	// write the position of the particles
    int noParticles = gadgetHeader.totalParticles();
    buffer = noParticles*3*sizeof(float);	//number particles * 3 floats * 4 bytes/float
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
    WRITE_BUFFER;
	
    outputFile.close();
    cout << "Done\n" << flush;
}



/* This function writes in a Gadget2 snapshot only the dark matter (DM) particles with position, velocity, mass and tag.
*/
void writeGadgetFile_DM(string &outputFileName,
		        Gadget_header gadgetHeader,
                        Particle p[])
{
#define WRITE_BUFFER outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) )
	
            string const del(3, static_cast<char>( 8) );	// string that deletes 3 characters from output
    cout << "Writing the Gadget snapshot file. Done:  0%" << flush;
    int buffer;	//variable storing how many bytes has each block of data
    
    // Open the file for writting
    fstream outputFile;
    openOutputBinaryFile( outputFile, outputFileName );
    
    // Change header to reflect that there are only DM particles
    gadgetHeader.darkMatterOnly();
    
    // write the header
    buffer = 256;	//Gadget2 specifically tests that buffer=256 when reading the header
    WRITE_BUFFER;
    outputFile.write( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
    WRITE_BUFFER;
    
    // write the position of the particles
    int noParticles = gadgetHeader.npartTotal[1];
    buffer = noParticles*3*sizeof(float);	//number particles * 3 floats * 4 bytes/float
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
    WRITE_BUFFER;
    cout << del << setw(2) << int(3./7.*100.) << "%" << flush;
    
    // write the velocity of the particles
    buffer = noParticles*3*sizeof(float);
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
    WRITE_BUFFER;
    cout << del << setw(2) << int(6./7.*100.) << "%" << flush;
    
    // write the tags of the particles
    buffer = noParticles*sizeof(int);
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(&p[i].id), sizeof(int) );
    WRITE_BUFFER;
    
    
    // write masses if different
    if ( gadgetHeader.differentMass(1) )
    {
        buffer = noParticles*sizeof(float);
    
        WRITE_BUFFER;
        for (int i=0; i<noParticles; ++i)
            outputFile.write( reinterpret_cast<char *>(&p[i].id), sizeof(float) );
        WRITE_BUFFER;
    }
    
    outputFile.close();
    cout << del << "100%\n";
}


/* This function writes only the dark matter (DM) particle positions in a truncated Gadget2 snapshot file (file contains only the Gadget header and particle positions).
*/
void writeGadgetFile_DM(string &outputFileName,
                        Gadget_header gadgetHeader,
                        Particle_pm p[])
{
    cout << "Writing the dark matter (DM) particle positions in a truncated Gadget snapshot file (contains only DM particle positions)... " << flush;
    int buffer;	//variable storing how many bytes has each block of data
	
	// Open the file for writting
    fstream outputFile;
    openOutputBinaryFile( outputFile, outputFileName );
	
	// Change header to reflect that there are only DM particles
    gadgetHeader.darkMatterOnly();
	
	// write the header
    buffer = 256;	//Gadget2 specifically tests that buffer=256 when reading the header
    WRITE_BUFFER;
    outputFile.write( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
    WRITE_BUFFER;
	
	// write the position of the particles
    int noParticles = gadgetHeader.npartTotal[1];
    buffer = noParticles*3*sizeof(float);	//number particles * 3 floats * 4 bytes/float
    WRITE_BUFFER;
    for (int i=0; i<noParticles; ++i)
        outputFile.write( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
    WRITE_BUFFER;
	
    outputFile.close();
    cout << "Done\n" << flush;
}
