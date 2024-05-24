#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;


#include "gadgetFile.h"
#include <gadget_header.h>
#include <miscellaneous.h>




/* This function reads the header of a GADGET2 snapshot file.
*/
void readGadgetHeader(string const &inputFileName,
                      Gadget_header *gadgetHeader)
{
    cout << "\nReading the header of the Gadget snapshot file '" << inputFileName << "' ... " << flush;
    // check to see if there is only one input file or several
    string fileName = inputFileName;
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName += ".0";
        if ( not bfs::exists(fileName) )
            throwError( "The program could not open the input GADGET snapshot file/files. It cannot find the file/files." );
    }
    
    // open the file and read the Gadget header
    fstream inputFile;
    openInputBinaryFile( inputFile, fileName );
    int buffer1, buffer2;
    
    // read the header of the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(gadgetHeader), sizeof(*gadgetHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    inputFile.close();
    
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET snapshot file. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
    
    // if the GADGET2 snapshot is in multiple files, check to see if all are present
    if ( gadgetHeader->num_files!=1 )
        for (int i=0; i<gadgetHeader->num_files; ++i)
        {
            ostringstream buffer;
            buffer << inputFileName << "." << i;
            if ( not bfs::exists( buffer.str() ) )
                throwError( "Missing GADGET snapshot file. The program could not find the file '" + buffer.str() + "'." );
        }
    
    cout << "Done.\n";
}




/* This function reads the position, velocity, identity and masses of GAS and HALO particles, storing them in the same array 'p'.
The first 'gadgetHeader.npartTotal[0]' elements are the GAS particles while the remaining 'gadgetHeader.npartTotal[1]' elements are the HALO particles.
*/
void readGadgetData(Particle p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName)
{
    fstream inputFile;        // the input file
    string fileName;          // will store the name of the input file
    Gadget_header tempHeader; //temporary header that stores the header of each new file as it is read
    int buffer1, buffer2;     // int variable used for reading the integers before each data block and for error detection
    int noReadGas = 0;        // variable that keeps track of how many GAS particles were read
    int noReadHalo = 0;       // variable that keeps track of how many HALO particles were read
    
    // iterate over all the files and read the data in each of them
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        // get the file name
        if ( gadgetHeader.num_files==1 )
            fileName = inputFileName;
        else
        {
            ostringstream buffer;
                buffer << inputFileName << "." << i;
            fileName = buffer.str();
        }
        cout << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n";
        
        // open the file and read the header
        openInputBinaryFile( inputFile, fileName );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
        
        // read the position block
        cout << "\t reading positions of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // read the velocity block
        cout << "\t reading velocities of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle velocity data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // read the identity block
        cout << "\t reading identities of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(&p[i].id), sizeof(int) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(&p[i].id), sizeof(int) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle identity data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // read masses if different
        cout << "\t reading masses of the particles... " << flush;
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
            inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        if ( tempHeader.mass[0]==0. )
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                p[i].mass = tempHeader.mass[0];
        if ( tempHeader.mass[1]==0. )
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                p[i].mass = tempHeader.mass[1];
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
        {
            inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
            if ( buffer1!=buffer2 )
                throwError( "The integers before and after the particle mass data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        }
        cout << "Done\n";
        
        
        inputFile.close();
        
        noReadGas += tempHeader.npart[0];
        noReadHalo += tempHeader.npart[1];
    }
       
       // check that indeed we read all the particles
    if ( noParticles != noReadGas+noReadHalo )
    {
        ostringstream buffer;
        buffer << "Reading the particles in the GADGET snapshot. Expected to read " << noParticles << " particles, but the code read " << noReadGas+noReadHalo << ".";
        throwError( buffer.str() );
    }
}



/* This function reads the position, velocity and masses of GAS and HALO particles, storing them in the same array 'p'.
The first 'gadgetHeader.npartTotal[0]' elements are the GAS particles while the remaining 'gadgetHeader.npartTotal[1]' elements are the HALO particles.
*/
void readGadgetData(Particle_pvm p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName)
{
    fstream inputFile;        // the input file
    string fileName;          // will store the name of the input file
    Gadget_header tempHeader; //temporary header that stores the header of each new file as it is read
    int buffer1, buffer2;     // int variable used for reading the integers before each data block and for error detection
    int noReadGas = 0;        // variable that keeps track of how many GAS particles were read
    int noReadHalo = 0;       // variable that keeps track of how many HALO particles were read
    
    // iterate over all the files and read the data in each of them
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        // get the file name
        if ( gadgetHeader.num_files==1 )
            fileName = inputFileName;
        else
        {
            ostringstream buffer;
            buffer << inputFileName << "." << i;
            fileName = buffer.str();
        }
        cout << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n";
        
        // open the file and read the header
        openInputBinaryFile( inputFile, fileName );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
        
        // read the position block
        cout << "\t reading positions of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // read the velocity block
        cout << "\t reading velocities of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].vel), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle velocity data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // jump over particle identity block
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.seekg( buffer1+sizeof(int), ios::cur );
        
        // read masses if different, otherwise assign header value
        cout << "\t reading masses of the particles... " << flush;
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
            inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        if ( tempHeader.mass[0]==0. )
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                p[i].mass = tempHeader.mass[0];
        if ( tempHeader.mass[1]==0. )
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                p[i].mass = tempHeader.mass[1];
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
        {
            inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
            if ( buffer1!=buffer2 )
                throwError( "The integers before and after the particle mass data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        }
        cout << "Done\n";
        
        
        inputFile.close();
        
        noReadGas += tempHeader.npart[0];
        noReadHalo += tempHeader.npart[1];
    }
       
       // check that indeed we read all the particles
    if ( noParticles != noReadGas+noReadHalo )
    {
        ostringstream buffer;
        buffer << "Reading the particles in the GADGET snapshot. Expected to read " << noParticles << " particles, but the code read " << noReadGas+noReadHalo << ".";
        throwError( buffer.str() );
    }
}



/* This function reads the position and masses of GAS and HALO particles, storing them in the same array 'p'.
The first 'gadgetHeader.npartTotal[0]' elements are the GAS particles while the remaining 'gadgetHeader.npartTotal[1]' elements are the HALO particles.
*/
void readGadgetData(Particle_pm p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName)
{
    fstream inputFile;        // the input file
    string fileName;          // will store the name of the input file
    Gadget_header tempHeader; //temporary header that stores the header of each new file as it is read
    int buffer1, buffer2;     // int variable used for reading the integers before each data block and for error detection
    int noReadGas = 0;        // variable that keeps track of how many GAS particles were read
    int noReadHalo = 0;       // variable that keeps track of how many HALO particles were read
    
    // iterate over all the files and read the data in each of them
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        // get the file name
        if ( gadgetHeader.num_files==1 )
            fileName = inputFileName;
        else
        {
            ostringstream buffer;
            buffer << inputFileName << "." << i;
            fileName = buffer.str();
        }
        cout << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n";
        
        // open the file and read the header
        openInputBinaryFile( inputFile, fileName );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
        
        // read the position block
        cout << "\t reading positions of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        // jump over velocity block
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.seekg( buffer1+sizeof(int), ios::cur );
        
        // jump over particle identity block
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.seekg( buffer1+sizeof(int), ios::cur );
        
        // read masses if different, otherwise assign header value
        cout << "\t reading masses of the particles... " << flush;
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
            inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        if ( tempHeader.mass[0]==0. )
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
                p[i].mass = tempHeader.mass[0];
        if ( tempHeader.mass[1]==0. )
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                inputFile.read( reinterpret_cast<char *>(&p[i].mass), sizeof(float) );
        else
            for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
                p[i].mass = tempHeader.mass[1];
        if ( (tempHeader.mass[0]==0. and tempHeader.npart[0]!=0) or (tempHeader.mass[1]==0. and tempHeader.npart[1]!=0) )
        {
            inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
            if ( buffer1!=buffer2 )
                throwError( "The integers before and after the particle mass data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        }
        cout << "Done\n";
        
        
        inputFile.close();
        
        noReadGas += tempHeader.npart[0];
        noReadHalo += tempHeader.npart[1];
    }
       
       // check that indeed we read all the particles
    if ( noParticles != noReadGas+noReadHalo )
    {
        ostringstream buffer;
        buffer << "Reading the particles in the GADGET snapshot. Expected to read " << noParticles << " particles, but the code read " << noReadGas+noReadHalo << ".";
        throwError( buffer.str() );
    }
}




/* This function reads the position of GAS and HALO particles, storing them in the same array 'p'.
The first 'gadgetHeader.npartTotal[0]' elements are the GAS particles while the remaining 'gadgetHeader.npartTotal[1]' elements are the HALO particles.
*/
void readGadgetData(Particle_p p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName)
{
    fstream inputFile;        // the input file
    string fileName;          // will store the name of the input file
    Gadget_header tempHeader; //temporary header that stores the header of each new file as it is read
    int buffer1, buffer2;     // int variable used for reading the integers before each data block and for error detection
    int noReadGas = 0;        // variable that keeps track of how many GAS particles were read
    int noReadHalo = 0;       // variable that keeps track of how many HALO particles were read
    
    // iterate over all the files and read the data in each of them
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        // get the file name
        if ( gadgetHeader.num_files==1 )
            fileName = inputFileName;
        else
        {
            ostringstream buffer;
            buffer << inputFileName << "." << i;
            fileName = buffer.str();
        }
        cout << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n";
        
        // open the file and read the header
        openInputBinaryFile( inputFile, fileName );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
        
        // read the position block
        cout << "\t reading positions of the particles... " << flush;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        for (int i=noReadGas; i<noReadGas+tempHeader.npart[0]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        for (int i=gadgetHeader.npartTotal[0]+noReadHalo; i<gadgetHeader.npartTotal[0]+noReadHalo+tempHeader.npart[1]; ++i)
            inputFile.read( reinterpret_cast<char *>(p[i].pos), 3*sizeof(float) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." );
        cout << "Done\n";
        
        
        inputFile.close();
        
        noReadGas += tempHeader.npart[0];
        noReadHalo += tempHeader.npart[1];
    }
       
       // check that indeed we read all the particles
    if ( noParticles != noReadGas+noReadHalo )
    {
        ostringstream buffer;
        buffer << "Reading the particles in the GADGET snapshot. Expected to read " << noParticles << " particles, but the code read " << noReadGas+noReadHalo << ".";
        throwError( buffer.str() );
    }
}



