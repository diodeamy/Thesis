#ifndef DENSITY_FILE_HEADER
#define DENSITY_FILE_HEADER


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;

#include <defines.h>
#include <density_header.h>
#include <miscellaneous.h>
#include <array.h>


//! The following section contains functions for reading a density file (i.e. a binary file with the header "Density_header" followed by binary data for each grid cell - the two blocks are preceded and followed by 'size_t' blocks giving the size of the binary data written in each block)


/* This function reads the density header of a density file. */
template <typename T_Density_header>
void readDensityHeader(T_Density_header *densityHeader,
                       std::string &inputFileName,
                       bool verbose = true)
{
    if ( verbose )
        std::cout << "\nReading the header of the density file '" << inputFileName << "' ... " << std::flush;
    
    // check if the data is in one file or in multiple ones
    std::string fileName = inputFileName;
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName = densityHeader->filename( fileName, 0 );
        if ( not bfs::exists(fileName) )
            throwError( "The program could not open the input density file/files. It cannot find the file/files." );
    }
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );
    size_t buffer1, buffer2;
    
    // read the header of the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=sizeof(*densityHeader) )
        throwError( "There was an error while reading the header of the density file in function 'readDensityHeader'. The size of the binary data (given by 'buffer1') does not match with the size of the density header." );
    inputFile.read( reinterpret_cast<char *>(densityHeader), sizeof(*densityHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error while reading the header of the density file in function 'readDensityHeader'. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    if ( densityHeader->totalGrid!=int(densityHeader->gridSize[0]*densityHeader->gridSize[1]*densityHeader->gridSize[2]) )
        throwError( "There was an error while reading the header of the density file. The total number of grid points does not match with the product of the grid points along each direction." );
    
    inputFile.close();
    if ( verbose )
        std::cout << "Done.\n";
}

#include <typeinfo>
/* This function reads the data from one single binary density file. */
template <typename T>
void readDensityData_singleFile(T *density,
                                size_t gridSize,
                                std::string &inputFileName)
{
    std::cout << "Reading the data in the binary density file '" << inputFileName << "' ... " << std::flush;
    
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, inputFileName );
    size_t buffer1, buffer2;
    
    // skip over the header part
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.seekg( buffer1+sizeof(buffer2), std::ios::cur );
    
    // read the data in the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=gridSize*sizeof(density[0]) )
        throwError( "There was an error while reading the data from the density file in function 'readDensityData_singleFile'. The length of the binary data (given by 'buffer1') does not match with the number of grid cells in the density header." );
    inputFile.read( reinterpret_cast<char *>(density), buffer1 );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error while reading the data from the density file in function 'readDensityData_singleFile'. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    inputFile.close();
    std::cout << "Done.\n";
}



/* Reads the density data splitted in multiple binary density files. */
template <typename T>
void readDensityData_multipleFile(Density_header &densityHeader,
                                  T *density,
                                  size_t gridSize,
                                  std::string &rootFilename)
{
    std::fstream inputFile;
    size_t buffer1, buffer2;
    
    if ( gridSize!=densityHeader.dataSize() )
        throwError( "The size of the array where to read the density is different than the total size of the density grid." );
    std::cout << "Reading the density data saved within multiple files. The density grid is saved in " << densityHeader.noDensityFiles << " binary files:\n";
    
    
    for (int i=0; i<densityHeader.noDensityFiles; ++i)
    {
        std::string fileName = densityHeader.filename( rootFilename, i );
        std::cout << "\t reading file " << i << ": '" << fileName << "' ... " << std::flush;
        openInputBinaryFile( inputFile, fileName );
        Density_header tempHeader;
        
        // read the header of the file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 or buffer1!=sizeof(tempHeader) )
            throwError( "There was an error while reading the header of the density file. 'buffer1' and 'buffer2' do not have the same value or their values are different from the size of the density header. The file may be corrupt." );
        if ( i!=tempHeader.indexDensityFile )
            throwError( "The index of the density file does not match with the expected index in function 'readDensityDataMultipleFiles'.");
        
        // read the data in the file
        size_t indices[6];
        tempHeader.boundariesIndices( i, indices );
        size_t const tempGridSize = (indices[1]-indices[0]) * (indices[3]-indices[2]) * (indices[5]-indices[4]);
        size_t const tempSize[] = { densityHeader.gridSize[1]*densityHeader.gridSize[2], densityHeader.gridSize[2] };
        
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        if ( buffer1!=size_t(tempGridSize*sizeof(density[0])) )
            throwError( "There was an error while reading the data from the density file. The length of the binary data (given by 'buffer1') does not match with the number of grid cells in the density header." );
        for (size_t i1=indices[0]; i1<indices[1]; ++i1)
            for (size_t i2=indices[2]; i2<indices[3]; ++i2)
            {
                size_t const index = i1 * tempSize[0] + i2 * tempSize[1] + indices[4];
                inputFile.read( reinterpret_cast<char *>(&(density[index])), sizeof(density[0])*(indices[5]-indices[4]) );
            }
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( buffer1!=buffer2 )
            throwError( "There was an error while reading the data in the binary density file. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
        
        inputFile.close();
        std::cout << "Done.\n";
    }
}


/* This function reads the data stored in a binary density file. */
template <typename T>
void readDensityData(T *density,
                     size_t gridSize,
                     std::string &inputFileName)
{
    Density_header densityHeader;
    readDensityHeader( &densityHeader, inputFileName, false );
    
    
    // check in how many files the data is written and ask for the corresponding function
    if ( densityHeader.noDensityFiles>1 )
        readDensityData_multipleFile( densityHeader, density, gridSize, inputFileName );
    else
        readDensityData_singleFile( density, gridSize, inputFileName );
}
/* Does the same as the previous function but for an "Array<T,N>" element.*/
template <typename T, size_t N>
void readDensityData(Array<T,N> *density,
                     std::string &inputFileName)
{
    size_t gridSize = density->size();
    readDensityData<T>( density->ptrData(), gridSize, inputFileName );
}





//! The following functions write the data into binary density files.


/* This function writes the data to a single binary density file. */
template <typename T>
void writeDensityData_singleFile(T *density,
                                size_t const gridSize,
                                Density_header &densityHeader,
                                std::string outputFileName)
{
    std::fstream outputFile;	// the output file
    openOutputBinaryFile( outputFile, outputFileName );   //open the output file
    
    std::cout << "Writing the data to the binary file '" << outputFileName << "' ...  " << std::flush;
    size_t buffer;	//variable storing how many bytes has each block of data
    if ( gridSize!=densityHeader.dataSize() and densityHeader.noDensityFiles==1 )
        throwError( "The density data to be written to file in function 'writeDensityData' has a different size than the one indicated in the density header supplied as an argument to the function." );
    else if ( gridSize!=densityHeader.secondaryTotalSize( densityHeader.indexDensityFile ) )
        throwError( "You selected to write the density in multiple binary files. But the density data to be written in one file in function 'writeDensityData' has a different size than the one indicated in the density header supplied as an argument to the function." );
    
    // write the header
    buffer = sizeof(densityHeader);
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&densityHeader), sizeof(densityHeader) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the position of the particles
    buffer = gridSize*sizeof(density[0]);	//total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(density), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    outputFile.close();
    std::cout << "Done.\n";
}


/* This function writes the data to multiple binary density files. */
template <typename T>
void writeDensityData_multipleFile(T *density,
                                   size_t const gridSize,
                                   Density_header &densityHeader,
                                   std::string rootFilename)
{
    std::cout << "Writing the data into " << densityHeader.noDensityFiles << " binary density files:\n";
    if ( gridSize!=densityHeader.dataSize() )
        throwError( "The data to be written to multiple files has a different size than the one indicated in the density header supplied as an argument to the function." );
    Density_header tempHeader = densityHeader;
    
    for (int i=0; i<densityHeader.noDensityFiles; ++i)
    {
        std::string fileName = densityHeader.filename( rootFilename, i );
        std::cout << "\t writting file " << i << ": '" << fileName << "' ... " << std::flush;
        std::fstream outputFile;	// the output file
        openOutputBinaryFile( outputFile, fileName );
        
        //write the header
        tempHeader.indexDensityFile = i;
        size_t buffer = sizeof( tempHeader );
        outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
        outputFile.write( reinterpret_cast<char *>(&densityHeader), sizeof(densityHeader) );
        outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
        
        //write the data
        size_t indices[6];
        tempHeader.boundariesIndices( i, indices );
        size_t const tempGridSize = (indices[1]-indices[0]) * (indices[3]-indices[2]) * (indices[5]-indices[4]);
        size_t const tempSize[] = { densityHeader.gridSize[1]*densityHeader.gridSize[2], densityHeader.gridSize[2] };
        
        buffer = tempGridSize * sizeof(density[0]);
        outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
        for (size_t i1=indices[0]; i1<indices[1]; ++i1)
            for (size_t i2=indices[2]; i2<indices[3]; ++i2)
            {
                size_t const index = i1 * tempSize[0] + i2 * tempSize[1] + indices[4];
                outputFile.write( reinterpret_cast<char *>(&(density[index])), sizeof(density[0])*(indices[5]-indices[4]) );
            }
        outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
        
        outputFile.close();
        std::cout << "Done.\n";
    }
}


/* This function writes the data to a binary density file.
 NOTE: If "fileIndex>=0 or fileIndex<densityHeader.noDensityFiles" than it means that the "density" data corresponds to file 'fileIndex' from a multiple binary density files batch. */
template <typename T>
void writeDensityData(T *density,
                      size_t const gridSize,
                      Density_header &densityHeader,
                      std::string outputFileName,
                      int fileIndex = -1)
{
    if ( densityHeader.noDensityFiles>1 and not (fileIndex>=0 or fileIndex<densityHeader.noDensityFiles) )
        writeDensityData_multipleFile( density, gridSize, densityHeader, outputFileName );
    else
        writeDensityData_singleFile( density, gridSize, densityHeader, outputFileName );
}
/* Does the same as the previous function but for an "Array<T,N>" element.*/
template <typename T, size_t N>
void writeDensityData(Array<T,N> &density,
                      Density_header &densityHeader,
                      std::string outputFileName,
                      int fileIndex = -1)
{
    size_t const gridSize = density.size();
    writeDensityData<T>( density.ptrData(), gridSize, densityHeader, outputFileName, fileIndex );
}


#endif




