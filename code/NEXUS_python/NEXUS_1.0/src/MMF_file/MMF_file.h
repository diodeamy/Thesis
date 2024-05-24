#ifndef MMF_FILE_HEADER
#define MMF_FILE_HEADER

#include <string>
#include <fstream>

#include <defines.h>
#include "MMF_header.h"
#include <miscellaneous.h>
#include <array.h>
#include <vector.h>


typedef Array<Real,NO_DIM>             ArrayReal3D;
typedef Vector<Real,NO_DIM>            VectorReal3D;
typedef Vector<Real,NO_DIM*(NO_DIM-1)> VectorReal9D;
typedef Array<VectorReal3D,NO_DIM>     ArrayVector3D;
typedef Array<VectorReal9D,NO_DIM>     ArrayVector9D;




//! Functions used to read the data in MMF files

/* This function reads the MMF header of a MMF file. */
template <typename T_MMF_header>
void readMMF_Header(T_MMF_header *mmfHeader,
                    std::string &inputFileName,
                    bool verbose = true)
{
    if ( verbose )
        std::cout << "\nReading the header of the MMF file '" << inputFileName << "' ... " << std::flush;
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, inputFileName );
    size_t buffer1, buffer2;
    
    // read the header of the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(mmfHeader), sizeof(*mmfHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    
    if ( buffer1!=buffer2 or buffer1!=sizeof(*mmfHeader) )
        throwError( "While reading the header of the MMF file. 'buffer1='", buffer1, " while 'buffer2'=", buffer2, " when both should be 1024." );
    if ( mmfHeader->totalGrid!=(mmfHeader->gridSize[0]*mmfHeader->gridSize[1]*mmfHeader->gridSize[2]) )
        throwError( "While reading the header of the MMF file. The total number of grid points does not match with the product of the grid points along each direction." );
    
    inputFile.close();
    if ( verbose )
        std::cout << "Done.\n";
}

template <typename T_MMF_header>
void readMMF_Header(T_MMF_header *mmfHeader,
                    std::fstream &inputFile)
{
    // open the file
    size_t buffer1, buffer2;
    
    // read the header of the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(mmfHeader), sizeof(*mmfHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    
    if ( buffer1!=buffer2 or buffer1!=sizeof(*mmfHeader) )
        throwError( "While reading the header of the MMF file. 'buffer1='", buffer1, " while 'buffer2'=", buffer2, " when both should be 1024." );
    if ( mmfHeader->totalGrid!=(mmfHeader->gridSize[0]*mmfHeader->gridSize[1]*mmfHeader->gridSize[2]) )
        throwError( "While reading the header of the MMF file. The total number of grid points does not match with the product of the grid points along each direction." );
}



/* This function is used to read the data from an MMF file. */
template <typename T>
void readMMF_file(T *data,
                  std::string &inputFileName,
                  size_t skipSize = 1040)
{
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, inputFileName );
    size_t buffer1, buffer2;
    
    //skip the first part of the file
    inputFile.seekg( skipSize, std::ios::cur );
    
    //read the data
    size_t expectedBuffer = data->size() * sizeof((*data)[0]);
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer before the data block does not agree with the expected size of the MMF data." );
    inputFile.read( reinterpret_cast<char *>(data->ptrData()), expectedBuffer );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer2!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer after the data block does not agree with the expected size of the MMF data." );
    
    inputFile.close();
}
template <typename T1, typename T2>
void readMMF_file(T1 *data1,
                  T2 *data2,
                  std::string &inputFileName,
                  size_t skipSize = 1040)
{
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, inputFileName );
    size_t buffer1, buffer2;
    
    //skip the first part of the file
    inputFile.seekg( skipSize, std::ios::cur );
    
    //read the 1st data
    size_t expectedBuffer = data1->size() * sizeof((*data1)[0]);
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer before the data block does not agree with the expected size of the MMF data." );
    inputFile.read( reinterpret_cast<char *>(data1->ptrData()), expectedBuffer );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer2!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer after the data block does not agree with the expected size of the MMF data." );
    
    //read the 2nd data
    expectedBuffer = data2->size() * sizeof((*data2)[0]);
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer before the data block does not agree with the expected size of the MMF data." );
    inputFile.read( reinterpret_cast<char *>(data2->ptrData()), expectedBuffer );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer2!=expectedBuffer )
        throwError( "There was an error while reading the MMF binary file. The value of the buffer after the data block does not agree with the expected size of the MMF data." );
    
    inputFile.close();
}



/* Functions to read the different components and variables of the MMF data. */
template <typename T>
void readMMF_response(T *response,
                      std::string &inputFileName)
{
    std::cout << "Reading the MMF response from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_RESPONSE and mmfHeader.fileType!=MMF_MAX_RESPONSE )
        throwError( "Your are trying to read the MMF response data from a file that does not store the MMF response. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( response, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T1, typename T2>
void readMMF_maxResponse(T1 *response,
                         T2 *scale,
                         std::string &inputFileName)
{
    std::cout << "Reading the MMF maximum response and scale from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_MAX_RESPONSE_SCALE  )
        throwError( "Your are trying to read the MMF maximum response and scale data from a file that does not store the MMF maximum response and its corresponding scale. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( response, scale, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T>
void readMMF_scale(T *scale,
                   std::string &inputFileName)
{
    std::cout << "Reading the MMF scale from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_MAX_RESPONSE_SCALE )
        throwError( "Your are trying to read the scale data for the maximum MMF response from a file that does not store the scale of the maximum MMF response. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader) + 2*sizeof(size_t) + mmfHeader.totalGrid*sizeof(Real);
    readMMF_file( scale, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T1, typename T2>
void readMMF_eigen(T1 *eigenvalues,
                   T2 *eigenvectors,
                   std::string &inputFileName)
{
    std::cout << "Reading the MMF eigenvalues and eigenvectors from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_EIGEN and mmfHeader.fileType!=MMF_MAX_EIGEN )
        throwError( "Your are trying to read the MMF eigenvalues and eigenvectors data from a file that does not store this type of data. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( eigenvalues, eigenvectors, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T>
void readMMF_eigenvalues(T *eigenvalues,
                         std::string &inputFileName)
{
    std::cout << "Reading the MMF eigenvalues from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_EIGEN and mmfHeader.fileType!=MMF_MAX_EIGEN )
        throwError( "Your are trying to read the MMF eigenvalues data from a file that does not store this type of data. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( eigenvalues, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T>
void readMMF_eigenvectors(T *eigenvectors,
                          std::string &inputFileName)
{
    std::cout << "Reading the MMF eigenvectors from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_EIGEN and mmfHeader.fileType!=MMF_MAX_EIGEN )
        throwError( "Your are trying to read the MMF eigenvectors data from a file that does not store this type of data. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader) + 2*sizeof(size_t) + mmfHeader.totalGrid*NO_DIM*sizeof(Real);
    readMMF_file( eigenvectors, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T>
void readMMF_cleanResponse(T *cleanResponse,
                           std::string &inputFileName)
{
    std::cout << "Reading the clean MMF response from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    if ( mmfHeader.fileType!=MMF_CLEAN_RESPONSE )
        throwError( "Your are trying to read the clean MMF response data from a file that does not store this type of data. Curent file type is '", mmfHeader.fileType, "'." );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( cleanResponse, inputFileName, skipSize );
    std::cout << "Done.\n";
}
template <typename T>
void readMMF(T *data,
             std::string &inputFileName,
             std::string variableName = "data")
{
    std::cout << "Reading the " << variableName << " from the binary MMF file '" << inputFileName << "' ... " << std::flush;
    MMF_header mmfHeader;
    readMMF_Header( &mmfHeader, inputFileName, false );
    size_t skipSize = 2*sizeof(size_t) + sizeof(mmfHeader);
    readMMF_file( data, inputFileName, skipSize );
    std::cout << "Done.\n";
}







//! Functions used to write data to a MMF files

/* Writes a general MMF binary file. */
template <typename T>
void writeMMF_file(T &data,
                   MMF_header &mmfHeader,
                   std::string outputFileName)
{
    std::fstream outputFile;    // the output file
    openOutputBinaryFile( outputFile, outputFileName );   //open the output file
    
    if ( mmfHeader.totalGrid!=data.size() )
        throwError( "In function 'writeMMF_file' the MMF header and the data to be writen do not have the same size." );
    
    // write the header
    size_t buffer;  //variable storing how many bytes has each block of data
    buffer = sizeof(mmfHeader); //specific test that buffer=1024 when reading the header
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&mmfHeader), sizeof(mmfHeader) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the 1st data
    size_t gridSize = mmfHeader.totalSize();
    buffer = gridSize * sizeof(data[0]);    //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(data.ptrData()), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    outputFile.close();
}
template <typename T1, typename T2>
void writeMMF_file(T1 &data1,
                   T2 &data2,
                   MMF_header &mmfHeader,
                   std::string outputFileName)
{
    std::fstream outputFile;    // the output file
    openOutputBinaryFile( outputFile, outputFileName );   //open the output file
    
    if ( mmfHeader.totalGrid!=data1.size() or mmfHeader.totalGrid!=data2.size() )
        throwError( "In function 'writeMMF_file' the MMF header and the data to be writen do not have the same size." );
    
    // write the header
    size_t buffer;  //variable storing how many bytes has each block of data
    buffer = sizeof(mmfHeader); //specific test that buffer=1024 when reading the header
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&mmfHeader), sizeof(mmfHeader) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the 1st data
    size_t gridSize = mmfHeader.totalSize();
    buffer = gridSize * sizeof(data1[0]);   //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(data1.ptrData()), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the 2nd data
    buffer = gridSize * sizeof(data2[0]);   //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(data2.ptrData()), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    outputFile.close();
}


/* Functions to write to a binary file the different components and variables of the MMF data. */
template <typename T>
void writeMMF_response(T &response,
                       MMF_header &mmfHeader,
                       std::string outputFileName)
{
    std::cout << "Writing the MMF response to the binary file '" << outputFileName << "' ...  " << std::flush;
    writeMMF_file( response, mmfHeader, outputFileName );
    std::cout << "Done.\n";
}
template <typename T1, typename T2>
void writeMMF_maxResponse(T1 &response,
                          T2 &scale,
                          MMF_header &mmfHeader,
                          std::string outputFileName)
{
    std::cout << "Writing the MMF maximum response to the binary file '" << outputFileName << "' ...  " << std::flush;
    writeMMF_file( response, scale, mmfHeader, outputFileName );
    std::cout << "Done.\n";
}
template <typename T1, typename T2>
void writeMMF_eigen(T1 &eigenvalues,
                    T2 &eigenvectors,
                    MMF_header &mmfHeader,
                    std::string outputFileName)
{
    std::cout << "Writing the MMF eigenvalues and eigenvectors to the binary file '" << outputFileName << "' ...  " << std::flush;
    writeMMF_file( eigenvalues, eigenvectors, mmfHeader, outputFileName );
    std::cout << "Done.\n";
}
template <typename T>
void writeMMF(T &data,
              MMF_header &mmfHeader,
              std::string outputFileName,
              std::string variableName = "data")
{
    std::cout << "Writing the " << variableName << " to the binary MMF file '" << outputFileName << "' ...  " << std::flush;
    writeMMF_file( data, mmfHeader, outputFileName );
    std::cout << "Done.\n";
}


#endif
