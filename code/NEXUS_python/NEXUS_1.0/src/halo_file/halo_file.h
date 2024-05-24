#ifndef HALO_FILE_HEADER
#define HALO_FILE_HEADER


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;

#include <defines.h>
#include <halo_header.h>
#include <miscellaneous.h>
#include <array.h>


//! The following section contains functions for reading a halo file (i.e. a binary file with the header "Halo_header" followed by binary data for each grid cell - the two blocks are preceded and followed by 'size_t' blocks giving the size of the binary data written in each block)


/* This function reads the header of a halo binary file. */
template <typename T_Halo_header>
void readHaloHeader(T_Halo_header *haloHeader,
                    Array<char,2> *integerColumns,
                    Array<char,2> *floatColumns,
                    std::string &inputFileName,
                    bool verbose = true)
{
    if ( verbose )
        std::cout << "\nReading the header of the halo file '" << inputFileName << "' ... " << std::flush;
    
    // check if the data is in one file or in multiple ones
    std::string fileName = inputFileName;
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
        throwError( "The program could not open the input halo file/files. It cannot find the file/files." );
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );
    size_t buffer1, buffer2;
    
    // read the header of the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=sizeof(*haloHeader) )
        throwError( "There was an error while reading the header of the halo file in function 'readHaloHeader'. The size of the binary data (given by 'buffer1') does not match with the size of the halo header." );
    inputFile.read( reinterpret_cast<char *>(haloHeader), sizeof(*haloHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error while reading the header of the halo file in function 'readHaloHeader'. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    if ( haloHeader->noColumns!=haloHeader->noColumnsIntegers+haloHeader->noColumnsFloats )
        throwError( "There was an error while reading the header of the halo file. The total number of data columns is not the sum of the columns in the integer and floating point data sets." );
    
    // read the column names
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=haloHeader->noColumns*columnNameLength*sizeof(char) )
        throwError( "There was an error in function 'readHaloHeader' while reading the data block storing the column names of a halo binary file. The size of the binary data (given by 'buffer1') does not match with the size of the expected data." );
    inputFile.read( reinterpret_cast<char *>(integerColumns->ptrData()), haloHeader->noColumnsIntegers*columnNameLength*sizeof(char) );
    inputFile.read( reinterpret_cast<char *>(floatColumns->ptrData()), haloHeader->noColumnsFloats*columnNameLength*sizeof(char) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error in function 'readHaloHeader' while reading the data block storing the column names of a halo binary file. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    inputFile.close();
    if ( verbose )
        std::cout << "Done.\n";
}


/* This function reads the data from one binary halo file. */
template <typename T1, size_t N1, typename T2, size_t N2>
void readHaloData(Array<T1,N1> *dataIntegers,
                  Array<T2,N2> *dataFloats,
                  std::string &inputFileName,
                  bool verbose = true)
{
    if ( verbose )
        std::cout << "Reading the data in the binary halo file '" << inputFileName << "' ... " << std::flush;
    
    
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, inputFileName );
    size_t buffer1, buffer2;
    
    // skip over the header part
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.seekg( buffer1+sizeof(buffer2), std::ios::cur );
    
    // skip over the column name part
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.seekg( buffer1+sizeof(buffer2), std::ios::cur );
    
    // read the integer data in the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=dataIntegers->size()*sizeof((*dataIntegers)[0]) )
        throwError( "There was an error in function 'readHaloData' while reading the integer data from the halo binary file . The length of the binary data (given by 'buffer1') does not match with the size of the data expected from the halo header." );
    inputFile.read( reinterpret_cast<char *>(&((*dataIntegers)[0])), buffer1 );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error in function 'readHaloData' while reading the integer data from the binary halo file. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    // read the floating point data in the file
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1!=dataFloats->size()*sizeof((*dataFloats)[0]) )
        throwError( "There was an error in function 'readHaloData' while reading the floating point data from the halo binary file . The length of the binary data (given by 'buffer1') does not match with the size of the data expected from the halo header." );
    inputFile.read( reinterpret_cast<char *>(&((*dataFloats)[0])), buffer1 );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    if ( buffer1!=buffer2 )
        throwError( "There was an error in function 'readHaloData' while reading the floating point data from the binary halo file. 'buffer1' and 'buffer2' do not have the same value. The file may be corrupt." );
    
    inputFile.close();
    std::cout << "Done.\n";
}


/* This function writes the data to one binary halo file. */
template <typename T1, size_t N1, typename T2, size_t N2>
void writeHaloData(Halo_header &haloHeader,
                   Array<char,2> &integerColumns,
                   Array<char,2> &floatColumns,
                   Array<T1,N1> &dataIntegers,
                   Array<T2,N2> &dataFloats,
                   std::string &outputFileName,
                   bool verbose = true)
{
    if ( verbose )
        std::cout << "Writing the data to the binary halo file '" << outputFileName << "' ... " << std::flush;
    
    // do some error checking
    haloHeader.noHalos = dataIntegers.axisSize(0);
    if ( haloHeader.noHalos!=dataFloats.axisSize(0) )
        throwError( "There was an error in function 'writeHaloData'. The number of halos in the header does not match the number of halos in the floating point data set." );
    haloHeader.noColumnsIntegers = dataIntegers.axisSize(1);
    haloHeader.noColumnsFloats = dataFloats.axisSize(1);
    haloHeader.noColumns = haloHeader.noColumnsIntegers + haloHeader.noColumnsFloats;
    
    // find the mass range in the data
    haloHeader.massRange[0] = dataFloats(0,haloHeader.massColumn);
    haloHeader.massRange[1] = dataFloats(0,haloHeader.massColumn);
    for (Int i=1; i<haloHeader.noHalos; ++i)
    {
        if ( haloHeader.massRange[0]>dataFloats(i,haloHeader.massColumn) ) haloHeader.massRange[0] = dataFloats(i,haloHeader.massColumn);
        if ( haloHeader.massRange[1]<dataFloats(i,haloHeader.massColumn) ) haloHeader.massRange[1] = dataFloats(i,haloHeader.massColumn);
    }
    
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, outputFileName );
    
    // write the header part
    size_t buffer = sizeof(haloHeader);
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&haloHeader), sizeof(haloHeader) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the column names part
    buffer = haloHeader.noColumns * columnNameLength * sizeof(char);
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(integerColumns.ptrData()), haloHeader.noColumnsIntegers * columnNameLength * sizeof(char) );
    outputFile.write( reinterpret_cast<char *>(floatColumns.ptrData()), haloHeader.noColumnsFloats * columnNameLength * sizeof(char) );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the integer data
    buffer = dataIntegers.size()*sizeof(dataIntegers[0]);   //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&(dataIntegers[0])), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    // write the floating point data
    buffer = dataFloats.size()*sizeof(dataFloats[0]);   //total number of data bytes written
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    outputFile.write( reinterpret_cast<char *>(&(dataFloats[0])), buffer );
    outputFile.write( reinterpret_cast<char *>(&buffer), sizeof(buffer) );
    
    outputFile.close();
    std::cout << "Done.\n";
}



#endif
