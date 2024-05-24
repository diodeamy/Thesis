#ifndef RAW_HEADER
#define RAW_HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <miscellaneous.h>


/* Save the quantity as a raw binary file. */
template <typename T>
void outputRaw(T *delta,
               size_t const totalSize,
               std::string outputFilename)
{
    std::cout << "Writing the raw binary file '" << outputFilename << "' ... " << std::flush;
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, outputFilename );
    outputFile.write( reinterpret_cast<char *>(delta), totalSize*sizeof(delta[0]) );
    outputFile.close();
    std::cout << "Done.\n";
}


/* Save the results to a binary file with a descriptive text file. */
template <typename T>
void outputBinaryFileWithInfo(std::string outputFilename,
                              std::string description,
                              T *data,
                              size_t const totalSize)
{
    std::string infoFile = outputFilename + ".info";
    std::cout << "Writing the binary file '" << outputFilename << "' with an additional information file '" << infoFile << "' ... " << std::flush;
    
    // first write the information file
    std::fstream outputFile;
    openOutputTextFile( outputFile, infoFile );
    outputFile << description;
    outputFile.close();
    
    // now write the data file
    openOutputBinaryFile( outputFile, outputFilename );
    outputFile.write( reinterpret_cast<char *>(data), totalSize*sizeof(data[0]) );
    outputFile.close();
    std::cout << "Done.\n";
}



#endif
