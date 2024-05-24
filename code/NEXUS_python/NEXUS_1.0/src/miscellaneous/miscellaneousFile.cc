#include <iostream>
#include <string>
#include <fstream>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;


#include "miscellaneous.h"
using namespace std;




/* This function opens a binary file for input (to read from the file). It checks that the operation was done succesfully.
*/
void openInputBinaryFile(fstream & inputFile,
                         string & fileName)
{
    char const *temp = fileName.c_str();
    inputFile.open( temp, ios::in | ios::binary );
    
    if ( not inputFile.is_open() )
    {
        cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for reading! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}




/* This function opens a binary file for output (to write data to the file). It checks that the operation was done succesfully.
*/
void openOutputBinaryFile(fstream & outputFile,
                          string & fileName)
{
    char const *temp = fileName.c_str();
    outputFile.open( temp, ios::out | ios::binary );
    
    if ( not outputFile.is_open() )
    {
        cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for writing! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}


/* This function opens a text file for input (to read from the file). It checks that the operation was done succesfully.
*/
void openInputTextFile(fstream & inputFile,
                       string & fileName)
{
    char const *temp = fileName.c_str();
    inputFile.open( temp, ios::in );
    
    if ( not inputFile.is_open() )
    {
        cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for reading! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}




/* This function opens a text file for output (to write data to the file). It checks that the operation was done succesfully.
*/
void openOutputTextFile(fstream & outputFile,
                        string & fileName)
{
    char const *temp = fileName.c_str();
    outputFile.open( temp, ios::out );
    
    if ( not outputFile.is_open() )
    {
        cout << "~~~ ERROR ~~~ The file '" << fileName << "' could not be opened for writing! The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}



/* This function checks to see if a file exists.
*/
void existentFile(string const & fileName)
{
    if ( not bfs::exists(fileName) )
    {
        cout << "~~~ ERROR ~~~ The file '" << fileName << "' does not exist! Choose a file that exists. The program ended unsucessfully!\n";
        exit( EXIT_FAILURE );
    }
}



/* This function finds the full path of a file/directory and writes it to a char array. If the name is longer than 'maxSize' than it writes only the last 'maxSize' characters.
*/
void getFullPath(string &fileName,
                 char output[],
                 int const maxSize)
{
    bfs::path full_path( bfs::initial_path<bfs::path>() );
    full_path = bfs::system_complete( bfs::path( fileName ) );
    string fullOutput = full_path.string();
    
    // now check that the path name is not longer than maxSize
    int offset = 0;
    if ( fullOutput.length()>size_t(maxSize-1) ) offset = fullOutput.length() - maxSize +1;  //if path too long, copy only the last 'maxSize' characters of the path
    
    size_t i;
    for ( i=0; i<fullOutput.length()-offset; ++i)
        output[i] = fullOutput[i+offset];
    output[i] = char(0);
}


/* Returns the current directory. */
string getCurrentPath()
{
    bfs::path full_path( bfs::initial_path<bfs::path>() );
    return full_path.string();
}

/* Returns the current filename. */
string getFilename(string fileName)
{
    return bfs::path( fileName ).filename().string();
}








