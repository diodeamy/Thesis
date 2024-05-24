#include <iostream>
#include <string>
#include <vector>
#include <cstring>


#include "miscellaneous.h"
using namespace std;



/* This functions reads an input char array and extracts a vector of integers from it.
The integers can be separated by 2 different chars:
    interval - when present between two integers shows that all integers in the range between the two were supplied (2nd number must be larger than first)
    split - separate two integers that aren't the limits of an interval
*/
void getIntegerVector(char input[],
                      size_t const length,
                      vector<int> *scales,
                      char interval,
                      char split)
{
    bool checkInterval = false;   //variable to keep track if there is an interval in the input char array
    bool nonNumber = true;        // true if an 'interval' or a 'split' was just read, false when reading a number
    int temp = 0;                 //intermediate variable to store the scale or part of it
    
    for (size_t i=0; i<length; ++i)
    {
        if ( input[i]>='0' and input[i]<='9' )
        {
            temp = temp*10 + input[i] - '0';
            nonNumber = false;
        }
        else if ( nonNumber )
        {
             cout << "~~~ ERROR ~~~ Each '" << interval << "' and/or '" << split << "' has to be preceded and followed by an integer number. Your entry did not satisfy this condition. The program ended unsuccesfully!\n";
             exit( EXIT_FAILURE );
        }
        else if ( input[i]==interval )
        {
            nonNumber = true;
            if ( checkInterval )     //if interval insert, insert interval
            {
                 int temp2 = scales->back();
                 if ( temp2 > temp )
                 {
                      cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
                      exit( EXIT_FAILURE );
                 }
                 for (int j=temp2+1; j<=temp; ++j)
                     scales->push_back( j );
                 temp = 0;
            }
            else //insert the number and rememeber that an interval insert started
            {
                checkInterval = true;
                scales->push_back( temp );
                temp = 0;
            }
        }
        else if ( input[i]==split )
        {
            nonNumber = true;
            if ( checkInterval )     //if interval insert, insert interval
            {
                 int temp2 = scales->back();
                 if ( temp2 > temp )
                 {
                      cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
                      exit( EXIT_FAILURE );
                 }
                 for (int j=temp2+1; j<=temp; ++j)
                     scales->push_back( j );
                 temp = 0;
                 checkInterval = false;
            }
            else //insert the number and rememeber that an interval insert started
            {
                scales->push_back( temp );
                temp = 0;
            }
        }
        else
        {
            cout << "~~~ ERROR ~~~ When inserting an integer interval, the only delimiters you are allowed to use are: '" << interval << "' and '" << split << "'. The program ended unsuccesfully!\n";
            exit( EXIT_FAILURE );
        }
    }
    
    
    // still need to save the last integer/sequence of integer number/s
    if ( checkInterval )     //if interval insert, insert interval
    {
        int temp2 = scales->back();
        if ( temp2 > temp )
        {
            cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
            exit( EXIT_FAILURE );
        }
        for (int j=temp2+1; j<=temp; ++j)
            scales->push_back( j );
        temp = 0;
        checkInterval = false;
    }
    else //insert the number
    {
        scales->push_back( temp );
        temp = 0;
    }
}
void getIntegerVector(string input,
                      vector<int> *scales,
                      char interval,
                      char split)
{
    bool checkInterval = false;   //variable to keep track if there is an interval in the input char array
    bool nonNumber = true;        // true if an 'interval' or a 'split' was just read, false when reading a number
    int temp = 0;                 //intermediate variable to store the scale or part of it
    
    for (size_t i=0; i<input.length(); ++i)
    {
        if ( input[i]>='0' and input[i]<='9' )
        {
            temp = temp*10 + input[i] - '0';
            nonNumber = false;
        }
        else if ( nonNumber )
        {
             cout << "~~~ ERROR ~~~ Each '" << interval << "' and/or '" << split << "' has to be preceded and followed by an integer number. Your entry did not satisfy this condition. The program ended unsuccesfully!\n";
             exit( EXIT_FAILURE );
        }
        else if ( input[i]==interval )
        {
            nonNumber = true;
            if ( checkInterval )     //if interval insert, insert interval
            {
                 int temp2 = scales->back();
                 if ( temp2 > temp )
                 {
                      cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
                      exit( EXIT_FAILURE );
                 }
                 for (int j=temp2+1; j<=temp; ++j)
                     scales->push_back( j );
                 temp = 0;
            }
            else //insert the number and rememeber that an interval insert started
            {
                checkInterval = true;
                scales->push_back( temp );
                temp = 0;
            }
        }
        else if ( input[i]==split )
        {
            nonNumber = true;
            if ( checkInterval )     //if interval insert, insert interval
            {
                 int temp2 = scales->back();
                 if ( temp2 > temp )
                 {
                      cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
                      exit( EXIT_FAILURE );
                 }
                 for (int j=temp2+1; j<=temp; ++j)
                     scales->push_back( j );
                 temp = 0;
                 checkInterval = false;
            }
            else //insert the number and remember that an interval insert started
            {
                scales->push_back( temp );
                temp = 0;
            }
        }
        else
        {
            cout << "~~~ ERROR ~~~ When inserting an integer interval, the only delimiters you are allowed to use are: '" << interval << "' and '" << split << "'. The program ended unsuccesfully!\n";
            exit( EXIT_FAILURE );
        }
    }
    
    
    // still need to save the last integer/sequence of integer number/s
    if ( checkInterval )     //if interval insert, insert interval
    {
        int temp2 = scales->back();
        if ( temp2 > temp )
        {
            cout << "~~~ ERROR ~~~ When inserting an integer interval, the left interval limit must be <= than right interval limit. You inserted interval '" << temp2 << interval << temp << "' which does not satisfy the previous condition. The program ended unsuccesfully!\n";
            exit( EXIT_FAILURE );
        }
        for (int j=temp2+1; j<=temp; ++j)
            scales->push_back( j );
        temp = 0;
        checkInterval = false;
    }
    else //insert the number
    {
        scales->push_back( temp );
        temp = 0;
    }
}



/* Copy from a char array to another. It check if the array to be coppied to is smaller than the initial one, than it copies only the first characters that fit into the array. */
void copyStringFront(char const *oldString,
                     size_t const oldSize,
                     char *newString,
                     size_t const newSize,
                     string message)
{
    if ( oldSize<newSize )
       strcpy( newString, oldString );
    else
    {
        throwWarning( message );
        strncpy( newString, oldString, newSize-1 );
        strcpy( &(newString[newSize-1]), "" );
    }
}


/* Copy from a char array to another. It check if the array to be coppied to is smaller than the initial one, than it copies only the last characters that fit into the array.
*/
void copyStringEnd(char const *oldString,
                   size_t const oldSize,
                   char *newString,
                   size_t const newSize,
                   string message)
{
    if ( oldSize<newSize )
       strcpy( newString, oldString );
    else
    {
        throwWarning( message );
        int temp = oldSize - newSize;
        strcpy( newString, &(oldString[temp]) );
    }
}

/* Copy from a char array to another. It check if the array to be coppied to is smaller than the initial one, than it copies only the first characters that fit into the array.
*/
void copyStringFront(string const old,
                     char *newString,
                     size_t const newSize,
                     string message)
{
    size_t const oldSize = old.length();
    if ( oldSize<newSize )
       strcpy( newString, old.c_str() );
    else
    {
        throwWarning( message );
        strncpy( newString, old.c_str(), newSize-1 );
        strcpy( &(newString[newSize-1]), "" );
    }
}


/* Copy from a char array to another. It check if the array to be coppied to is smaller than the initial one, than it copies only the last characters that fit into the array.
*/
void copyStringEnd(string const old,
                   char *newString,
                   size_t const newSize,
                   string message)
{
    size_t const oldSize = old.length();
    if ( oldSize<newSize )
       strcpy( newString, old.c_str() );
    else
    {
        throwWarning( message );
        int temp = oldSize - newSize -1;
        strcpy( newString, &((old.c_str())[temp]) );
        strcpy( &(newString[newSize-1]), "" );
    }
}


bool duplicateEntries(int *array,
                      int const size)
{
    for (int i=0; i<size; ++i)
        for (int j=i+1; j<size; ++j)
            if ( array[i]==array[j] )
                return true;
    return false;
}


/* Gets the full file path. */
string fullFilePath(string filename)
{
    string currentPath = getCurrentPath();
    char *temp = (char*) currentPath.c_str();
    size_t tempLength = currentPath.length();
    
    // loop over the input string and find how many instances of '../' there are
    size_t i=0, leftDirectories = 0, j=0;
    bool copy = false;
    char tempFilename[filename.length()+1];
    while (i<filename.length() )
    {
        if ( filename[i-1]=='/' and filename[i]=='/' ) // if finding '//' keep only one
            ++i;
        else if ( filename[i]=='.' and filename[i+1]=='/' and not copy )  // skip the './' characters which show the curent directory
            i += 2;
        else if ( filename[i]=='.' and filename[i+1]=='.' and filename[i+2]=='/' and not copy ) // skip the '../' characters - but keep track to remove the left directory from curent path
        {
            i += 3;
            ++leftDirectories;
        }
        else
        {
            copy = true;
            tempFilename[j++] = filename[i++];
        }
    }
    tempFilename[j] = '\0';
    
    // delete directories from curent path if '../' was found
    while (leftDirectories>0)
    {
        size_t i1 = tempLength-1;
        while (temp[i1]!='/' and i1>0)
            temp[i1--] = '\0';
        temp[i1] = '\0';  // delete also the '/' character
        tempLength = i1;
        --leftDirectories;
    }
    
    // get the result
    string result = temp;
    result = result + '/' + tempFilename;
    return result;
}

/* Initiates a char array to '' default values for all the elements. */
void initializeCharArray(char *array, size_t const size)
{
    for (size_t i=0; i<size; ++i)
        array[i] = '\0';
}

