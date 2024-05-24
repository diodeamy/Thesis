#ifndef MEMORYALLOCATIONCHECK_HEADER
#define MEMORYALLOCATIONCHECK_HEADER

#include <iostream>


/* Check to see if memory allocation was succesful (if the pointer points to a memory location) or unsuccesful (if the pointer is null).
The variable message holds a message to be sent to the user in case the pointer is null.
*/
template <typename T1, typename STRING>
void memoryAllocationCheck(T1 *pointer,
                           STRING const message)
{
    if ( pointer==0 )
    {
        std::cout << "\n\n~~~ ERROR ~~~ The program could not initialize the memory needed to store the data. '" << message << "' Try again with a smaller number of particles or on a machine with a bigger RAM memory!\n";
        exit( EXIT_FAILURE );
    }
}



#endif
