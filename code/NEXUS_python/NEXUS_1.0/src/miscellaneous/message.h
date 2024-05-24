#ifndef MESSAGE_HEADER
#define MESSAGE_HEADER

#include <iostream>
#include <iomanip>

// function to output progress message
struct ProgressMessage
{
    int _percentangeDone;
    // Class constructor
    ProgressMessage() {  _percentangeDone = -1; }
    
    // Overaloads the << which is used to send messages to the user
    template <typename T>
    inline ProgressMessage& operator << (T right)
    {
        std::cout << right;
        return *this;
    }
    
    // Show progress to the user
    inline void updateProgress(int const percentageDone)
    {
        if (_percentangeDone==percentageDone )
            return;
        _percentangeDone = percentageDone;
        char const del[4] = "\b\b\b";
        std::cout << std::setw(2) << _percentangeDone << "\%" << std::flush << del;
    }
};


#endif



