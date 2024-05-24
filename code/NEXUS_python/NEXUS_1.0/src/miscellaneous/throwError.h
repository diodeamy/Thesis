#ifndef THROWERROR_HEADER
#define THROWERROR_HEADER

#include <iostream>
#include <string>
#include <cstdlib>

// function to output progress message
struct Error
{
    int iterationNo;
    Error(void) { iterationNo = 0; }
    template <typename T>
    inline Error& operator << (T right)
    {
        if ( (iterationNo++)==0 )
            std::cout << "\n\n~~~ ERROR ~~~ ";
        std::cout << right;
        return *this;
    }
    void EndErrorMessage()
    {
        std::cout << " The program ended unsuccessfully!\n\n";
        exit( EXIT_FAILURE );
    }
};



/* Throw error and stop the program. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3,
                T4 message_4,
                T5 message_5,
                T6 message_6,
                T6 message_7)
{
    Error error;
    error << message_1 << message_2 << message_3 << message_4 << message_5 << message_6 << message_7;
    error.EndErrorMessage();
}
template <typename T1, typename T2, typename T3, typename T4>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3,
                T4 message_4)
{
    Error error;
    error << message_1 << message_2 << message_3 << message_4;
    error.EndErrorMessage();
}
template <typename T1, typename T2, typename T3>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3)
{
    throwError<T1,T2,T3,std::string>( message_1, message_2, message_3, "" );
}
template <typename T1, typename T2, typename T3>
void throwError(T1 message_1,
                T2 message_2)
{
    throwError<T1,T2,std::string,std::string>( message_1, message_2, "", "" );
}
template <typename T1>
void throwError(T1 message_1)
{
    throwError<T1,std::string,std::string,std::string>( message_1, "", "", "" );
}
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3,
                T4 message_4,
                T5 message_5)
{
    throwError<T1,T2,T3,T4,T5,std::string,std::string>( message_1, message_2, message_3, message_4, message_5, "", "" );
}
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void throwError(T1 message_1,
                T2 message_2,
                T3 message_3,
                T4 message_4,
                T5 message_5,
                T6 message_6)
{
    throwError<T1,T2,T3,T4,T5,T6,std::string>( message_1, message_2, message_3, message_4, message_5, message_6, "" );
}


#endif
