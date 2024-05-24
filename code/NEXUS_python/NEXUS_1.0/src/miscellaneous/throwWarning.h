#ifndef THROWWARNING_HEADER
#define THROWWARNING_HEADER

#include <iostream>
#include <string>


/* Throw warning, but do not stop the program. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
void throwWarning(T1 message_1,
                    T2 message_2,
                    T3 message_3,
                    T4 message_4,
                    T5 message_5,
                    T6 message_6,
                    T6 message_7)
{
    std::cout << "\n~~~ WARNING ~~~ " << message_1 << message_2 << message_3 << message_4 << message_5 << message_6 << message_7 << "\n\n";
}
template <typename T1, typename T2, typename T3, typename T4>
void throwWarning(T1 message_1,
                  T2 message_2,
                  T3 message_3,
                  T4 message_4)
{
    std::cout << "\n~~~ WARNING ~~~ " << message_1 << message_2 << message_3 << message_4 << "\n\n";
}
template <typename T1, typename T2, typename T3>
void throwWarning(T1 message_1,
                    T2 message_2,
                    T3 message_3)
{
    throwWarning<T1,T2,T3,std::string>( message_1, message_2, message_3, "" );
}
template <typename T1, typename T2>
void throwWarning(T1 message_1,
                    T2 message_2)
{
    throwWarning<T1,T2,std::string,std::string>( message_1, message_2, "", "" );
}
template <typename T1>
void throwWarning(T1 message_1)
{
    throwWarning<T1,std::string,std::string,std::string>( message_1, "", "", "" );
}
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void throwWarning(T1 message_1,
                    T2 message_2,
                    T3 message_3,
                    T4 message_4,
                    T5 message_5)
{
    throwWarning<T1,T2,T3,T4,T5,std::string,std::string>( message_1, message_2, message_3, message_4, message_5, "", "" );
}
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void throwWarning(T1 message_1,
                    T2 message_2,
                    T3 message_3,
                    T4 message_4,
                    T5 message_5,
                    T6 message_6)
{
    throwWarning<T1,T2,T3,T4,T5,T6,std::string>( message_1, message_2, message_3, message_4, message_5, message_6, "" );
}


#endif
