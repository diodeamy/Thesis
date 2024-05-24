
#ifndef OPTIONSCHECK_HEADER
#define OPTIONSCHECK_HEADER

#include "throwError.h"

/*
This header can be use only in conjunction with the 'boost/program_options.hpp' library.
*/



/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
template <typename T>
void conflicting_options(const T &vm,
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted())
    {
        std::cout << "Conflicting options '" << opt1 << "' and '" << opt2 << "'!\n";
        exit( EXIT_FAILURE );
    }
}


/* This function checks that at most only one of the options given in the option array was inserted by the user. */
template <typename T>
void conflicting_options(const T &vm,
                         int const noOptions, const char** options)
{
    int optionsOn = 0;  // keep track of the number of options from the option array supplied by the user
    for (int i=0; i<noOptions; ++i)
        if ( vm.count(options[i]) && !vm[options[i]].defaulted() )
            ++optionsOn;
    
    if ( optionsOn>1 )
    {
        std::cout << "You can select at most one option from the following option list: ";
        for (int i=0; i<noOptions; ++i)
            std::cout << " '" << options[i] << "'";
        std::cout << ". No two option from the list can be used at the same time!\n";
        exit( EXIT_FAILURE );
    }
}


/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
template <typename T>
void option_dependency(const T &vm,
                       const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
        {
            std::cout << "Option '" << for_what << "' requires option '" << required_option << "'!\n";
            exit( EXIT_FAILURE );
        }
}


/* Function used to check that 'opt1' is not supplied in the absence of 'opt2'. */
template <typename T>
void superfluous_options(const T &vm,
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted())
        if (vm.count(opt2) == 0 || vm[opt2].defaulted())
        {
            std::cout << "The option '" << opt1 << "' can be used only in the presence of '" << opt2 << "'. It does not make sense to use '" << opt1 << "' otherwise!\n";
            exit( EXIT_FAILURE );
        }
}
/* Function used to check that 'opt1' is not supplied if 'optionsOn' is false.
NOTE: 'OptionOn' should be true if one or several options were suplied. Give the name of the options in 'optionsName'. */
template <typename T>
void superfluous_options(const T &vm,
                         const char* opt1,
                         const bool optionsOn,
                         std::string const optionsName)
{
    if (vm.count(opt1) && !vm[opt1].defaulted())
        if ( not optionsOn )
        {
            std::cout << "The option '" << opt1 << "' can be used only in the presence of option/s: " << optionsName << ". It does not make sense to use '" << opt1 << "' otherwise!\n";
            exit( EXIT_FAILURE );
        }
}


#endif
