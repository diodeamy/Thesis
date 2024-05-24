
#ifndef TIMER_CLASS_HEADER
#define TIMER_CLASS_HEADER

#include <boost/timer.hpp>
#include <string>
#include <ctime>


struct Timer
{
    boost::timer t;
    std::time_t userT;
    
    // functions
    void start()        // start the timer
    {
        t.restart();
        std::time( &userT );
    }
    void printTime(std::string variableName = "")       // print the time passed
    {
        std::time_t end;
        std::time( &end );
        std::cout << "\t <<< '" << variableName << "': user time = " << std::difftime(end,userT) << " sec. and CPU time = " << t.elapsed() << " sec.\n";
    }
    double cpuTime()    // return the CPU time
    {
        return t.elapsed();
    }
    double userTime()   // return the time passed by the user
    {
        std::time_t end;
        std::time( &end );
        return std::difftime( end,userT );
    }
    void limits()
    {
        std::cout << "timer::elapsed_min() reports " << t.elapsed_min() << " seconds\n";
        std::cout << "timer::elapsed_max() reports " << t.elapsed_max()
                << " seconds, which is " << t.elapsed_max()/3600.0 << " hours\n";
    }
};


#endif
