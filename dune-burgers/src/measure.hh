#ifndef MEASURE_HH
#define MEASURE_HH


// The following is taken from
// http://stackoverflow.com/questions/2808398/easily-measure-elapsed-time

#include <iostream>
#include <chrono>

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
    template<typename F, typename ...Args>
    static typename TimeT::rep execution(F func, Args&&... args)
    {
        auto start = std::chrono::system_clock::now();
        func(std::forward<Args>(args)...);
        auto duration = std::chrono::duration_cast< TimeT>
                            (std::chrono::system_clock::now() - start);
        return duration.count();
    }
};


#endif // MEASURE_HH
