//
// tutorial.cpp - A sample library to go along with the tutorial.
//

#include "tutorial.hpp"

namespace tutorial
{

void Function1(void)
{
    return;
}

// integer and real arguments
double Function2(double arg1, int arg2)
{
    return arg1 + arg2;
}

// bool arguments
bool Function3(bool arg)
{
    return ! arg;
}

// string
const std::string& Function4(const std::string& arg1, const std::string& arg2)
{
    static std::string rv(arg1 + arg2);
    return rv;
}


}  // namespace
