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
const std::string& Function4a(const std::string& arg1, const std::string& arg2)
{
    static std::string rv(arg1 + arg2);
    return rv;
}

const std::string& Function4b(const std::string& arg1, const std::string& arg2)
{
    static std::string rv(arg1 + arg2);
    return rv;
}

double Function5(double arg1, int arg2)
{
    return arg1 + arg2;
}

void Function6(const std::string &name)
{
    return;
}

void Function6(int indx)
{
    return;
}

void Function9(double arg)
{
}

int Sum(int len, int *values)
{
    int sum;
    for (int i=0; i < len; i++) {
	sum += values[i];
    }
    return sum;
}

std::string global_string;
const std::string& LastFunctionCalled()
{
    return global_string;
}



}  // namespace
