//
// tutorial.hpp - wrapped routines
//

#include "tutorial.hpp"

// These variables exist to avoid warning errors
static std::string rv_str;
static int rv_int;
static double rv_double;


void Function1()
{
    return;
}

double Function2(double arg1, int arg2)
{
    return arg1 + arg2;
}

bool Function3(bool arg)
{
    return ! arg;
}

const std::string& Function4a(const std::string& arg1, const std::string& arg2)
{
    rv_str = arg1 + arg2;
    return rv_str;
}
const std::string& Function4b(const std::string& arg1, const std::string& arg2)
{
    rv_str = arg1 + arg2;
    return rv_str;
}

double Function5(double arg1, int arg2)
{
    return arg1 + arg2;
}

void Function6(const std::string& name)
{
    rv_str = name;  // avoid unused-parameter
    return;
}
void Function6(int indx)
{
    rv_int = indx;
    return;
}

template<>
void Function7<int>(int arg)
{
    rv_int = arg;
}

template<>
void Function7<double>(double arg)
{
    rv_double = arg;
}

template<>
int Function8<int>()
{
    return rv_int;
}

template<>
double Function8<double>()
{
    return rv_double;
}

void Function9(double arg)
{
    rv_double = arg;
    return;
}

void Class1::Method1()
{
    return;
}
