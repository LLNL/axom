//
// tutorial.hpp - wrapped routines
//

#include "tutorial.hpp"

static std::string last_function_called;

// These variables exist to avoid warning errors
static std::string rv_str;
static int rv_int;
static double rv_double;




void Function1()
{
    last_function_called = "Function1";
    return;
}

double Function2(double arg1, int arg2)
{
    last_function_called = "Function2";
    return arg1 + arg2;
}

bool Function3(bool arg)
{
    last_function_called = "Function3";
    return ! arg;
}

const std::string& Function4a(const std::string& arg1, const std::string& arg2)
{
    last_function_called = "Function4a";
    rv_str = arg1 + arg2;
    return rv_str;
}
const std::string& Function4b(const std::string& arg1, const std::string& arg2)
{
    last_function_called = "Function4b";
    rv_str = arg1 + arg2;
    return rv_str;
}

double Function5(double arg1, int arg2)
{
    last_function_called = "Function5";
    return arg1 + arg2;
}

void Function6(const std::string& name)
{
    last_function_called = "Function6(string)";
    rv_str = name;
    return;
}
void Function6(int indx)
{
    last_function_called = "Function6(int)";
    rv_int = indx;
    return;
}

template<>
void Function7<int>(int arg)
{
    last_function_called = "Function7<int>";
    rv_int = arg;
}

template<>
void Function7<double>(double arg)
{
    last_function_called = "Function7<double>";
    rv_double = arg;
}

template<>
int Function8<int>()
{
    last_function_called = "Function8<int>";
    return rv_int;
}

template<>
double Function8<double>()
{
    last_function_called = "Function8<double>";
    return rv_double;
}

void Function9(double arg)
{
    last_function_called = "Function9";
    rv_double = arg;
    return;
}

//----------------------------------------------------------------------

void Class1::Method1()
{
    last_function_called = "Class1::Method1";
    return;
}

//----------------------------------------------------------------------

const std::string& LastFunctionCalled()
{
    return last_function_called;
}

