//
// tutorial.hpp - wrapped routines
//

#include "tutorial.hpp"

namespace tutorial
{

static std::string last_function_called;

// These variables exist to avoid warning errors
static std::string global_str;
static int global_int;
static double global_double;




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
    global_str = arg1 + arg2;
    return global_str;
}
const std::string Function4b(const std::string& arg1, const std::string& arg2)
{
    last_function_called = "Function4b";
    return arg1 + arg2;
}

double Function5(double arg1, bool arg2)
{
    last_function_called = "Function5";
    if (arg2) {
	return arg1 + 10.0;
    } else {
	return arg1;
    }
}

void Function6(const std::string& name)
{
    last_function_called = "Function6(string)";
    global_str = name;
    return;
}
void Function6(int indx)
{
    last_function_called = "Function6(int)";
    global_int = indx;
    return;
}

template<>
void Function7<int>(int arg)
{
    last_function_called = "Function7<int>";
    global_int = arg;
}

template<>
void Function7<double>(double arg)
{
    last_function_called = "Function7<double>";
    global_double = arg;
}

template<>
int Function8<int>()
{
    last_function_called = "Function8<int>";
    return global_int;
}

template<>
double Function8<double>()
{
    last_function_called = "Function8<double>";
    return global_double;
}

void Function9(double arg)
{
    last_function_called = "Function9";
    global_double = arg;
    return;
}

void Function10()
{
    last_function_called = "Function10_0";
}

void Function10(const std::string &name, double arg2)
{
    last_function_called = "Function10_1";
    global_str = name;
    global_double = arg2;
}

void Sum(int len, int *values, int *result)
{
    last_function_called = "Sum";

    int sum = 0;
    for (int i=0; i < len; i++) {
	sum += values[i];
    }
    *result = sum;
    return;
}

int overload1(int num, int offset, int stride)
{
    last_function_called = "overload1_0";
    return num + offset * stride;
    
}

int overload1(double type, int num, int offset, int stride)
{
    last_function_called = "overload1_1";
    global_double = type;
    return num + offset * stride;
}

TypeID typefunc(TypeID arg)
{
    last_function_called = "typefunc";
    return static_cast<int>(arg);
}

EnumTypeID enumfunc(EnumTypeID arg)
{
    last_function_called = "enumfunc";
    switch (arg) {
    default:
	return ENUM2;
    }
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

} /* end namespace tutorial */
