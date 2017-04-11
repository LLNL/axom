//
// tutorial.hpp - wrapped routines
//

#ifndef TUTORIAL_HPP
#define TUTORIAL_HPP

#include <string>

namespace tutorial
{

enum EnumTypeID {
    ENUM0,
    ENUM1,
    ENUM2
};

typedef int TypeID;

void Function1();

double Function2(double arg1, int arg2);

bool Function3(bool arg);
void Function3b(const bool arg1, bool *arg2, bool *arg3);

const std::string Function4a(const std::string& arg1, const std::string& arg2);
const std::string& Function4b(const std::string& arg1, const std::string& arg2);

double Function5(double arg1 = 3.1415, bool arg2 = true);

void Function6(const std::string& name);
void Function6(int indx);

// specialize for int and double in tutorial.cpp
template<typename ArgType>
void Function7(ArgType arg);

// specialize for int and double in tutorial.cpp
template<typename RetType>
RetType Function8();

void Function9(double arg);

void Function10();
void Function10(const std::string &name, double arg2);

void Sum(int len, int * values, int *result);

int overload1(int num, int offset = 0, int stride = 1);
int overload1(double type, int num, int offset = 0, int stride = 1);

TypeID typefunc(TypeID arg);

EnumTypeID enumfunc(EnumTypeID arg);

const std::string& LastFunctionCalled();

class Class1
{
public:
    void Method1();
};

void useclass(const Class1 *arg);
void getclass(const Class1 **arg);

} /* end namespace tutorial */

#endif // TUTORIAL_HPP
