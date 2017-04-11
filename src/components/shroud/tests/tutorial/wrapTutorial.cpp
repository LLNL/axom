// wrapTutorial.cpp
// This is generated code, do not edit
// wrapTutorial.cpp
#include "wrapTutorial.h"
#include <string>
#include "shroudrt.hpp"
#include "tutorial.hpp"

extern "C" {
namespace tutorial {

// void Function1()
// function_index=3
void TUT_function1()
{
// splicer begin function.function1
    Function1();
    return;
// splicer end function.function1
}

// double Function2(double arg1+intent(in)+value, int arg2+intent(in)+value)
// function_index=4
double TUT_function2(double arg1, int arg2)
{
// splicer begin function.function2
    double SH_rv = Function2(arg1, arg2);
    return SH_rv;
// splicer end function.function2
}

// void Sum(int len+intent(in)+value, int * values+dimension((*))+intent(in), int * result+intent(out))
// function_index=5
void TUT_sum(int len, int * values, int * result)
{
// splicer begin function.sum
    Sum(len, values, result);
    return;
// splicer end function.sum
}

// bool Function3(bool arg+intent(in)+value)
// function_index=6
bool TUT_function3(bool arg)
{
// splicer begin function.function3
    bool SH_rv = Function3(arg);
    return SH_rv;
// splicer end function.function3
}

// void Function3b(const bool arg1+intent(in)+value, bool * arg2+intent(out), bool * arg3+intent(inout))
// function_index=7
void TUT_function3b(const bool arg1, bool * arg2, bool * arg3)
{
// splicer begin function.function3b
    Function3b(arg1, arg2, arg3);
    return;
// splicer end function.function3b
}

// void Function4a(const std::string & arg1+intent(in)+len_trim(Larg1), const std::string & arg2+intent(in)+len_trim(Larg2), std::string * SH_F_rv+intent(out)+len(LSH_F_rv))
// function_index=34
void TUT_function4a_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2, char * SH_F_rv, int LSH_F_rv)
{
// splicer begin function.function4a_bufferify
    const std::string SH_arg1(arg1, Larg1);
    const std::string SH_arg2(arg2, Larg2);
    const std::string SH_rv = Function4a(SH_arg1, SH_arg2);
    shroud_FccCopy(SH_F_rv, LSH_F_rv, SH_rv.c_str());
    return;
// splicer end function.function4a_bufferify
}

// const std::string & Function4b(const std::string & arg1+intent(in), const std::string & arg2+intent(in))
// function_index=9
const char * TUT_function4b(const char * arg1, const char * arg2)
{
// splicer begin function.function4b
    const std::string SH_arg1(arg1);
    const std::string SH_arg2(arg2);
    const std::string & SH_rv = Function4b(SH_arg1, SH_arg2);
    return SH_rv.c_str();
// splicer end function.function4b
}

// void Function4b(const std::string & arg1+intent(in)+len_trim(Larg1), const std::string & arg2+intent(in)+len_trim(Larg2), std::string & output+intent(out)+len(Loutput))
// function_index=35
void TUT_function4b_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2, char * output, int Loutput)
{
// splicer begin function.function4b_bufferify
    const std::string SH_arg1(arg1, Larg1);
    const std::string SH_arg2(arg2, Larg2);
    const std::string & SH_rv = Function4b(SH_arg1, SH_arg2);
    shroud_FccCopy(output, Loutput, SH_rv.c_str());
    return;
// splicer end function.function4b_bufferify
}

// double Function5()
// function_index=24
double TUT_function5()
{
// splicer begin function.function5
    double SH_rv = Function5();
    return SH_rv;
// splicer end function.function5
}

// double Function5(double arg1+default(3.1415)+intent(in)+value)
// function_index=25
double TUT_function5_arg1(double arg1)
{
// splicer begin function.function5_arg1
    double SH_rv = Function5(arg1);
    return SH_rv;
// splicer end function.function5_arg1
}

// double Function5(double arg1+default(3.1415)+intent(in)+value, bool arg2+default(true)+intent(in)+value)
// function_index=10
double TUT_function5_arg1_arg2(double arg1, bool arg2)
{
// splicer begin function.function5_arg1_arg2
    double SH_rv = Function5(arg1, arg2);
    return SH_rv;
// splicer end function.function5_arg1_arg2
}

// void Function6(const std::string & name+intent(in))
// function_index=11
void TUT_function6_from_name(const char * name)
{
// splicer begin function.function6_from_name
    const std::string SH_name(name);
    Function6(SH_name);
    return;
// splicer end function.function6_from_name
}

// void Function6(const std::string & name+intent(in)+len_trim(Lname))
// function_index=37
void TUT_function6_from_name_bufferify(const char * name, int Lname)
{
// splicer begin function.function6_from_name_bufferify
    const std::string SH_name(name, Lname);
    Function6(SH_name);
    return;
// splicer end function.function6_from_name_bufferify
}

// void Function6(int indx+intent(in)+value)
// function_index=12
void TUT_function6_from_index(int indx)
{
// splicer begin function.function6_from_index
    Function6(indx);
    return;
// splicer end function.function6_from_index
}

// void Function7(int arg+intent(in)+value)
// function_index=26
void TUT_function7_int(int arg)
{
// splicer begin function.function7_int
    Function7<int>(arg);
    return;
// splicer end function.function7_int
}

// void Function7(double arg+intent(in)+value)
// function_index=27
void TUT_function7_double(double arg)
{
// splicer begin function.function7_double
    Function7<double>(arg);
    return;
// splicer end function.function7_double
}

// int Function8()
// function_index=28
int TUT_function8_int()
{
// splicer begin function.function8_int
    int SH_rv = Function8<int>();
    return SH_rv;
// splicer end function.function8_int
}

// double Function8()
// function_index=29
double TUT_function8_double()
{
// splicer begin function.function8_double
    double SH_rv = Function8<double>();
    return SH_rv;
// splicer end function.function8_double
}

// void Function9(double arg+intent(in)+value)
// function_index=15
void TUT_function9(double arg)
{
// splicer begin function.function9
    Function9(arg);
    return;
// splicer end function.function9
}

// void Function10()
// function_index=16
void TUT_function10_0()
{
// splicer begin function.function10_0
    Function10();
    return;
// splicer end function.function10_0
}

// void Function10(const std::string & name+intent(in), double arg2+intent(in)+value)
// function_index=17
void TUT_function10_1(const char * name, double arg2)
{
// splicer begin function.function10_1
    const std::string SH_name(name);
    Function10(SH_name, arg2);
    return;
// splicer end function.function10_1
}

// void Function10(const std::string & name+intent(in)+len_trim(Lname), double arg2+intent(in)+value)
// function_index=38
void TUT_function10_1_bufferify(const char * name, int Lname, double arg2)
{
// splicer begin function.function10_1_bufferify
    const std::string SH_name(name, Lname);
    Function10(SH_name, arg2);
    return;
// splicer end function.function10_1_bufferify
}

// int overload1(int num+intent(in)+value)
// function_index=30
int TUT_overload1_num(int num)
{
// splicer begin function.overload1_num
    int SH_rv = overload1(num);
    return SH_rv;
// splicer end function.overload1_num
}

// int overload1(int num+intent(in)+value, int offset+default(0)+intent(in)+value)
// function_index=31
int TUT_overload1_num_offset(int num, int offset)
{
// splicer begin function.overload1_num_offset
    int SH_rv = overload1(num, offset);
    return SH_rv;
// splicer end function.overload1_num_offset
}

// int overload1(int num+intent(in)+value, int offset+default(0)+intent(in)+value, int stride+default(1)+intent(in)+value)
// function_index=18
int TUT_overload1_num_offset_stride(int num, int offset, int stride)
{
// splicer begin function.overload1_num_offset_stride
    int SH_rv = overload1(num, offset, stride);
    return SH_rv;
// splicer end function.overload1_num_offset_stride
}

// int overload1(double type+intent(in)+value, int num+intent(in)+value)
// function_index=32
int TUT_overload1_3(double type, int num)
{
// splicer begin function.overload1_3
    int SH_rv = overload1(type, num);
    return SH_rv;
// splicer end function.overload1_3
}

// int overload1(double type+intent(in)+value, int num+intent(in)+value, int offset+default(0)+intent(in)+value)
// function_index=33
int TUT_overload1_4(double type, int num, int offset)
{
// splicer begin function.overload1_4
    int SH_rv = overload1(type, num, offset);
    return SH_rv;
// splicer end function.overload1_4
}

// int overload1(double type+intent(in)+value, int num+intent(in)+value, int offset+default(0)+intent(in)+value, int stride+default(1)+intent(in)+value)
// function_index=19
int TUT_overload1_5(double type, int num, int offset, int stride)
{
// splicer begin function.overload1_5
    int SH_rv = overload1(type, num, offset, stride);
    return SH_rv;
// splicer end function.overload1_5
}

// TypeID typefunc(TypeID arg+intent(in)+value)
// function_index=20
int TUT_typefunc(int arg)
{
// splicer begin function.typefunc
    TypeID SH_rv = typefunc(arg);
    return SH_rv;
// splicer end function.typefunc
}

// EnumTypeID enumfunc(EnumTypeID arg+intent(in)+value)
// function_index=21
int TUT_enumfunc(int arg)
{
// splicer begin function.enumfunc
    EnumTypeID SH_rv = enumfunc(static_cast<EnumTypeID>(arg));
    return static_cast<int>(SH_rv);
// splicer end function.enumfunc
}

// void useclass(const Class1 * arg1+intent(in)+value)
// function_index=22
void TUT_useclass(const TUT_class1 * arg1)
{
// splicer begin function.useclass
    useclass(static_cast<const Class1 *>(static_cast<const void *>(arg1)));
    return;
// splicer end function.useclass
}

// const std::string & LastFunctionCalled()+pure
// function_index=23
const char * TUT_last_function_called()
{
// splicer begin function.last_function_called
    const std::string & SH_rv = LastFunctionCalled();
    return SH_rv.c_str();
// splicer end function.last_function_called
}

// void LastFunctionCalled(std::string & SH_F_rv+intent(out)+len(LSH_F_rv))+pure
// function_index=39
void TUT_last_function_called_bufferify(char * SH_F_rv, int LSH_F_rv)
{
// splicer begin function.last_function_called_bufferify
    const std::string & SH_rv = LastFunctionCalled();
    shroud_FccCopy(SH_F_rv, LSH_F_rv, SH_rv.c_str());
    return;
// splicer end function.last_function_called_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace tutorial
}  // extern "C"
