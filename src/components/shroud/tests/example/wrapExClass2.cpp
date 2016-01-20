// wrapExClass2.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapExClass2.cpp
#include "wrapExClass2.h"
#include "ExClass2.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace example {
namespace nested {

// ExClass2 * ExClass2(const string * name+intent(in))+constructor
// function_index=17
AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{

// splicer begin class.ExClass2.method.ex_class2
ExClass2 * rv = new ExClass2(name);
return static_cast<AA_exclass2 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.ex_class2
}

// ExClass2 * ExClass2(const string * name+intent(in)+len_trim(Lname))+constructor
// function_index=34
AA_exclass2 * AA_exclass2_ex_class2_bufferify(const char * name, int Lname)
{

// splicer begin class.ExClass2.method.ex_class2_bufferify
ExClass2 * rv = new ExClass2(std::string(name, Lname));
return static_cast<AA_exclass2 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.ex_class2_bufferify
}

// void delete()+destructor
// function_index=18
void AA_exclass2_delete(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.delete
delete selfobj;
// splicer end class.ExClass2.method.delete
}

// const string & getName() const
// function_index=19
const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.ExClass2.method.get_name
}

// const int GetNameLength()
// function_index=20
/**
 * \brief helper function for Fortran
 *
 */
const int AA_exclass2_get_name_length(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name_length
return selfobj->getName().length();
// splicer end class.ExClass2.method.get_name_length
}

// ExClass1 * get_class1(const ExClass1 * in+intent(in)+value)
// function_index=21
AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, const AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_class1
ExClass1 * rv = selfobj->get_class1(static_cast<ExClass1 *>(static_cast<void *>(in)));
return static_cast<AA_exclass1 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.get_class1
}

// void * declare(TypeID type+intent(in)+value)
// function_index=27
void AA_exclass2_declare_0(AA_exclass2 * self, int type)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.declare_0
selfobj->declare(getTypeID(type));
return;
// splicer end class.ExClass2.method.declare_0
}

// void * declare(TypeID type+intent(in)+value, SidreLength len+default(1)+intent(in)+value)
// function_index=22
void AA_exclass2_declare_1(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.declare_1
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.ExClass2.method.declare_1
}

// void destroyall()
// function_index=23
void AA_exclass2_destroyall(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.destroyall
selfobj->destroyall();
return;
// splicer end class.ExClass2.method.destroyall
}

// TypeID getTypeID() const
// function_index=24
int AA_exclass2_get_type_id(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.ExClass2.method.get_type_id
}

// void setValue(int value+intent(in)+value)
// function_index=28
void AA_exclass2_set_value_int(AA_exclass2 * self, int value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_int
selfobj->setValue<int>(value);
return;
// splicer end class.ExClass2.method.set_value_int
}

// void setValue(long value+intent(in)+value)
// function_index=29
void AA_exclass2_set_value_long(AA_exclass2 * self, long value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_long
selfobj->setValue<long>(value);
return;
// splicer end class.ExClass2.method.set_value_long
}

// void setValue(float value+intent(in)+value)
// function_index=30
void AA_exclass2_set_value_float(AA_exclass2 * self, float value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_float
selfobj->setValue<float>(value);
return;
// splicer end class.ExClass2.method.set_value_float
}

// void setValue(double value+intent(in)+value)
// function_index=31
void AA_exclass2_set_value_double(AA_exclass2 * self, double value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_double
selfobj->setValue<double>(value);
return;
// splicer end class.ExClass2.method.set_value_double
}

// int getValue()
// function_index=32
int AA_exclass2_get_value_int(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_value_int
int rv = selfobj->getValue<int>();
return rv;
// splicer end class.ExClass2.method.get_value_int
}

// double getValue()
// function_index=33
double AA_exclass2_get_value_double(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_value_double
double rv = selfobj->getValue<double>();
return rv;
// splicer end class.ExClass2.method.get_value_double
}

// splicer begin class.ExClass2.additional_functions
// splicer end class.ExClass2.additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"
