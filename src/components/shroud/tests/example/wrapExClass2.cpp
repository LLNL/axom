// wrapExClass2.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapExClass2.cpp
#include "wrapExClass2.h"
#include <string>
#include "ExClass2.hpp"
#include "shroudrt.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace example {
namespace nested {

// ExClass2 * ExClass2(const string * name+intent(in))+constructor
// function_index=18
AA_exclass2 * AA_exclass2_ex_class2(const char * name)
{

// splicer begin class.ExClass2.method.ex_class2
const std::string SH_name(name);
ExClass2 * rv = new ExClass2(SH_name);
return static_cast<AA_exclass2 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.ex_class2
}

// ExClass2 * ExClass2(const string * name+intent(in)+len_trim(Lname))+constructor
// function_index=38
AA_exclass2 * AA_exclass2_ex_class2_bufferify(const char * name, int Lname)
{

// splicer begin class.ExClass2.method.ex_class2_bufferify
const std::string SH_name(name, Lname);
ExClass2 * rv = new ExClass2(SH_name);
return static_cast<AA_exclass2 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.ex_class2_bufferify
}

// void delete()+destructor
// function_index=19
void AA_exclass2_delete(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.delete
delete selfobj;
// splicer end class.ExClass2.method.delete
}

// const string & getName() const
// function_index=20
const char * AA_exclass2_get_name(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.ExClass2.method.get_name
}

// void getName(string & SH_F_rv+intent(out)+len(LSH_F_rv)) const
// function_index=39
void AA_exclass2_get_name_bufferify(const AA_exclass2 * self, char * SH_F_rv, int LSH_F_rv)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name_bufferify
const std::string & rv = selfobj->getName();
asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end class.ExClass2.method.get_name_bufferify
}

// const string & getName2()
// function_index=21
const char * AA_exclass2_get_name2(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name2
const std::string & rv = selfobj->getName2();
return rv.c_str();
// splicer end class.ExClass2.method.get_name2
}

// void getName2(string & SH_F_rv+intent(out)+len(LSH_F_rv))
// function_index=40
void AA_exclass2_get_name2_bufferify(AA_exclass2 * self, char * SH_F_rv, int LSH_F_rv)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name2_bufferify
const std::string & rv = selfobj->getName2();
asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end class.ExClass2.method.get_name2_bufferify
}

// string & getName3() const
// function_index=22
char * AA_exclass2_get_name3(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name3
std::string & rv = selfobj->getName3();
return rv.c_str();
// splicer end class.ExClass2.method.get_name3
}

// void getName3(string & SH_F_rv+intent(out)+len(LSH_F_rv)) const
// function_index=41
void AA_exclass2_get_name3_bufferify(const AA_exclass2 * self, char * SH_F_rv, int LSH_F_rv)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_name3_bufferify
std::string & rv = selfobj->getName3();
asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end class.ExClass2.method.get_name3_bufferify
}

// string & getName4()
// function_index=23
char * AA_exclass2_get_name4(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name4
std::string & rv = selfobj->getName4();
return rv.c_str();
// splicer end class.ExClass2.method.get_name4
}

// void getName4(string & SH_F_rv+intent(out)+len(LSH_F_rv))
// function_index=42
void AA_exclass2_get_name4_bufferify(AA_exclass2 * self, char * SH_F_rv, int LSH_F_rv)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name4_bufferify
std::string & rv = selfobj->getName4();
asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end class.ExClass2.method.get_name4_bufferify
}

// const int GetNameLength()
// function_index=24
/**
 * \brief helper function for Fortran
 *
 */
const int AA_exclass2_get_name_length(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_name_length
return selfobj->getName().length();
// splicer end class.ExClass2.method.get_name_length
}

// ExClass1 * get_class1(const ExClass1 * in+intent(in)+value)
// function_index=25
AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, const AA_exclass1 * in)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_class1
ExClass1 * rv = selfobj->get_class1(static_cast<const ExClass1 *>(static_cast<const void *>(in)));
return static_cast<AA_exclass1 *>(static_cast<void *>(rv));
// splicer end class.ExClass2.method.get_class1
}

// void * declare(TypeID type+intent(in)+value)
// function_index=31
void AA_exclass2_declare_0(AA_exclass2 * self, int type)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.declare_0
selfobj->declare(getTypeID(type));
return;
// splicer end class.ExClass2.method.declare_0
}

// void * declare(TypeID type+intent(in)+value, SidreLength len+default(1)+intent(in)+value)
// function_index=26
void AA_exclass2_declare_1(AA_exclass2 * self, int type, ATK_SidreLength len)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.declare_1
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.ExClass2.method.declare_1
}

// void destroyall()
// function_index=27
void AA_exclass2_destroyall(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.destroyall
selfobj->destroyall();
return;
// splicer end class.ExClass2.method.destroyall
}

// TypeID getTypeID() const
// function_index=28
int AA_exclass2_get_type_id(const AA_exclass2 * self)
{
const ExClass2 *selfobj = static_cast<const ExClass2 *>(static_cast<const void *>(self));
// splicer begin class.ExClass2.method.get_type_id
TypeID rv = selfobj->getTypeID();
return static_cast<int>(rv);
// splicer end class.ExClass2.method.get_type_id
}

// void setValue(int value+intent(in)+value)
// function_index=32
void AA_exclass2_set_value_int(AA_exclass2 * self, int value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_int
selfobj->setValue<int>(value);
return;
// splicer end class.ExClass2.method.set_value_int
}

// void setValue(long value+intent(in)+value)
// function_index=33
void AA_exclass2_set_value_long(AA_exclass2 * self, long value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_long
selfobj->setValue<long>(value);
return;
// splicer end class.ExClass2.method.set_value_long
}

// void setValue(float value+intent(in)+value)
// function_index=34
void AA_exclass2_set_value_float(AA_exclass2 * self, float value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_float
selfobj->setValue<float>(value);
return;
// splicer end class.ExClass2.method.set_value_float
}

// void setValue(double value+intent(in)+value)
// function_index=35
void AA_exclass2_set_value_double(AA_exclass2 * self, double value)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.set_value_double
selfobj->setValue<double>(value);
return;
// splicer end class.ExClass2.method.set_value_double
}

// int getValue()
// function_index=36
int AA_exclass2_get_value_int(AA_exclass2 * self)
{
ExClass2 *selfobj = static_cast<ExClass2 *>(static_cast<void *>(self));
// splicer begin class.ExClass2.method.get_value_int
int rv = selfobj->getValue<int>();
return rv;
// splicer end class.ExClass2.method.get_value_int
}

// double getValue()
// function_index=37
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
